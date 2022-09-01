from PyAutoDock.utils import calc_ddd_Mehler_Solmajer


class GridMap:
    def __init__(self,num_receptor_maps=None,*args,**kwargs):
        self.num_receptor_maps = num_receptor_maps if num_receptor_maps else 0
        self.atomtype = None
        self.is_covalent = None
        self.energy_max = 0.0
        self.energy_min = 0.0
        self.energy = 0.0
        self.vol = 0.0
        self.sol = 0.0
        self.Rij = 0.0
        self.epsij = 0.0
        self.hbond = 0
        self.Rij_hb = 0.0
        self.epsij_hb = 0.0
        self.num_cmaps = num_receptor_maps      # used for number of valid maps
        self.catomtypes = ['X' for i in range(num_receptor_maps)]
        self.cA = [0.0 for i in range(num_receptor_maps)]
        self.cB = [0.0 for i in range(num_receptor_maps)]
        self.nbp_r = [0.0 for i in range(num_receptor_maps)]
        self.nbp_eps = [0.0 for i in range(num_receptor_maps)]
        self.xA = [0.0 for i in range(num_receptor_maps)]
        self.xB = [0.0 for i in range(num_receptor_maps)]
        self.hbonder = [0 for i in range(num_receptor_maps)]    # int


def setup_atomtypes_interaction_maps(
        referatomtypes,referpar,againstatomtypes=None,againstpar=None,loglevel=None
    ):
        """Setup atomtypes interaction maps

        Maps can be calculated for self or cross interactions, for example:
        If we have atomtypes:
            referatomtypes = ['A','C','N']
            againstatomtypes = ['O', 'S', 'P']

        For cross-interaction maps:
            ['A-O', 'A-S', 'A-P', 'C-O', 'C-S', 'C-P', 'N-O', 'N-S', 'N-P']

        For self-interaction maps:
            ['A-A', 'A-N', 'C-C', 'C-N', 'N-N']
            Or: ['O-O', 'O-S', 'O-P', 'S-S', 'S-P', 'P-P']

        Sequences are correspondent with the inputs
        Detail info will be printout if loglevel<=10
        """
        if againstatomtypes:
            num_maps = len(againstatomtypes)
            fitatomtypes = [againstatomtypes for i in range(num_maps)]
            fitpar = [againstpar for i in range(num_maps)]
        else:
            # means self
            num_maps = len(referatomtypes)
            fitatomtypes = [referatomtypes[i:] for i in range(num_maps)]
            fitpar = [referpar[i:] for i in range(num_maps)]

        MAPS = []
        for i in range(len(referatomtypes)):
            m = GridMap(num_receptor_maps=num_maps)
            m.atomtype = referatomtypes[i]
            par = referpar[i]
            m.sol = par.sol
            m.vol = par.vol
            m.Rij = par.Rij
            m.epsij = par.epsij
            m.hbond = par.hbond
            m.Rij_hb = par.Rij_hb
            m.epsij_hb = par.epsij_hb
            m.num_cmaps = len(fitatomtypes[i])
            for j in range(len(fitatomtypes[i])):
                p = fitpar[i][j]
                m.catomtypes[j] = fitatomtypes[i][j]
                m.nbp_r[j] = (m.Rij + p.Rij) / 2.0
                m.nbp_eps[j] = pow(m.epsij*p.epsij, 0.5)
                m.xA[j] = 12
                m.xB[j] = 6
                if m.hbond > 2 and p.hbond in [1,2]:
                    m.xB[j] = 10
                    m.hbonder[j] = True
                    m.nbp_r[j] = m.Rij_hb
                    m.nbp_eps[j] = m.epsij_hb
                elif m.hbond in [1,2] and p.hbond > 2:
                    m.xB[j] = 10
                    m.hbonder[j] = True
                    m.nbp_r[j] = p.Rij_hb
                    m.nbp_eps[j] = p.epsij_hb
            MAPS.append(m)

        if loglevel and loglevel <= 10:
            print()
            for i,a in enumerate(referatomtypes):
                m = MAPS[i]
                print('\nFor ligand type: {:}'.format(a))
                print(
                    'hbond={:}, sol={:.8f}, vol={:.6f}, Rij={:.3f}, '
                    'epsij={:.4f}, Rij_hb={:.4f}, epsij_hb={:.4f}'.format(
                        m.hbond, m.sol, m.vol, m.Rij, m.epsij,
                        m.Rij_hb, m.epsij_hb
                    )
                )
                for j in range(num_maps):
                    print(
                        '>> receptor_atomtype={:}, hbonder={:}, xA={:}, xB={:}, '
                        'nbp_r={:.4f}, nbp_eps={:.4f}'.format(m.catomtypes[j],
                            m.hbonder[j], m.xA[j], m.xB[j], m.nbp_r[j], m.nbp_eps[j]
                        )
                    )
            print()

        return MAPS


def setup_eps_table(num_points, gridsize=None,is_distance=True,loglevel=None):
    if not gridsize: gridsize = 100
    if not is_distance:
        return [1.0 for i in range(num_points*gridsize)]

    EPSTable = [calc_ddd_Mehler_Solmajer(i/gridsize) for i in range(num_points)]
    EPSTable[0] = 1.0   # important
    if loglevel and loglevel <= 10:
        print('\nUsing *distance-dependent* dielectric function of Mehler and Solmajer, Prot.Eng.4, 903-910.')
        print('   d     Dielectric')
        for i in range(0,min(310,num_points),10):
            print('{:>5.1f} {:>9.2f}'.format(i/gridsize,EPSTable[i]))
        print('\n')
    return EPSTable


def setup_vdW_hbond_table(
        num_points,referatomtypes,referpar,againstatomtypes=None,againstpar=None,
        gridsize=None,num_smooth=None,eclamp=None,loglevel=None
    ):
    if not gridsize: gridsize = 100
    if not eclamp: eclamp = 100000

    MAPS = setup_atomtypes_interaction_maps(referatomtypes,referpar,againstatomtypes,againstpar,loglevel)
    if not MAPS: return []

    bigest = max([t.num_cmaps for t in MAPS])
    EvdWHBondTable = [
        [[0.0 for k in range(len(MAPS))] for j in range(bigest)]
        for i in range(num_points)
    ]

    for i in range(len(MAPS)):       # index i is super important, be careful
        if MAPS[i].is_covalent:
            print('Note: covalent map: any internal non-bonded parameters will be ignored')
            print('Note: covalent map: ligand_type={:}'.format(MAPS[i].atomtype))
            continue
        m = MAPS[i]
        for j in range(MAPS[i].num_cmaps):
            tmp = m.nbp_eps[j] / (m.xA[j]-m.xB[j])
            m.cA[j] = tmp * pow(m.nbp_r[j],m.xA[j]) * m.xB[j]
            m.cB[j] = tmp * pow(m.nbp_r[j],m.xB[j]) * m.xA[j]
            for k in range(1,num_points):   # care in here, starting from 1
                r = k / gridsize
                rA = pow(r,m.xA[j])
                rB = pow(r,m.xB[j])
                EvdWHBondTable[k][j][i] = min(eclamp,m.cA[j]/rA-m.cB[j]/rB)
            EvdWHBondTable[0][j][i] = eclamp
            EvdWHBondTable[num_points-1][j][i] = 0.0

    if loglevel and loglevel <= 10:
        txt = 'Before Smooth: ' if num_smooth and num_smooth > 0 else ''
        for i in range(len(MAPS)):       # index i is super important, be careful
            print('\n{:}Pairwise interactions:'.format(txt))
            out = '>> {:>3}-'.format(MAPS[i].atomtype)
            cn = MAPS[i].num_cmaps
            for j in range(cn):
                tmpa = 'cA={:.2f} / {:}'.format(MAPS[i].cA[j],MAPS[i].xA[j])
                tmpb = 'cB={:.2f} / {:}'.format(MAPS[i].cB[j],MAPS[i].xB[j])
                print(out+'{:<3}   {:30}   {:}'.format(MAPS[i].catomtypes[j],tmpa,tmpb))

            print('>> Atomtype: {:}'.format(MAPS[i].atomtype))
            out = '  r  ' + ' '.join(['{:^9}'.format(v) for v in MAPS[i].catomtypes[:cn]])
            out += '\n' + ' ___ ' + ' '.join(['_'*9 for v in range(cn)])
            print(out)
            for t in range(0,min(num_points,510),10):
                out = '{:4.1f} '.format(t/gridsize)
                for j in range(cn):
                    out += '{:9.2f} '.format(EvdWHBondTable[t][j][i])
                print(out)
            print()

    if num_smooth and num_smooth > 0:
        # do smooth
        for i in range(len(MAPS)):       # index i is super important, be careful
            etmp = [0.0 for t in range(num_points)]
            for j in range(MAPS[i].num_cmaps):
                for k in range(num_points):
                    etmp[k] = eclamp
                    for s in range(max(0,k-num_smooth), min(num_points,k+num_smooth+1)):
                        etmp[k] = min(etmp[k], EvdWHBondTable[s][j][i])
                for k in range(num_points):
                    EvdWHBondTable[k][j][i] = etmp[k]

        if loglevel and loglevel <= 10:
            for i in range(len(MAPS)):       # index i is super important, be careful
                print('\nAfter Smooth: Pairwise interactions:')
                print('Lowest pairwise interaction energy within {:.3f} Angstrom'.format(num_smooth*2.0/gridsize))
                out = '>> {:>3}-'.format(MAPS[i].atomtype)
                cn = MAPS[i].num_cmaps
                for j in range(cn):
                    tmpa = 'cA={:.2f} / {:}'.format(MAPS[i].cA[j],MAPS[i].xA[j])
                    tmpb = 'cB={:.2f} / {:}'.format(MAPS[i].cB[j],MAPS[i].xB[j])
                    print(out+'{:<3}   {:30}   {:}'.format(MAPS[i].catomtypes[j],tmpa,tmpb))

                print('>> Atomtype: {:}'.format(MAPS[i].atomtype))
                out = '  r  ' + ' '.join(['{:^9}'.format(v) for v in MAPS[i].catomtypes[:cn]])
                out += '\n' + ' ___ ' + ' '.join(['_'*9 for v in range(cn)])
                print(out)
                for t in range(0,min(num_points,510),10):
                    out = '{:4.1f} '.format(t/gridsize)
                    for j in range(cn):
                        out += '{:9.2f} '.format(EvdWHBondTable[t][j][i])
                    print(out)
                print()

    return EvdWHBondTable


