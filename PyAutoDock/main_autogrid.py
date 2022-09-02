from PyAutoDock.logger import mylogger
from PyAutoDock.read_gpf import ReadGPF
from PyAutoDock.read_pdbqt import ReadPDBQT
from PyAutoDock.setup_par_library import SetupParLibrary
from PyAutoDock.utils import file_gen_new, calc_ddd_Mehler_Solmajer

import math

logger = mylogger()
if logger.level == 0: logger.level = 20     # NOTSET


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
        self.cA = [0.0 for i in range(num_receptor_maps)]
        self.cB = [0.0 for i in range(num_receptor_maps)]
        self.nbp_r = [0.0 for i in range(num_receptor_maps)]
        self.nbp_eps = [0.0 for i in range(num_receptor_maps)]
        self.xA = [0.0 for i in range(num_receptor_maps)]
        self.xB = [0.0 for i in range(num_receptor_maps)]
        self.hbonder = [0 for i in range(num_receptor_maps)]    # int


class SetupGridMaps:
    """setup receptor/ligand interactions grid maps

    Args:
        library_filename(str)       : library file
        ligand_gpf_filename(str)    : grid parameter file
        receptor_mol_filename(str)  : file for receptor, it takes highest precedence

        gpf (collections.OrderedDict):
        EnergyMaps ( List(x(y(z(num_maps)))) ):
    """
    def __init__(
            self,library_filename,ligand_gpf_filename,receptor_mol_filename=None,
            nbc=None,eclamp=None,tol=None,gsize=None,eintcutoff=None,
            *args,**kwargs
        ):
        self.library_filename = library_filename
        self.receptor_mol_filename = receptor_mol_filename if receptor_mol_filename else None
        self.ligand_gpf_filename = ligand_gpf_filename
        self.nbc = nbc if nbc else 8
        self.eclamp = eclamp if eclamp else 100000
        self.tol = tol if tol else 0.00000001
        self.gsize = gsize if gsize else 100
        self.eintcutoff = eintcutoff if eintcutoff else 2047
        self.gpf = None
        self.EnergyMaps = None
        self.gen(tol=self.tol)
    
    def setup_library(self,library_filename=None):
        if not library_filename: library_filename = self.library_filename
        return SetupParLibrary(library_filename)

    def gen(self,tol=None):
        tol = tol if tol else 0.00000001
        LIB = SetupParLibrary(self.library_filename)
        if not len(LIB.atoms): return
        GPF = ReadGPF(self.ligand_gpf_filename,loglevel=100)

        if self.receptor_mol_filename:
            MOL = ReadPDBQT(self.receptor_mol_filename)
            if not len(MOL.atoms): return
        elif GPF.gpf['receptor']:
            MOL = ReadPDBQT(self.GPF.gpf['receptor'])
            if not len(MOL.atoms): return
        else:
            logger.critical('required: PDBQT file: receptor')
            return

        atomtypes = set([a[0] for a in MOL.atoms])
        for a in atomtypes:
            if not LIB.get_atom_par(a):
                logger.critical('not defined: receptor atomtype: {:}'.format(a))

        #TODO guessing
        if not GPF.gpf['ligand_types']:
            logger.critical('required: gpf({:}): ligand_types'.format(self.ligand_gpf_filename))

        if GPF.gpf['npts']:
            num_grids = 1
            for i in GPF.gpf['npts']: num_grids = num_grids * (i+1)
        else:
            # guess npts from receptor file
            xn = int((MOL.xmax - MOL.xmin + 0.6) / GPF.gpf['spacing'])
            xn = xn if xn%2 == 0 else xn+1
            yn = int((MOL.ymax - MOL.ymin + 0.6) / GPF.gpf['spacing'])
            yn = yn if yn%2 == 0 else yn+1
            zn = int((MOL.zmax - MOL.zmin + 0.6) / GPF.gpf['spacing'])
            zn = zn if zn%2 == 0 else zn+1
            GPF.gpf['npts'] = [xn, yn, zn]
            num_grids = (xn+1) * (yn+1) * (zn+1)
        logger.info('number of total grids: {:}'.format(num_grids))

        if GPF.gpf['gridcenter']:
            # based on new gridcenter, recalculate box coordinates
            MOL.translate(*GPF.gpf['gridcenter'])
        else:
            GPF.gpf['gridcenter'] = MOL.center

        if GPF.gpf['receptor_types']:
            defined = set(GPF.gpf['receptor_types'])
            if MOL.atomset != defined:
                left = defined.difference(MOL.atomset)
                more = MOL.atomset.difference(defined)
                if left:
                    logger.warning('not defined: receptor_types: {:}'.format(left))
                if more:
                    logger.critical('not included: receptor_types: {:}'.format(more))
        else:
            GPF.gpf['receptor_types'] = list(MOL.atomset)
            GPF.gpf['map'] = []     # reset
            for a in GPF.gpf['receptor_types']:
                GPF.gpf['map'].append(file_gen_new('receptor.{:}.map'.format(a)))
            logger.info('receptor_types: {:}'.format(GPF.gpf['receptor_types']))
            logger.info('map: {:}'.format(GPF.gpf['map']))
        MOLATOMSINDEX = [GPF.gpf['receptor_types'].index(i[0]) for i in MOL.atoms]

        # check receptor_types and its maps
        if len(GPF.gpf['ligand_types']) != len(GPF.gpf['map']):
            logger.critical('not equal: number of entries in: ligand_types & map')

        # epsilon distance lookup table
        # we do not need build the very big table, slightly bigger than molecule box is OK
        dx = MOL.xmax - MOL.xmin
        dy = MOL.ymax - MOL.ymin
        dz = MOL.zmax - MOL.zmin
        rmax = pow(dx*dx+dy*dy+dz*dz,0.5)
        points = int(rmax+0.6) * self.gsize
        if GPF.gpf['dielectric'] > 0:
            #TODO factor not known
            logger.info('Using a *constant* dielectric: {:}'.format(GPF.gpf['dielectric']))
            EPSTable = [1.0 for i in range(points)]
        else:
            EPSTable = [calc_ddd_Mehler_Solmajer(i/self.gsize) for i in range(points)]
            EPSTable[0] = 1.0   # important
            logger.info('Using *distance-dependent* dielectric function of Mehler and Solmajer, Prot.Eng.4, 903-910.')
            if logger.level <= 10:
                logger.debug('printout distance table\n')
                print('   d     Dielectric')
                for i in range(0,min(310,points),10):
                    print('{:>5.1f} {:>9.2f}'.format(i/self.gsize,EPSTable[i]))
            for i in range(len(EPSTable)):
                EPSTable[i] = 332 / EPSTable[i]     # waiting on further check 332.06363

        prec = []
        for a in GPF.gpf['receptor_types']:
            par = LIB.get_atom_par(a)
            if par:
                prec.append(par)
            else:
                logger.critical('not defined: receptor_types: {:}'.format(a))

        plig = []
        for a in GPF.gpf['ligand_types']:
            par = LIB.get_atom_par(a)
            if par:
                plig.append(par)
            else:
                logger.critical('not defined: ligand_types: {:}'.format(a))

        # calculate receptor atoms distance map
        #TODO awaits negotiation, takes too much memory

        MAPS = []
        for i in range(len(GPF.gpf['ligand_types'])):
            m = GridMap(num_receptor_maps=len(GPF.gpf['receptor_types']))
            m.atomtype = GPF.gpf['ligand_types'][i]
            par = plig[i]
            m.sol = par.sol
            m.vol = par.vol
            m.Rij = par.Rij
            m.epsij = par.epsij
            m.hbond = par.hbond
            m.Rij_hb = par.Rij_hb
            m.epsij_hb = par.epsij_hb
            for j in range(len(GPF.gpf['receptor_types'])):
                p = prec[j]
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
        # additional electric and dsolvation
        m = GridMap(num_receptor_maps=len(GPF.gpf['receptor_types']))
        m.atomtype = 'electric'
        MAPS.append(m)
        m = GridMap(num_receptor_maps=len(GPF.gpf['receptor_types']))
        m.atomtype = 'dsolvation'
        MAPS.append(m)

        if logger.level <= 10:
            print()
            for i,a in enumerate(GPF.gpf['ligand_types']):
                m = MAPS[i]
                print('\nFor ligand type: {:}'.format(a))
                print(
                    'hbond={:}, sol={:.8f}, vol={:.6f}, Rij={:.3f}, '
                    'epsij={:.4f}, Rij_hb={:.4f}, epsij_hb={:.4f}'.format(
                        m.hbond, m.sol, m.vol, m.Rij, m.epsij,
                        m.Rij_hb, m.epsij_hb
                    )
                )
                for j,r in enumerate(GPF.gpf['receptor_types']):
                    print(
                        '>> receptor_atomtype={:}, hbonder={:}, xA={:}, xB={:}, '
                        'nbp_r={:.4f}, nbp_eps={:.4f}'.format(
                            r, m.hbonder[j], m.xA[j], m.xB[j], m.nbp_r[j], m.nbp_eps[j]
                        )
                    )

        EvdWHBondTable = [
            [
                [0.0 for i in range(len(GPF.gpf['ligand_types']))]
                for j in range(len(GPF.gpf['receptor_types']))
            ]
            for k in range(points)
        ]
        for i in range(len(GPF.gpf['ligand_types'])):       # index i is super important, be careful
            if MAPS[i].is_covalent:
                logger.info('For covalent map: any internal non-bonded parameters will be ignored')
                logger.info('covalent map: ligand_type={:}'.format(GPF.gpf['ligand_types'][i]))
                continue
            m = MAPS[i]
            for j in range(len(GPF.gpf['receptor_types'])):
                tmp = m.nbp_eps[j] / (m.xA[j]-m.xB[j])
                m.cA[j] = tmp * pow(m.nbp_r[j],m.xA[j]) * m.xB[j]
                m.cB[j] = tmp * pow(m.nbp_r[j],m.xB[j]) * m.xA[j]
                for k in range(1,points):   # care in here, starting from 1
                    r = k / self.gsize
                    rA = pow(r,m.xA[j])
                    rB = pow(r,m.xB[j])
                    EvdWHBondTable[k][j][i] = min(self.eclamp,m.cA[j]/rA-m.cB[j]/rB)
                EvdWHBondTable[0][j][i] = self.eclamp
                EvdWHBondTable[points-1][j][i] = 0.0

            # smooth
            n = int(GPF.gpf['smooth'] * self.gsize / 2.0 + 0.6)

            if logger.level <= 10 and n > 0:
                tmp = 'Before Smooth: ' if n > 0 else ''
                print('\n{:}Pairwise interactions:'.format(tmp))
                out = '>> {:>3}-'.format(GPF.gpf['ligand_types'][i])
                for j,v in enumerate(GPF.gpf['receptor_types']):
                    tmpa = 'cA={:.2f} / {:}'.format(m.cA[j],m.xA[j])
                    tmpb = 'cB={:.2f} / {:}'.format(m.cB[j],m.xB[j])
                    print(out+'{:<3}   {:30}   {:}'.format(v,tmpa,tmpb))

                print('\n\n{:}Lowest pairwise interaction energy within {:} Angstrom'.format(tmp,GPF.gpf['smooth']))
                print('>> Atomtype: {:}'.format(GPF.gpf['ligand_types'][i]))
                out = '  r  ' + ' '.join(['{:^9}'.format(v) for v in GPF.gpf['receptor_types']])
                out += '\n' + ' ___ ' + ' '.join(['_'*9 for v in GPF.gpf['receptor_types']])
                print(out)
                #print('  r      A         C         H         HD        N        NA        OA        SA')
                #print(' ___ _________ _________ _________ _________ _________ _________ _________ _________')
                for t in range(0,min(points,510),10):
                    out = '{:4.1f} '.format(t/self.gsize)
                    for j in range(len(GPF.gpf['receptor_types'])):
                        out += '{:9.2f} '.format(EvdWHBondTable[t][j][i])
                    print(out)
                print()

            if n > 0:
                # do smooth
                etmp = [0.0 for t in range(points)]
                for j in range(len(GPF.gpf['receptor_types'])):
                    for k in range(points):
                        etmp[k] = self.eclamp
                        for s in range(max(0,k-n), min(points,k+n+1)):
                            etmp[k] = min(etmp[k], EvdWHBondTable[s][j][i])
                    for k in range(points):
                        EvdWHBondTable[k][j][i] = etmp[k]

            if logger.level <= 10 and n > 0:
                tmp = 'After Smooth: ' if n > 0 else ''
                print('\n\n{:}Pairwise interactions:'.format(tmp))
                out = '>> {:>3}-'.format(GPF.gpf['ligand_types'][i])
                for j,v in enumerate(GPF.gpf['receptor_types']):
                    tmpa = 'cA={:.2f} / {:}'.format(m.cA[j],m.xA[j])
                    tmpb = 'cB={:.2f} / {:}'.format(m.cB[j],m.xB[j])
                    print(out+'{:<3}   {:30}   {:}'.format(v,tmpa,tmpb))

                print('\n\n{:}Lowest pairwise interaction energy within {:} Angstrom'.format(tmp,GPF.gpf['smooth']))
                print('>> Atomtype: {:}'.format(GPF.gpf['ligand_types'][i]))
                out = '  r  ' + ' '.join(['{:^9}'.format(v) for v in GPF.gpf['receptor_types']])
                out += '\n' + ' ___ ' + ' '.join(['_'*9 for v in GPF.gpf['receptor_types']])
                print(out)
                #print('  r      A         C         H         HD        N        NA        OA        SA')
                #print(' ___ _________ _________ _________ _________ _________ _________ _________ _________')
                for t in range(0,min(points,510),10):
                    out = '{:4.1f} '.format(t/self.gsize)
                    for j in range(len(GPF.gpf['receptor_types'])):
                        out += '{:9.2f} '.format(EvdWHBondTable[t][j][i])
                    print(out)
                print()

        # setup solvent maps, sigma = 3.6 Angstrom
        ESolTable = [0.0 for i in range(points)]
        for t in range(1,points):       # starts from 1
            ESolTable[t] = LIB.e.FE_coeff_desolv * math.exp(-t*t/10000.0/(2*3.6*3.6))

        #TODO parameter inside GPF
        Disorder_h = False

        # for hbond maps
        # format:
        #      vector     rexp   normvector
        #   [(rx,ry,rz), ... ]                  # for ha == 4 or `left`
        #   [((rx,ry,rz), rexp), ... ]          # for ha == 2
        #   [((rx,ry,rz), (nx,ny,nz)), ... ]    # for ha == 5
        #
        # important:
        #   a) values are in different formats
        #   b) corresponding to the order
        HBMAPS = []
        HBONDTYPES = []
        for ia,a in enumerate(MOL.atoms):
            ha = LIB.get_atom_hbond(a[0])
            HBONDTYPES.append(ha)
            if ha == 2:     # D1
                boset = False
                for jb,b in enumerate(MOL.atoms):
                    if ia == jb: continue
                    v = [a[t]-b[t] for t in range(2,5)]
                    dd = sum([t*t for t in v])
                    if dd < 1.90:
                        boset = True
                        btype = b[0].upper()
                        if btype == 'OA' or btype == 'SA':
                            rexp = 4
                        else:
                            rexp = 2
                        dd = max(dd, tol)
                        dinv = 1.0 / pow(dd,0.5)
                        v = tuple([t*dinv for t in v])
                        HBMAPS.append((v,rexp))
                        # Question: only find one of them??
                        break
                if not boset:
                    logger.warning(f'Hydrogen bond: Receptor atom index: {ia+1}')
                    HBMAPS.append(((0.0,0.0,0.0),0))
            elif ha == 4:   # A1
                nbond = 0
                j1 = 0
                j2 = 0
                j3 = 0
                # to follow the original codes
                bfrom = max(ia-20,0)
                bto = min(ia+20, len(MOL.atoms)-1)
                for jb in range(bfrom,bto):
                    b = MOL.atoms[jb]
                    if ia == jb: continue
                    v = [a[t]-b[t] for t in range(2,5)]
                    dd = sum([t*t for t in v])
                    btype = b[0].upper()
                    if (dd < 2.89 and (btype != 'HD' and btype != 'H')) or \
                        (dd < 1.69 and (btype == 'HD' or btype == 'H')):
                        if nbond == 0:
                            nbond = 1
                            j1 = jb
                        elif nbond == 1:
                            nbond = 2
                            j2 = jb
                        elif nbond == 2:
                            nbond = 3
                            j3 = jb
                        else:
                            logger.warning(
                                'Found an N-bonding atom with four bonded atoms, '
                                'atom serial number {:}'.format(ia+1)
                            )
                v = (0.0, 0.0, 0.0)
                if nbond == 0:
                    logger.warning(f'Nitrogen atom found with no bonded atoms, atom serial number: {ia+1}')
                elif nbond == 1:    # one bond: Azide Nitrogen :N=C-X
                    v = self.normVpp(MOL.atoms[j1][2:5], a[2:5])
                elif nbond == 2:    # two bonds: X1-N=X2
                    #TODO test
                    tmpv = [(MOL.atoms[j1][t]+MOL.atoms[j2][t])/2.0 for t in range(2,5)]
                    v = self.normVpp(tmpv, a[2:5])
                elif nbond == 3:    # three bonds: X1,X2,X3
                    #TODO test
                    tmpv = [(MOL.atoms[j1][t]+MOL.atoms[j2][t]+MOL.atoms[j3][t])/3.0 for t in range(2,5)]
                    v = self.normVpp(tmpv, a[2:5])
                else:
                    logger.critical('ha=4: How can it be??')
                HBMAPS.append(tuple(v))
            elif ha == 5:   # A2
                nbond = 0
                j1 = 0
                j2 = 0
                for jb,b in enumerate(MOL.atoms):
                    if ia == jb: continue
                    v = [a[t]-b[t] for t in range(2,5)]
                    dd = max(sum([t*t for t in v]), tol)
                    btype = b[0].upper()
                    if (dd < 3.61 and (btype != 'HD' and btype != 'H')) or \
                        (dd < 1.69 and (btype == 'HD' or btype == 'H')):
                        if nbond == 0:
                            nbond = 1
                            j1 = jb
                        elif nbond == 1:
                            nbond = 2
                            j2 = jb
                        else:
                            logger.warning(
                                'Found an H-bonding atom with three bonded atoms, '
                                'atom serial number {:}'.format(ia+1)
                            )
                if nbond == 0:
                    logger.warning(f'Oxygen atom found with no bonded atoms, atom serial number: {ia+1}')
                    HBMAPS.append(((0.0,0.0,0.0),(0.0,0.0,0.0)))
                elif nbond == 1:        # one bond: carbonyl oxygen: O=C-X
                    v1 = [a[t]-MOL.atoms[j1][t] for t in range(2,5)]
                    dd1 = max(sum([t*t for t in v1]), tol)
                    d1inv = 1.0 / pow(dd1,0.5)
                    v1 = tuple([t*d1inv for t in v1])
                    v2 = None
                    # try to find second
                    for j2,c in enumerate(MOL.atoms):
                        if j2 == ia or j2 == j1: continue
                        vd = [c[t]-MOL.atoms[j1][t] for t in range(2,5)]
                        dd2 = max(sum([t*t for t in vd]), tol)
                        btype = c[0].upper()
                        if (dd2 < 2.89 and btype != 'HD') or (dd2 < 1.69 and btype == 'HD'):
                            d2inv = 1.0 / pow(dd2,0.5)
                            vd = [t*d2inv for t in vd]
                            # C=O cross C-X to get long pair plane norm vector
                            x = v1[1]*vd[2] - v1[2]*vd[1]
                            y = -v1[0]*vd[2] + v1[2]*vd[0]
                            z = v1[0]*vd[1] - v1[1]*vd[0]
                            d = pow(max(x*x+y*y+z*z, tol), 0.5)
                            v2 = (x/d, y/d, z/d)
                            #TODO why not break?
                            #break
                    if not v2:
                        logger.warning(f'Oxygen atom found only one bonded atom, atom serial number: {ia+1}')
                        v2 = (0.0, 0.0, 0.0)
                    HBMAPS.append((v1,v2))
                elif nbond == 2:
                    btype1 = MOL.atoms[j1][0].upper()
                    btype2 = MOL.atoms[j2][0].upper()
                    if Disorder_h and (btype1 == 'HD' or btype2 == 'HD') and (btype1 != btype2):
                        #TODO test
                        ndx = None
                        if btype1 == 'A' or btype1 == 'C': ndx = j1
                        if btype1 == 'A' or btype2 == 'C': ndx = j2
                        if not ndx:
                            ndx = 0
                            logger.warning('parameters error')
                        v = [a[t]-MOL.atoms[ndx][t] for t in range(2,5)]
                        dd = max(sum([t*t for t in v]), tol)
                        dinv = 1.0 / pow(dd,0.5)
                        v = tuple([t*dinv for t in v])
                        HBMAPS.append((v,(0.0,0.0,0.0)))
                    else:
                        v2 = [MOL.atoms[j2][t]-MOL.atoms[j1][t] for t in range(2,5)]
                        dd2 = max(sum([t*t for t in v2]), tol)
                        d2inv = 1.0 / pow(dd2,0.5)
                        v2 = tuple([t*d2inv for t in v2])

                        dot = sum([(a[t+2]-MOL.atoms[j1][t+2])*v2[t] for t in range(3)])
                        v1 = [a[t+2] - (dot*v2[t]+MOL.atoms[j1][t+2]) for t in range(3)]
                        dd1 = max(sum([t*t for t in v1]), tol)
                        d1inv = 1.0 / pow(dd1,0.5)
                        v1 = tuple([t*d1inv for t in v1])
                        HBMAPS.append((v1,v2))
                else:
                    logger.critical('ha=5: How can it be??')
            else:
                # important to append an empty tuple
                HBMAPS.append((0.0,0.0,0.0))

        totmaps = len(GPF.gpf['ligand_types'])+2
        logger.info('Beginning grid calculation')
        logger.info(
            'Calculating {:} grids over {:} elements, around {:} receptor atoms'.format(
                totmaps, num_grids, len(MOL.atoms)
            )
        )

        elecmap = totmaps - 1
        dsolmap = totmaps - 2
        estat = LIB.e.FE_coeff_estat
        hnpts = [t//2 for t in GPF.gpf['npts']]
        # format: [ x(y(z(num_maps(energy)))), ... ]
        EnergyMaps = [
            [
                [ [0.0 for t in range(totmaps)] for k in range(2*hnpts[2]+1) ]
                for j in range(2*hnpts[1]+1)
            ]
            for i in range(2*hnpts[0]+1)
        ]
        for iz in range(-hnpts[2], hnpts[2]+1):
            c = [0.0, 0.0, iz * GPF.gpf['spacing']]
            zndx = iz + hnpts[2]
            for iy in range(-hnpts[1], hnpts[1]+1):
                c[1] = iy * GPF.gpf['spacing']
                yndx = iy + hnpts[1]
                for ix in range(-hnpts[0],hnpts[0]+1):
                    c[0] = ix * GPF.gpf['spacing']
                    xndx = ix + hnpts[0]

                    Emap = EnergyMaps[xndx][yndx][zndx]

                    hbondmin = [999999.0 for t in range(totmaps-2)]
                    hbondmax = [-999999.0 for t in range(totmaps-2)]
                    hbondflag = [False for t in range(totmaps-2)]

                    # find closest Hbond
                    rmin = 999999.0
                    closestH = 0
                    for i,a in enumerate(MOL.atoms):
                        ha = HBONDTYPES[i]
                        if ha == 1 or ha == 2:
                            d = [a[t+2]-c[t] for t in range(3)]
                            vv = sum([t*t for t in d])
                            if vv < rmin:
                                rmin = vv
                                closestH = i

                    for ia,a in enumerate(MOL.atoms):
                        d = [a[t+2]-c[t] for t in range(3)]
                        vv = max(sum([t*t for t in d]),tol)
                        r = pow(vv,0.5)
                        d = [t/r for t in d]
                        inv_rmax = 1.0 / max(r,0.5)

                        index_r = int(r*self.gsize)
                        index_n = min(int(r*self.gsize), 2047)         #NEINT

                        if GPF.gpf['dielectric'] > 0:
                            #TODO test
                            Emap[elecmap] += a[5] * inv_rmax * 1.0/40 * estat
                        else:
                            Emap[elecmap] += a[5] * inv_rmax * EPSTable[index_r] * estat

                        # NBC
                        if r > self.nbc: continue

                        if Disorder_h and a[0].upper() == 'HD': continue

                        ha = HBONDTYPES[ia]
                        racc = 1.0
                        rdon = 1.0
                        Hramp = 1.0
                        if ha == 2:
                            v,rexp = HBMAPS[ia]
                            costheta = 0.0 - sum([d[t]*v[t] for t in range(3)])
                            if costheta <= 0:
                                racc = 0.0
                            else:
                                c2 = costheta * costheta
                                if rexp == 2:
                                    racc = c2
                                elif rexp == 4:
                                    racc = c2 * c2
                                else:
                                    racc = costheta
                                if ia == closestH:
                                    Hramp = 1.0
                                else:
                                    vc = HBMAPS[closestH][0]
                                    costheta = sum([v[t]*vc[t] for t in range(3)])
                                    if costheta > 1.0:
                                        costheta = 1.0
                                    elif costheta < -1.0:
                                        costheta = -1.0
                                    theta = math.acos(costheta)
                                    Hramp = 0.5 - 0.5*math.cos(theta*120.0/90.0)
                        elif ha == 4:
                            v = HBMAPS[ia]
                            costheta = 0.0 - sum([d[t]*v[t] for t in range(3)])
                            if costheta <= 0.0:
                                rdon = 0.0
                            else:
                                rdon = costheta * costheta
                        elif ha == 5:
                            v, v2 = HBMAPS[ia]
                            costheta = 0.0 - sum([d[t]*v[t] for t in range(3)])
                            t0 = sum([d[t]*v2[t] for t in range(3)])
                            if t0 > 1.0:
                                t0 = 1.0
                            elif t0 < -1.0:
                                t0 = -1.0
                            t0 = math.pi/2.0 - math.acos(t0)
                            cx = d[1]*v2[2] - d[2]*v2[1]
                            cy = d[2]*v2[0] - d[0]*v2[2]
                            cz = d[0]*v2[1] - d[1]*v2[0]
                            rd2 = cx*cx + cy*cy + cz*cz
                            rd2 = max(rd2,tol)
                            ti = (cx*v[0]+cy*v[1]+cz*v[2]) / pow(rd2,0.5)
                            rdon = 0.0
                            if costheta >= 0.0:
                                if ti > 1.0:
                                    ti = 1.0
                                elif ti < -1.0:
                                    ti = -1.0
                                ti = math.acos(ti) - math.pi/2.0
                                if ti < 0.0:
                                    ti = -ti
                                rdon = (0.9+0.1*math.sin(ti+ti)) * math.cos(t0)
                            elif costheta >= -0.34202:
                                rdon = 562.25 * pow(0.116978-costheta*costheta, 3) * math.cos(t0)

                        par = LIB.get_atom_par(a[0])
                        qasp = GPF.gpf['qasp']
                        for idx in range(totmaps-2):
                            gmap = MAPS[idx]
                            if gmap.is_covalent: continue
                            atype = MOLATOMSINDEX[ia]
                            if gmap.hbond > 0:
                                rsph = EvdWHBondTable[index_n][atype][idx] / self.gsize
                                if rsph < 0.0:
                                    rsph = 0.0
                                elif rsph > 1.0:
                                    rsph = 1.0
                                if (gmap.hbond == 3 or gmap.hbond == 5) and (ha == 1 or ha == 2):
                                    Emap[idx] += EvdWHBondTable[index_n][atype][idx] * Hramp * (racc+(1.0-racc)*rsph)
                                elif gmap.hbond == 4 and (ha == 1 or ha == 2):
                                    ene = EvdWHBondTable[index_n][atype][idx] * (racc+(1.0-racc)*rsph)
                                    hbondmin[idx] = min(hbondmin[idx],ene)
                                    hbondmax[idx] = max(hbondmax[idx],ene)
                                    hbondflag[idx] = True
                                elif (gmap.hbond == 1 or gmap.hbond == 2) and ha > 2:
                                    ene = EvdWHBondTable[index_n][atype][idx] * (rdon+(1.0-rdon)*rsph)
                                    hbondmin[idx] = min(hbondmin[idx],ene)
                                    hbondmax[idx] = max(hbondmax[idx],ene)
                                    hbondflag[idx] = True
                                else:
                                    Emap[idx] += EvdWHBondTable[index_n][atype][idx]
                            else:
                                Emap[idx] += EvdWHBondTable[index_n][atype][idx]

                            Emap[idx] += gmap.sol * par.vol * ESolTable[index_r] + \
                                    (par.sol+qasp*abs(a[5]))*gmap.vol*ESolTable[index_r]

                        # exclusively for dsolvation
                        Emap[dsolmap] += qasp * par.vol * ESolTable[index_r]

                    # adjust maps of hydrogen-bonding atoms by adding largest and
                    # smallest interaction of all 'pair-wise' interactions with receptor atoms
                    for idx in range(totmaps-2):
                        if hbondflag[idx]:
                            Emap[idx] += hbondmin[idx]
                            Emap[idx] += hbondmax[idx]

                    for idx in range(totmaps):
                        MAPS[idx].energy_max = max(MAPS[idx].energy_max, MAPS[idx].energy)
                        MAPS[idx].energy_min = max(MAPS[idx].energy_min, MAPS[idx].energy)

        # get the things that we want
        self.gpf = GPF.gpf
        self.EnergyMaps = EnergyMaps

    def normVpp(self,p1,p2,tol=None):
        """calculate norm vector for points in direction p1->p2"""
        v = [p2[t]-p1[t] for t in range(3)]
        return self.normV(v,tol)

    def normV(self,v,tol=None):
        """normalize vector v on given tolerance"""
        tol = tol if tol else self.tol
        dd = max(sum([t*t for t in v]),tol)
        d = pow(dd,0.5)
        return [t/d for t in v]




