from PyAutoDock.logger import mylogger
from PyAutoDock.read_gpf import ReadGPF
from PyAutoDock.read_pdbqt import ReadPDBQT
from PyAutoDock.setup_par_library import SetupParLibrary
from PyAutoDock.utils import file_gen_new, calc_ddd_Mehler_Solmajer

import os
import sys

logger = mylogger()


class GridMap:
    def __init__(self,num_receptor_maps=None,*args,**kwargs):
        self.num_receptor_maps = num_receptor_maps if num_receptor_maps else 0
        self.atomtype = None
        self.is_covalent = None
        self.is_hbonder = None
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


class SetupMolMaps:
    def __init__(self,library_filename,ligand_gpf_filename,receptor_mol_filename=None,*args,**kwargs):
        self.library_filename = library_filename
        self.receptor_mol_filename = receptor_mol_filename if receptor_mol_filename else None
        self.ligand_gpf_filename = ligand_gpf_filename
        self.maps = []

    def gen(self):
        LIB = SetupParLibrary(self.library_filename)
        if not len(LIB.atoms): return
        GPF = ReadGPF(self.ligand_gpf_filename)

        if self.receptor_mol_filename:
            MOL = ReadPDBQT(self.receptor_mol_filename)
            if not len(MOL.atoms): return
        elif GPF.gpf['receptor']:
            MOL = ReadPDBQT(self.GPF.gpf['receptor'])
            if not len(MOL.atoms): return
        else:
            logger.critical('required: PDBQT file: receptor')
            return
        
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
            MOL.centroid(*GPF.gpf['gridcenter'])
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
                    logger.warning('not included: receptor_types: {:}'.format(more))
        else:
            GPF.gpf['receptor_types'] = list(MOL.atomset)
            GPF.gpf['map'] = []     # reset
            for a in GPF.gpf['receptor_types']:
                GPF.gpf['map'].append(file_gen_new('receptor.{:}.map'.format(a)))
            logger.info('receptor_types: {:}'.format(GPF.gpf['receptor_types']))
            logger.info('map: {:}'.format(GPF.gpf['map']))

        # check receptor_types and its maps
        if len(GPF.gpf['receptor_types']) != (GPF.gpf['map']):
            logger.critical('not equal: number of entries in: receptor_types & map')

        # epsilon distance lookup table
        # we do not need build the very big table, slightly bigger than box is OK
        points = (max(GPF.gpf['npts'])+5) * 100
        if GPF.gpf['dielectric'] > 0:
            #TODO factor not known
            logger.info('Using a *constant* dielectric: {:}'.format(GPF.gpf['dielectric']))
            EPSTable = [1.0 for i in range(points)]
        else:
            EPSTable = [calc_ddd_Mehler_Solmajer(i/100.0) for i in range(points)]
            logger.info('Using *distance-dependent* dielectric function of Mehler and Solmajer, Prot.Eng.4, 903-910.')
            if logger.level < 10:
                logger.debug('printout distance table\n')
                print('  d   Dielectric')
                for i in range(0,min(300,points),10):
                    print('{:4.1f}{:9.2f}'.format(EPSTable[i]))
                print()
            for i in range(len(EPSTable)):
                EPSTable[i] = 332.06363 / EPSTable[i]

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
            m.sol = par.solpar
            m.vol = par.vol
            m.Rij = par.Rij
            m.epsij = par.epsij
            m.hbond = par.hbond
            m.Rij_hb = par.Rij_hb
            m.epsij_hb = par.epsij_hb
            if m.hbond > 0:
                m.is_hbonder = True
            for j in range(len(GPF.gpf['receptor_types'])):
                p = prec[j]
                m.nbp_r[j] = (m.Rij + p.Rij) / 2.0
                m.nbp_eps = pow(m.epsij*p.epsij, 0.5)
                m.xA[j] = 12
                m.xB[j] = 6
                if m.hbonder > 2 and p.hbond in [1,2]:
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

        if logger.level < 10:
            print('\n')
            for i,a in enumerate(GPF.gpf['ligand_types']):
                m = MAPS[i]
                print('\nFor ligand type: {:}'.format(a))
                print(
                    'is_hbonder={:}, hbond={:}, sol={:.8f}, vol={:.6f}, Rij={:.3f}, '
                    'epsij={:.3f}, Rij_hb={:.3f}, epsij_hb={:.3f}'.format(
                        m.is_hbonder, m.hbond, m.sol, m.vol, m.Rij, m.epsij,
                        m.Rij_hb, m.epsij_hb
                    )
                )
                for j,r in enumerate(GPF.gpf['receptor_types']):
                    print(
                        '>> receptor_atomtype={:}, hbonder={:}, xA={:}, xB={:}, '
                        'nbp_r={:.3f}, nbp_eps={:.3f}'.format(
                            r, m.hbonder[j], m.xA[j], m.xB[j], m.nbp_r[j], m.nbp_eps[j]
                        )
                    )
            print('\n')

        #TODO parameter inside GPF
        Disorder_h = False

        # for hbond maps
        # format:
        #        index         hbond         index      rexp    vector    disorder  normvector
        #   [(atom_i_index, atom_i_hbond, atom_j_index, rexp, (rx,ry,rz), disorder), ...]
        #   [(atom_i_index, atom_i_hbond, atom_j_index, rexp, (rx,ry,rz), disorder, (nx,ny,nz)), ...]
        #
        # important:
        #   a) some indexes are missing
        #   b) increasing firstly by atom_i_index, then by atom_j_index
        Tol = 0.00000001
        HBMAPS = []
        for ia,a in enumerate(MOL.atoms):
            ha = LIB.get_atom_hbond(a[0])
            if ha == 2:     # D1
                for jb,b in enumerate(MOL.atoms):
                    if ia == jb: continue
                    v = [a[t]-b[t] for t in range(2,5)]
                    dd = max(sum([t*t for t in v]), Tol)
                    if dd < 1.90:
                        btype = b[0].upper()
                        bo = False
                        if btype == 'OA' or btype == 'SA':
                            rexp = 4
                            if Disorder_h:
                                bo = True
                        else:
                            rexp = 2
                        dinv = 1.0 / pow(dd,0.5)
                        v = tuple([t*dinv for t in v])
                        HBMAPS.append((ia,ha,jb,rexp,v,bo))
                        #TODO only find one of them??
                        break
            elif ha == 5:   # A2
                nbond = 0
                j1 = 0
                j2 = 0
                for jb,b in enumerate(MOL.atoms):
                    if ia == jb: continue
                    v = [a[t]-b[t] for t in range(2,5)]
                    dd = max(sum([t*t for t in v]), Tol)
                    btype = b[0].upper()
                    if (dd < 3.61 and (btype != 'HD' and btype != 'H')) or \
                        (dd < 1.69 and (btype == 'HD' or btype == 'H')):
                        if nbond == 0:
                            nbond == 1
                            j1 = jb
                        elif nbond == 1:
                            nbond == 2
                            j2 = jb
                        else:
                            logger.warning(
                                'Found an H-bonding atom with three bonded atoms, '
                                'atom serial number {:}'.format(ia+1)
                            )
                if nbond == 0:
                    logger.warning(f'Oxygen atom found with no bonded atoms, atom serial number: {ia+1}')
                elif nbond == 1:        # one bond: carbonyl oxygen: O=C-X
                    v1 = [a[t]-MOL.atoms[j1][t] for t in range(2,5)]
                    dd1 = max(sum([t*t for t in v1]), Tol)
                    d1inv = 1.0 / pow(dd1,0.5)
                    v1 = tuple([t*d1inv for t in v1])
                    v2 = None
                    # try to find second
                    for j2,c in enumerate(MOL.atoms):
                        if j2 == ia or j2 == j1: continue
                        vd = [c[t]-b[t] for t in range(2,5)]
                        dd2 = max(sum([t*t for t in vd]), Tol)
                        btype = c[0].upper()
                        if (dd2 < 2.89 and btype != 'HD') or (dd2 < 1.69 and btype == 'HD'):
                            d2inv = 1.0 / pow(dd2,0.5) 
                            vd = [t*d2inv for t in vd]
                            # C=O cross C-X to get long pair plane norm vector
                            x = v1[1]*vd[2] - v1[2]*vd[1]
                            y = -v1[0]*vd[2] + v1[2]*vd[0]
                            z = v1[0]*vd[1] - v1[1]*vd[0]
                            d = pow(max(x*x+y*y+z*z, Tol), 0.5)
                            v2 = (x/d, y/d, z/d)
                            break
                    if not v2:
                        logger.warning(f'Oxygen atom found only one bonded atom, atom serial number: {ia+1}')
                        v2 = (0.0, 0.0, 0.0)
                    HBMAPS.append((ia,ha,jb,rexp,v1,Disorder_h,v2))
                elif nbond == 2:
                    btype1 = MOL.atoms[j1][0].upper()
                    btype2 = MOL.atoms[j2][0].upper()
                    if Disorder_h and (btype1 == 'HD' or btype2 == 'HD') and (btype1 != btype2):
                        ndx = None
                        if btype1 == 'A' or btype1 == 'C': ndx = j1
                        if btype1 == 'A' or btype2 == 'C': ndx = j2
                        if not ndx:
                            ndx = 0
                            logger.warning('parameters error')
                        v = [a[t]-MOL.atoms[ndx][t] for t in range(2,5)]
                        dd = max(sum([t*t for t in v]), Tol)
                        dinv = 1.0 / pow(dd,0.5)
                        v = tuple([t*dinv for t in v])
                        HBMAPS.append((ia,ha,jb,rexp,v,True))
                    else:
                        v2 = [MOL.atoms[j2][t]-MOL.atoms[j1][t] for t in range(2,5)]
                        dd2 = max(sum([t*t for t in v2]), Tol)
                        d2inv = 1.0 / pow(dd2,0.5)
                        v2 = tuple([t*dinv for t in v2])
                        
                        dot = sum([(a[t+2]-MOL.atom[j1][t+2])*v2[t] for t in range(3)])
                        v1 = [a[t+2] - (dot*v2[t]+MOL.atoms[j1][t+2]) for t in range(3)]
                        dd1 = max(sum([t*t for t in v1]), Tol)
                        d1inv = 1.0 / pow(dd1,0.5)
                        v1 = tuple([t*dinv for t in v1])
                        HBMAPS.append((ia,ha,jb,rexp,v1,False,v2))
                else:
                    logger.critical('How can it be??')



