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
            logger.info('Using a *constant* dielectric: {:}'.format(GPF.gpf['dielectric']))
            EPSTable = [1.0 for i in range(points)]
        else:
            EPSTable = [calc_ddd_Mehler_Solmajer(i/100.0) for i in range(points)]
            logger.info('Using *distance-dependent* dielectric function of Mehler and Solmajer, Prot.Eng.4, 903-910.')
            logger.debug('  d   Dielectric')
            for i in range(0,min(300,points),10):
                logger.debug('{:4.1f}{:9.2f}'.format(EPSTable[i]))

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






