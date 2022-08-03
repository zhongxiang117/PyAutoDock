from PyAutoDock.logger import mylogger
from PyAutoDock.read_gpf import ReadGPF
from PyAutoDock.read_pdbqt import ReadPDBQT
from PyAutoDock.setup_par_library import SetupParLibrary
from PyAutoDock.utils import file_gen_new

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
        self.vol_probe = 0.0
        self.solpar_probe = 0.0
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
        self.hbonder = [0.0 for i in range(num_receptor_maps)]


class SetupMolMaps:
    def __init__(self,library_filename,ligand_gpf_filename,receptor_mol_filename=None,*args,**kwargs):
        self.library_filename = library_filename
        self.receptor_mol_filename = receptor_mol_filename if receptor_mol_filename else None
        self.ligand_gpf_filename = ligand_gpf_filename
        self.map = []

    def gen(self):
        LIB = SetupParLibrary(self.library_filename)
        if not len(LIB.atoms): return
        GPF = ReadGPF(self.ligand_gpf_filename)

        if self.receptor_mol_filename:
            MOL = ReadPDBQT(self.receptor_mol_filename)
            if not len(MOL.atoms): return
        elif self.GPF.gpf['receptor']:
            MOL = ReadPDBQT(self.GPF.gpf['receptor'])
            if not len(MOL.atoms): return
        else:
            logger.critical('required: PDBQT file: receptor')
        
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

        if GPF.gpf['receptor_types']:
            defined = set(GPF.gpf['receptor_types'])
            if MOL.atomset.issubset(defined):
                left = defined - MOL.atomset
                logger.warning('not defined: receptor_types: {:}'.format(left))
            elif MOL.atomset.issuperset(defined):
                left = MOL.atomset - defined
                logger.warning('not included: receptor_types: {:}'.format(left))
        else:
            GPF.gpf['receptor_types'] = list(MOL.atomset)
            GPF.gpf['map'] = []
            for a in GPF.gpf['receptor_types']:
                GPF.gpf['map'].append(file_gen_new('receptor.{:}.map'.format(a)))
            logger.info('receptor_types: {:}'.format(GPF.gpf['receptor_types']))
            logger.info('map: {:}'.format(GPF.gpf['map']))







