from PyAutoDock.logger import mylogger
from PyAutoDock.read_dpf import ReadDPF
from PyAutoDock.read_pdbqt import ReadPDBQT
from PyAutoDock.setup_par_library import SetupParLibrary
from PyAutoDock.utils import file_gen_new, calc_ddd_Mehler_Solmajer
from PyAutoDock.setup_maps import GridMap
from PyAutoDock.main_autogrid import SetupGridMaps

import math

logger = mylogger()
if logger.level == 0: logger.level = 20     # NOTSET

SELECTION_MODE = {0:'Proportional', 1:'LinearRanking', 2:'Tournament', 3:'Boltzmann'}


def index2distance(i):
    return pow(i/32,0.5)

def setup_eps_square_table(num_points, gridsize=None,is_distance=True,loglevel=None):
    if not gridsize: gridsize = 100
    if not is_distance:
        return [1.0 for i in range(num_points*gridsize)]

    ds = [index2distance(i) for i in range(num_points)]     #TODO update number 32
    EPSTable = [calc_ddd_Mehler_Solmajer(ds[i]) for i in range(num_points)]
    EPSTable[0] = 1.0   # important
    REPSTable = [v*EPSTable[i] for i,v in enumerate(ds)]
    if loglevel and loglevel <= 10:
        print('\nUsing Square *distance-dependent* dielectric function of Mehler and Solmajer, Prot.Eng.4, 903-910.')
        print(' i   distance   Eeps     r*Epes')
        for i in range(0,min(310,num_points),10):
            print('{:<5} {:>5.3f} {:>9.3f} {:>9.3f}'.format(i,ds[i],EPSTable[i],REPSTable[i]))
        print('\n')
    return EPSTable,REPSTable


class SetupDockMaps:
    def __init__(
            self,library_filename,ligand_dpf_filename,
            nbc=None,eclamp=None,tol=None,gsize=None,eintcutoff=None,
            *args,**kwargs
        ):
        self.library_filename = library_filename
        self.ligand_dpf_filename = ligand_dpf_filename
        self.nbc = nbc if nbc else 8
        self.eclamp = eclamp if eclamp else 100000
        self.tol = tol if tol else 0.00000001
        self.gsize = gsize if gsize else 100
        self.eintcutoff = eintcutoff if eintcutoff else 2047
        self.gpf = None
        self.EnergyMaps = None
        self.gen(tol=self.tol)

    def gen(self,tol=None,loglevel=None):
        if not loglevel: loglevel = 10
        tol = tol if tol else 0.00000001
        LIB = SetupParLibrary(self.library_filename)
        if not len(LIB.atoms): return
        DPF = ReadDPF(self.ligand_dpf_filename,loglevel=100)

        for a in DPF.dpf['ligand_types']:
            if not LIB.get_atom_par(a):
                logger.critical('not defined: receptor atomtype: {:}'.format(a))

        num_points = 2000

        MAPS = []
        num_maps = len(DPF.dpf['ligand_types'])
        for i in range(num_maps):
            m = GridMap(num_receptor_maps=num_maps)
            m.atomtype = DPF.dpf['ligand_types'][i]
            par = LIB.get_atom_par(DPF.dpf['ligand_types'][i])
            m.sol = par.sol
            m.vol = par.vol
            m.Rij = par.Rij
            m.epsij = par.epsij
            m.hbond = par.hbond
            m.Rij_hb = par.Rij_hb
            m.epsij_hb = par.epsij_hb
            m.num_cmaps = num_maps - i
            for j in range(num_maps-i):
                m.catomtypes[j] = DPF.dpf['ligand_types'][j+i]      # index j+i
                p = LIB.get_atom_par(DPF.dpf['ligand_types'][j+i])  # index j+i
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
                else:
                    m.nbp_r[j] = (m.Rij + p.Rij) / 2.0
                    m.nbp_eps[j] = pow(m.epsij*p.epsij, 0.5)
                tmp = m.nbp_eps[j] / (m.xA[j]-m.xB[j])
                m.cA[j] = tmp * pow(m.nbp_r[j],m.xA[j]) * m.xB[j]
                m.cB[j] = tmp * pow(m.nbp_r[j],m.xB[j]) * m.xA[j]
            # one level up
            MAPS.append(m)

        gridsize = 100
        # smooth
        num_smooth = int(DPF.dpf['smooth'] * gridsize / 2.0 + 0.6)

        if loglevel and loglevel <= 10:
            txt = 'Before Smooth: ' if num_smooth and num_smooth > 0 else ''
            for i in range(len(MAPS)):       # index i is super important, be careful
                print('\n{:}Pairwise interactions:'.format(txt))
                out = '>> {:>3}-'.format(MAPS[i].atomtype)
                cn = MAPS[i].num_cmaps
                for j in range(cn):
                    tmpa = 'cA={:.4f} / {:}'.format(MAPS[i].cA[j],MAPS[i].xA[j])
                    tmpb = 'cB={:.4f} / {:}'.format(MAPS[i].cB[j],MAPS[i].xB[j])
                    print(out+'{:<3}   {:30}   {:}'.format(MAPS[i].catomtypes[j],tmpa,tmpb))

        EvdWHBondTable = [
            [[0.0 for k in range(num_maps)] for j in range(num_maps)]
            for i in range(num_points)
        ]
        EvdWHNonBondTable = [
            [[0.0 for k in range(num_maps)] for j in range(num_maps)]
            for i in range(num_points)
        ]
        for i in range(1,num_points):     # index i starts at 1
            r = index2distance(i)       #TODO
            for j,m in enumerate(MAPS):
                for n in range(m.num_cmaps):
                    rA = pow(r, m.xA[n])
                    rB = pow(r, m.xB[n])
                    EvdWHBondTable[i][n][j] = min(self.eclamp, m.cA[n]/rA-m.cB[n]/rB)
                    EvdWHBondTable[i][j][n] = EvdWHBondTable[i][n][j]
                    EvdWHNonBondTable[i][n][j] = min(self.eclamp, 392586.8/rA) - r      #TODO
                    EvdWHNonBondTable[i][j][n] = EvdWHNonBondTable[i][n][j]
                # clamp them
                EvdWHBondTable[0][j][n] = EvdWHBondTable[0][n][j] = self.eclamp
                EvdWHNonBondTable[0][j][n] = EvdWHNonBondTable[0][n][j] = self.eclamp
        
        infolist = self.printminvalue(EvdWHBondTable,MAPS)
        for i in infolist: print(i)



    def printminvalue(self,table,maps):
        infolist = []
        for i in range(len(maps)):
            for j in range(maps[i].num_cmaps):
                dmin = self.eclamp
                beg = -1
                for n in range(len(table)):
                    if table[n][j][i] < dmin:
                        dmin = table[n][j][i]
                        beg = n
                end = beg
                for m in range(len(table)-1,-1,-1):
                    if table[m][j][i] <= dmin:
                        end = m
                        break
                infolist.append(
                    'atomtypes: {:} - {:},  minimum value = {:.6f},  distance range: '
                    '{:.6f} ~ {:.6f}'.format(maps[i].atomtype, maps[i].catomtypes[j],
                    dmin,index2distance(beg),index2distance(end))
                )
        return infolist







