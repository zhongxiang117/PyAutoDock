from PyAutoDock.logger import mylogger
from PyAutoDock.read_gpf import ReadGPF
from PyAutoDock.read_pdbqt import ReadPDBQT
from PyAutoDock.setup_par_library import SetupParLibrary
from PyAutoDock.utils import file_gen_new, calc_ddd_Mehler_Solmajer

import math

logger = mylogger()
if logger.level == 0: logger.level = 20     # NOTSET

SELECTION_MODE = {0:'Proportional', 1:'LinearRanking', 2:'Tournament', 3:'Boltzmann'}
