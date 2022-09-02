from PyAutoDock.main_autogrid import SetupGridMaps
from PyAutoDock.utils import read_map_data
from utils import allclose

import unittest
import os


CWD = os.path.split(os.path.abspath(__file__))[0]


class TestMainAutoGrid(unittest.TestCase):

    def test_main_autogrid(self):
        maps = SetupGridMaps(
            library_filename=os.path.join(CWD,'Data/AD4.1_bound.dat'),
            ligand_gpf_filename=os.path.join(CWD,'Data/test_gridmaps/grid.gpf'),
            receptor_mol_filename=os.path.join(CWD,'Data/receptor-format-1-original.pdbqt'),
            loglevel=100
        )
        assert maps.EnergyMaps
        mine = read_map_data(os.path.join(CWD,'Data/test_gridmaps/maps.txt'))
        assert allclose(maps.EnergyMaps,mine,tol=0.001)


if __name__ == '__main__':
    unittest.main()




