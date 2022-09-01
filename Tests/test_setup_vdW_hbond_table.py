from PyAutoDock.setup_par_library import SetupParLibrary
from PyAutoDock.setup_maps import setup_vdW_hbond_table
from PyAutoDock.utils import read_table_vdW_hbond_data
from utils import allclose

import unittest
import os

CWD = os.path.split(os.path.abspath(__file__))[0]


class TestReadGPF(unittest.TestCase):

    def test_setup_vdW_hbond_table(self):
        library_filename = os.path.join(CWD,'Data/AD4.1_bound.dat')
        LIB = SetupParLibrary(filename=library_filename)
        ligand_types = ['H','HS','HD','N','C','A','OA','SA','OS','NS','NA','S']
        receptor_types = ['H','HS','HD','N','C','A','OA','SA','OS','NS','NA','S']
        ligand_par = [LIB.get_atom_par(i) for i in ligand_types]
        receptor_par = [LIB.get_atom_par(i) for i in receptor_types]
        table = setup_vdW_hbond_table(
            num_points=8300,
            referatomtypes=ligand_types,
            referpar=ligand_par,
            againstatomtypes=receptor_types,
            againstpar=receptor_par,
            gridsize=100,
            num_smooth=25,
            eclamp=100000,
            loglevel=100
        )
        assert table

        mine = read_table_vdW_hbond_data(os.path.join(CWD,'Data/vdw-hbond-table.txt'))
        assert allclose(table,mine)


if __name__ == '__main__':
    unittest.main()



