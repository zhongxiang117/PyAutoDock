from PyAutoDock.read_pdbqt import ReadPDBQT
from utils import allclose

import unittest
import os

CWD = os.path.split(os.path.abspath(__file__))[0]


class TestReadGPF(unittest.TestCase):

    def test_read_pdbqt(self):
        total_atomnum = 3635
        counter = {
            'undefined': 0,
            'A' : 111,
            'C' : 1019,
            'H' : 1387,
            'HD': 426,
            'N' : 343,
            'NA': 3,
            'OA': 343,
            'SA': 3 
        }
        xyzmax = [56.463, 37.979, 47.694]
        xyzmin = [12.118, -16.577, 4.354]
        qtot = -0.928
        qmax = 0.665
        qmin = -0.648

        f = ReadPDBQT(os.path.join(CWD,'Data/receptor-format-1-original.pdbqt'))
        self.assertEqual(total_atomnum,len(f.atoms))
        self.assertEqual(counter,f.counter)
        self.assertTrue(allclose([qtot,qmax,qmin],[f.total_charges,f.qmax,f.qmin],tol=0.001))
        self.assertTrue(allclose(xyzmax,[f.xmax,f.ymax,f.zmax],tol=0.001))
        self.assertTrue(allclose(xyzmin,[f.xmin,f.ymin,f.zmin],tol=0.001))

    def test_read_pdbqt_2(self):
        counter = {
            'undefined': 0,
            'N' : 1,
            'C' : 2,
            'OA': 1,
            'H' : 2,
            'HS': 4,
        }
        f = ReadPDBQT(os.path.join(CWD,'Data/receptor-format-2-good.pdbqt'))
        self.assertEqual(counter,f.counter)

    def test_read_pdbqt_centroid(self):
        f = ReadPDBQT(os.path.join(CWD,'Data/receptor-format-1-original.pdbqt'))
        f.centroid(32.961,7.860,27.277)
        self.assertTrue(allclose(f.center, [32.961,7.860,27.277], tol=0.001))

if __name__ == '__main__':
    unittest.main()


