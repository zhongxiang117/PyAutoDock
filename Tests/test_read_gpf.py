from PyAutoDock.read_gpf import ReadGPF

import unittest
import os

CWD = os.path.split(os.path.abspath(__file__))[0]

class TestReadGPF(unittest.TestCase):

    def test_read(self):
        GPF_0 = ReadGPF(os.path.join(CWD,'Data/grid-format-0-original.gpf'))
        GPF_1 = ReadGPF(os.path.join(CWD,'Data/grid-format-1-simplest.gpf'))
        GPF_2 = ReadGPF(os.path.join(CWD,'Data/grid-format-2-goodfull.gpf'))
        GPF_3 = ReadGPF(os.path.join(CWD,'Data/grid-format-3-goodpart.gpf'))
        # test on file 0 and 1
        self.assertEqual(GPF_0.gpf['gridfld'], GPF_1.gpf['gridfld'])
        # test on file 0 and 2
        self.assertEqual(GPF_0.gpf['npts'], GPF_2.gpf['npts'])
        self.assertEqual(GPF_0.gpf['spacing'], GPF_2.gpf['spacing'])
        self.assertEqual(GPF_0.gpf['receptor_types'], GPF_2.gpf['receptor_types'])
        self.assertEqual(GPF_0.gpf['map'], GPF_2.gpf['map'])
        self.assertEqual(GPF_0.gpf['dielectric'], GPF_2.gpf['dielectric'])
        # test on file 0 and 3
        self.assertEqual(GPF_0.gpf['spacing'], GPF_3.gpf['spacing'])
        self.assertEqual(GPF_3.gpf['gridcenter'], None)
        g_rt_0 = set(GPF_0.gpf['receptor_types'])
        g_rt_3 = set(GPF_3.gpf['receptor_types'])
        self.assertTrue(g_rt_0.issuperset(g_rt_3))
        g_map_0 = set(GPF_0.gpf['map'])
        g_map_3 = set(GPF_3.gpf['map'])
        self.assertTrue(g_map_0.issuperset(g_map_3))


if __name__ == '__main__':
    unittest.main()


