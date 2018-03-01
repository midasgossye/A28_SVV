from main import *
import unittest

# test
# print "stiff location print:", stif_loc(h, t_sk, n_st)
# print "torsional constant", torsional_constant(h, t_sk, C_a)
# testunits for unittests


class TestGeoPropFunctions(unittest.TestCase):
    def test_Xsection(self):
        self.assertEqual(cross_section(0, 0, 0, 0, 0, 0, 0, 0), (0, 0, 0, 0))  # zero test
        # self.assertAlmostEqual(sum(cross_section(h, C_a, t_sk, t_sp, 11, w_st, t_st, h_st)), 0.002, places=None,
        #                        msg=None, delta=(0.002 / 10))  # test data from catia model
        self.assertAlmostEqual(sum(cross_section(h, C_a, t_sk, t_sp, 0, w_st, t_st, h_st)), 0.002,
                               places=3)  # no stiffeners

    def test_enc_area(self):
        self.assertEqual(enc_area(0, 0, 0), (0, 0))  # zero test
        self.assertLess(sum(enc_area(1, 1, 0)), 1)  # selfdone test
        # self.assertAlmostEqual(sum(enc_area(h, C_a, t_sk)), ()) #verification

    # def test_stifloc(self):


suite = unittest.TestLoader().loadTestsFromTestCase(TestGeoPropFunctions)
unittest.TextTestRunner(verbosity=2).run(suite)