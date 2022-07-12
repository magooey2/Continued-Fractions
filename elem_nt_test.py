# -------------------------------------------------------
# elem_nt_test -- unit tests for elem_nbr_theory module
#
# Copyright 2007 Jesse I. Deutsch
# -------------------------------------------------------


import unittest
from elem_nt import *


class Elem_NT_Tests(unittest.TestCase):
    def testPhi(self):
        ans = phi(52961)
        self.failUnless(ans == 52500)

    def testPhi2(self):
        ans = phi2(52961)
        self.failUnless(ans == 52500)

    def testPhiFalse(self):
        ans = phi(97)
        self.failIf(ans == 99)

    def testPhi2False(self):
        ans = phi2(89)
        self.failIf(ans != 88)

    def testOrd_Modulo(self):
        ans = ord_modulo(2, 211)
        self.failUnless(ans == 210)
        ans = ord_modulo(3, 97)
        self.failUnless(ans == 48)
        ans = ord_modulo(7, 97)
        self.failUnless(ans == 96)

    def testEuclid_Alg(self):
        ans = euclid_alg(112, 211)
        self.failUnless(ans == 1)
        ans = euclid_alg(33, 96)
        self.failUnless(ans == 3)
        ans = euclid_alg(502, 52961)
        self.failUnless(ans == 251)

    def testSigma(self):
        ans = sigma(6)
        self.failUnless(ans == 12)
        ans = sigma(28)
        self.failUnless(ans == 56)
        ans = sigma(7)
        self.failUnless(ans == 8)
        ans = sigma(496)
        self.failUnless(ans == 2 * 496)
        ans = sigma(12)
        self.failUnless(ans == 28)
        ans = sigma(101)
        self.failUnless(ans == 102)

    def testSolve_Quad_Mod(self):
        ans = solve_quad_mod(1, 0, 1, 4)
        self.failUnless(ans == [])
        ans = solve_quad_mod(1, 0, -1, 4)
        self.failUnless(ans == [1, 3])
        ans = solve_quad_mod(1, 0, -1, 8)
        self.failUnless(ans == [1, 3, 5, 7])
        ans = solve_quad_mod(1, 0, 1, 7)
        self.failUnless(ans == [])
        ans = solve_quad_mod(1, 0, 1, 17)
        self.failUnless(ans == [4, 13])
        ans = solve_quad_mod(1, 0, 1, 97)
        self.failUnless(ans == [22, 75])

    def testLst_Bqform_Diag_Rpns(self):
        ans = lst_bqform_diag_rpns(3, 1, 1)
        self.failUnless(ans == [])
        ans = lst_bqform_diag_rpns(5, 1, 1)
        self.failUnless(ans == [[1, 2], [2, 1]])
        ans = lst_bqform_diag_rpns(29, 1, 1)
        self.failUnless(ans == [[2, 5], [5, 2]])
        ans = lst_bqform_diag_rpns(97, 1, 2)
        self.failUnless(ans == [[5, 6]])


def main():
    unittest.main()


if __name__ == "__main__":
    main()
