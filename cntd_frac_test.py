# -----------------------------------------------------------
# cntd_frac_test -- unit tests for continued fraction and
#                   surd module.
#
# Copyright 2007, 2017 Jesse I. Deutsch
# -----------------------------------------------------------


import unittest
from cntd_frac import *


class Surd_Tests(unittest.TestCase):
    def testIsqrt(self):
        ans = isqrt(170321)
        self.failUnless(ans == 412)
        ans = isqrt(2349870211234)
        self.failUnless(ans == 1532928)
        ans = isqrt(271828182845904523536314159265358979323846)
        self.failUnless(ans == 521371444217943838413)
        ans = isqrt(87872384752987548725983479287298472897)
        self.failUnless(ans == 9374027136347939372)

    def testCfFinite(self):
        tau = Surd(1, 1, 2, 5)
        ans = cf_finite(tau, 4)
        self.failUnless(ans == [1, 1, 1, 1])
        a = Surd(81, 0, 43, 2)
        ans = cf_finite(a, 4)
        self.failUnless(ans == [1, 1, 7, 1])
        ans = cf_finite(a, 6)
        self.failUnless(ans == [1, 1, 7, 1, 1, 2])
        a = Surd(2, 1, 1, 2)
        ans = cf_finite(a, 5)
        self.failUnless(ans == [3, 2, 2, 2, 2])
        a = Surd(0, 1, 1, 97)
        ans = cf_finite(a, 12)
        self.failUnless(ans == [9, 1, 5, 1, 1, 1, 1, 1, 1, 5, 1, 18])

    def testCntdFrac(self):
        tau = Surd(1, 1, 2, 5)
        ans = cntd_frac(tau)
        self.failUnless(ans == [[1]])
        a = Surd(81, 0, 43, 2)
        ans = cntd_frac(a)
        self.failUnless(ans == [1, 1, 7, 1, 1, 2])
        rt2 = Surd(0, 1, 1, 2)
        ans = cntd_frac(rt2)
        self.failUnless(ans == [1, [2]])
        rt97 = Surd(0, 1, 1, 97)
        ans = cntd_frac(rt97)
        self.failUnless(ans == [9, [1, 5, 1, 1, 1, 1, 1, 1, 5, 1, 18]])
        b = Surd(97, 0, 22, 2)
        ans = cntd_frac(b)
        self.failUnless(ans == [4, 2, 2, 4])

    def testCflistToRtnl(self):
        ans = cflist_to_rtnl([1, 1, 1, 1], 2)
        self.failUnless(str(ans) == str(Surd(5, 0, 3, 2)))
        ans = cflist_to_rtnl([1, 1, 1, 1, 1], 2)
        self.failUnless(str(ans) == str(Surd(8, 0, 5, 2)))

    def testCflistPurePeriodToSurd(self):
        ans = cflist_pureperiod_to_surd([1, 1, 1, 1])
        self.failUnless(str(ans) == str(Surd(1, 1, 2, 5)))
        ans = cflist_pureperiod_to_surd([1, 2])
        self.failUnless(str(ans) == str(Surd(1, 1, 2, 3)))
        ans = cflist_pureperiod_to_surd([3, 5, 2])
        self.failUnless(str(ans) == str(Surd(15, 1, 11, 401)))

    def testCflistToSurd(self):
        ans = cflist_to_surd([1, [2]])
        self.failUnless(str(ans) == str(Surd(0, 1, 1, 2)))
        ans = cflist_to_surd([[1, 2]])
        self.failUnless(str(ans) == str(Surd(1, 1, 2, 3)))
        ans = cflist_to_surd([9, [1, 5, 1, 1, 1, 1, 1, 1, 5, 1, 18]])
        self.failUnless(str(ans) == str(Surd(0, 1, 1, 97)))
        # 2 is default radix
        ans = cflist_to_surd([1, 1, 7, 1, 1, 2])
        self.failUnless(str(ans) == str(Surd(81, 0, 43, 2)))


def main():
    unittest.main()


if __name__ == "__main__":
    main()
