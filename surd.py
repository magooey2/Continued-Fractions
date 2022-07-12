# -------------------------------------------------------------
# Module for Quardratic Surds over $\Bbb Q$.
#
# Copyright 2006 Jesse I. Deutsch
#
#  $Id:$
# -------------------------------------------------------------


from elem_nt import euclid_alg, square_part


class Surd:
    """Surd, class for quadratic surds over $\Bbb Q$

    We are essentially dealing with the field Q(sqrt(r))
    with elements in the form
                                    a + b*sqrt(r)
                                    ------------- .
                                            d
    We implement common arithmetic functions, check for
    consistency in the value of r, and will attempt to
    construct the continued fraction expansion of positive
    surds.
    """

    def __init__(self, a, b, d, r):
        self.a = a
        self.b = b
        self.d = d
        self.r = r
        if d == 0:
            raise ZeroDivisionError
        Surd.normalize(self)

        # string version of object
        # useful for print

    def __str__(self):
        sgn = "+"
        self.b_abs = self.b
        if self.b < 0:
            sgn = "-"
            self.b_abs = -1 * self.b
            # case where denom is not 1
        if self.d != 1:
            return (
                "["
                + repr(self.a)
                + sgn
                + repr(self.b_abs)
                + "*rt("
                + repr(self.r)
                + ")] / "
                + repr(self.d)
            )
            # case where denominator is 1
        return repr(self.a) + sgn + repr(self.b_abs) + "*rt(" + repr(self.r) + ")"

    def normalize(self):
        if self.d < 0:
            self.d = -1 * self.d
            self.a = -1 * self.a
            self.b = -1 * self.b
            # now must clean out gcd
            # a and b are both zero
        if self.a == 0:
            if self.b == 0:
                self.d = 1
                return
                # a is zero, b not zero
            else:
                self.gcd = euclid_alg(self.b, self.d)
                self.b = self.b // self.gcd
                self.d = self.d // self.gcd
                return
                # b is zero, a not zero
        if self.b == 0:
            self.gcd = euclid_alg(self.a, self.d)
            self.a = self.a // self.gcd
            self.d = self.d // self.gcd
            return
            # both a and b are not zero
        self.gcd = euclid_alg(self.a, self.b)
        self.gcd = euclid_alg(self.gcd, self.d)
        self.a = self.a // self.gcd
        self.b = self.b // self.gcd
        self.d = self.d // self.gcd
        return

    def add(self, y):
        self.temp_a = self.a * y.d + self.d * y.a
        self.temp_b = self.b * y.d + self.d * y.b
        self.temp_d = self.d * y.d
        self.a = self.temp_a
        self.b = self.temp_b
        self.d = self.temp_d
        Surd.normalize(self)

    def add_replace(self, x, y):
        self.temp_a = x.a * y.d + x.d * y.a
        self.temp_b = x.b * y.d + x.d * y.b
        self.temp_d = x.d * y.d
        self.a = self.temp_a
        self.b = self.temp_b
        self.d = self.temp_d
        Surd.normalize(self)

    def sub(self, y):
        self.temp_a = self.a * y.d - self.d * y.a
        self.temp_b = self.b * y.d - self.d * y.b
        self.temp_d = self.d * y.d
        self.a = self.temp_a
        self.b = self.temp_b
        self.d = self.temp_d
        Surd.normalize(self)

    def sub_replace(self, x, y):
        self.temp_a = x.a * y.d - x.d * y.a
        self.temp_b = x.b * y.d - x.d * y.b
        self.temp_d = x.d * y.d
        self.a = self.temp_a
        self.b = self.temp_b
        self.d = self.temp_d
        Surd.normalize(self)

        # need temps

    def mult(self, y):
        self.temp_a = self.a * y.a + self.r * self.b * y.b
        self.temp_b = self.a * y.b + self.b * y.a
        self.temp_d = self.d * y.d
        self.a = self.temp_a
        self.b = self.temp_b
        self.d = self.temp_d
        Surd.normalize(self)

    def mult_replace(self, x, y):
        # need temp as x.rtnl used
        # on next line.  Resolves
        # a.mult_replace(a,a) bug.
        self.temp_a = x.a * y.a + x.r * x.b * y.b
        self.temp_b = x.a * y.b + x.b * y.a
        self.temp_d = x.d * y.d
        self.a = self.temp_a
        self.b = self.temp_b
        self.d = self.temp_d
        Surd.normalize(self)

        # need temps

    def div(self, y):
        self.temp_a = self.b * y.b * y.d * y.r - self.a * y.a * y.d
        self.temp_b = y.d * (self.a * y.b - y.a * self.b)
        self.temp_d = self.d * (y.b * y.b * y.r - y.a * y.a)
        self.a = self.temp_a
        self.b = self.temp_b
        self.d = self.temp_d
        Surd.normalize(self)

    def div_replace(self, x, y):
        # need temps.  Resolves
        # a.mult_replace(a,a) bug.
        self.temp_a = x.b * y.b * y.d * y.r - x.a * y.a * y.d
        self.temp_b = y.d * (x.a * y.b - y.a * x.b)
        self.temp_d = x.d * (y.b * y.b * y.r - y.a * y.a)
        self.a = self.temp_a
        self.b = self.temp_b
        self.d = self.temp_d
        Surd.normalize(self)

        # conjugate wrt Q(sqrt(r))

    def conj(self, x):
        self.a = x.a
        self.b = -1 * x.b
        self.d = x.d
        self.r = x.r

    def norm(self):
        self.temp_a = self.a * self.a - self.r * self.b * self.b
        self.temp_b = 0
        self.temp_d = self.d * self.d
        self.temp_norm = Surd(self.temp_a, self.temp_b, self.temp_d, self.r)
        return self.temp_norm

    def copy_in(self, x):
        self.a = x.a
        self.b = x.b
        self.c = x.c
        self.r = x.r

    def copy_out(self, x):
        x.a = self.a
        x.b = self.b
        x.c = self.c
        x.d = self.d

        # FIX!!! check if totally positive, >> 0

    def is_tot_pos(self):
        if self.rtnl <= 0:
            return False
        if self.norm() <= 0:
            return False
        return True


# ----- end of class -------


def isqrt(a):
    """Integer square root.

    Argument is assumed to be a positive integer.  We use the
    float square root when the argument is <= 10^8.  If larger
    we divide by an even power of 10, and then get the float
    square root of the integer part, and multiply it back by
    the appropriate power of 10.  Then we use the estimate
           sqrt(ap^2 + b) \approx  ap + b/(2*ap)
    which is actually equivalent to the ap <-- (ap + x/ap)/2
    iteration for the square root of x.  However, since we
    have the first few digits of the square root, convergence
    is very quick - 4 or 5 iterations for a 50 digit number.
    """

    if a < 1 + 10**8:
        return int(a**0.5)
    else:
        b = len(str(a))
        c = b // 8
        d = a // 10 ** (8 * c - 4)
        approx = int(d**0.5) * 10 ** (4 * c - 2)
        # 						----- Debugging -----
        # 	print b, c, d
        # 	print approx
        # 						------ end debugging -----
        while 1:
            b = a - approx * approx
            approx2 = approx + (b // (2 * approx))
            # 						----- Debugging -----
            # 	print approx2
            # 						----- end debugging -----
            if approx2 == approx:
                return approx
            approx = approx2


####----- end function -----


def fund_unit(m):
    """Returns fundamental unit of Q(sqrt(m)).

    Note that m must be square free and > 1.  We use the
    elementary 	algorithm in Marcus's book.
    Variables --
                            test_4_square --
    """

    if m <= 1:
        return "ERROR: m must be greater than 1."
    m_sq_fctr = square_part(m)
    if m_sq_fctr[0] > 1:
        return "ERROR: m must be square free"

        # Case of m is 2, 3 mod 4.  The
        # algorithm-- check mb^2 +-1 for
        # smallest positive b that makes
        # a perfect square, a^2.  Then
        # a + b sqrt(m) is the fund. unit.
    if m % 4 == 2 or m % 4 == 3:
        b = 1
        while 1:
            test_4_square = m * b * b - 1
            test_a = isqrt(test_4_square)
            if test_4_square - test_a**2 == 0:
                return Surd(test_a, b, 1, m)
            test_4_square = m * b * b + 1
            test_a = isqrt(test_4_square)
            if test_4_square - test_a**2 == 0:
                return Surd(test_a, b, 1, m)
            b += 1
            # Case of m = 1 mod 4
    else:
        b = 1
        while 1:
            test_4_square = m * b * b - 4
            test_a = isqrt(test_4_square)
            if test_4_square - test_a**2 == 0:
                return Surd(test_a, b, 2, m)
            test_4_square = m * b * b + 4
            test_a = isqrt(test_4_square)
            if test_4_square - test_a**2 == 0:
                return Surd(test_a, b, 2, m)
            b += 1


####----- end function -----
