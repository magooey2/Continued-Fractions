# -------------------------------------------------------------
# Module for Elementary Number Theory over $\Bbb Q$.
#
# Copyright 2006, 2007, 2008, 2018 Jesse I. Deutsch
# -------------------------------------------------------------


def ord_modulo(a, n):
    """Compute the order of a modulo n.

    Computes the order of rational int a modulo rational
    int n.
    """

    prod = a % n
    for count in range(1, n):
        if prod == 1:
            return count
        else:
            prod = (a * prod) % n

    return "Not relatively prime to n"


####----- end function -----


def solve_quad_mod(a, b, c, n):
    """Solve a quadratic equation modulo n.

    Find all solutions to the quadratic equation
        a*x^2 + b*x + c  = 0 mod n
    for integer n.  Here a, b, c are integers.
    """

    solutions = []
    for x in range(n):
        poly_val = (a * x * x + b * x + c) % n
        if poly_val == 0:
            solutions.append(x)

    return solutions


####----- end function -----


def solve2_quad_mod(a, b, c, d, n):
    """Solve a quadratic equation modulo n.

    Find all solutions to the quadratic equation
        a*x^2 + b*x*y + c*y^2 + d  = 0 mod n
    for integer n.  Here a, b, c, d are integers.
    """

    solutions = []
    nSols = 0
    for x in range(n):
        for y in range(n):
            poly_val = (a * x * x + b * x * y + c * y * y + d) % n
            if poly_val == 0:
                solutions.append((x, y))
                nSols = nSols + 1

    return (nSols, solutions)


####----- end function -----


def euclid_alg(a, b):
    """Euclidean algorithm for gcd.

    Finds the greatest common divisor of
    two integers.
    """

    a, b = abs(a), abs(b)
    while b != 0:
        r = a % b
        a, b = b, r

    return a


####----- end function -----


def ext_euclid_alg(m, n):
    """Extended Euclidean algorithm for gcd.

    Finds the greatest common divisor of
    two integers a and bm and solves for integer
    x, y such that ax + by = 1.  From Knuth, TAOCP
    vol. I, p. 14.
    Variables --
        q, r -- quotient and remainder
        a, b -- am + bn = gcd
        apme, bpme -- a prime and b prime
        t -- temporary
    """

    m, n = abs(m), abs(n)
    q, r = m // n, m % n
    apme = b = 1
    a = bpme = 0

    while r != 0:
        m, n = n, r
        t = apme
        apme = a
        a = t - q * a
        t = bpme
        bpme = b
        b = t - q * b
        """ reset q and r """
        q, r = m // n, m % n

    return (n, a, b)


####----- end function -----


def sigma(n):
    """Sum of divisors of integer n.

    Computes the sum of the positive divisors of the
    integer n.
    """

    n = abs(n)
    sum = n
    for d in range(1, 1 + (n // 2)):
        if n % d == 0:
            sum = sum + d

    return sum


####----- end function -----


def sigma_k(n, k):
    """Sum of divisors of integer n to the power k.

    Computes the sum of the positive divisors of the
    integer n raised to the power k.
    """

    n = abs(n)
    sum = n**k
    for d in range(1, 1 + (n // 2)):
        if n % d == 0:
            sum = sum + d**k

    return sum


####----- end function -----


def phi(n):
    """Euler phi function.

    A simplistic approach to the Euler totient phi function.
    We simply count the number of relatively prime integers between
    1 and n.
    """

    n = abs(n)
    sum = 0
    # the list [1, ..., n]
    for d in range(1, 1 + n):
        if euclid_alg(n, d) == 1:
            sum = sum + 1

    return sum


####----- end function -----


def phi2(n, p=2):
    """Euler phi function.  Second version.

    We use the fact that phi has a formula in terms of the prime
    decomposition of n.  We assume n is a positive integer.  The
    default initial value for prime p is 2.  Thus calls such as
    phi2 (52961) will succeed, though given n is odd we could also
    use phi2 (52961, 3).
    Variables --
                nextp -- the next prime dividing n
                expo -- exponent of nextp exactly dividing n
    """
    # stopping condition
    if n == 1:
        return 1

        # search for next prime dividing n
    nextp = p
    while n % nextp != 0:
        nextp = nextp + 1
        # find exact power of prime in n
    expo = 0
    while n % nextp == 0:
        n, expo = n // nextp, expo + 1

        # formula for phi(p**e)
    factor = nextp**expo - nextp ** (expo - 1)
    return factor * phi2(n, nextp)


####----- end function -----


def lst_bqform_diag_rpns(n, a, b):
    """List of representations of n by binary quad forms.

    Here we deal only with the diagonal case a*x**2 + b*y**2 = n.
    Also, we limit a and b to be > 0, and we are looking only
    for nonnegative solutions.
    Variables --
                nextp -- the next prime dividing n
    """

    # bounds for a and b with a
    # little extra room
    a_bnd = int((1 + n / a) ** 0.5)
    b_bnd = int((1 + n / a) ** 0.5)
    # initialize list
    lst_rpns = []
    # loops
    for x in range(a_bnd + 1):
        for y in range(b_bnd + 1):
            if n == a * x * x + b * y * y:
                lst_rpns.append([x, y])

    return lst_rpns


####----- end function -----


def square_part(n):
    """Returns factorization of n as a square times square-free part.

    Returns a list [nsq, nsq_free] so that
                    n = nsq_free * nsq**2.
    Variables --
                nextp -- the next prime dividing n
    """

    sgn_n = 1
    if n < 0:
        sgn_n = -1
    n = abs(n)
    nsq, nsq_free = 1, 1
    # next prime
    nextp = 2
    # loops
    while n != 1:
        # nextp does not divide n
        if n % nextp != 0:
            nextp += 1
            # nextp**2 divides n
        elif n % (nextp * nextp) == 0:
            nsq = nsq * nextp
            n = n // (nextp * nextp)
            # nextp || n
        else:
            n = n // nextp
            nsq_free = nsq_free * nextp

    return [nsq, sgn_n * nsq_free]


####----- end function -----
