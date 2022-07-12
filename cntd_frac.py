# -------------------------------------------------------------------------------
#  Module for classical Continued Fractions.
#
#  Copyright 2015 Jesse I. Deutsch
#
#   $Id:$
# -------------------------------------------------------------------------------


from elem_nt import euclid_alg, square_part
from surd import Surd
import linfractrans


def cf_finite(x, bnd):
    """Continued fraction expansion with bound.

    We take a Surd argument and compute its continued fraction
    expansion up to bnd convergents.  We test for a zero
    remainder as the input could be a rational with a small
    expansion.
    """

    cf_list = []
    cf = Surd(x.a, x.b, x.d, x.r)
    cnt = 0
    while cnt < bnd:
        # take floor, then reciprocal
        # r**0.5 won't work for large r
        cf_top_approx = cf.a + cf.b * cf.r**0.5
        # __floordiv__ is `//'
        cf_floor = cf_top_approx // cf.d
        cf_list.append(int(cf_floor))
        # get into range [0,1)
        cf.a = cf.a - cf_floor * cf.d
        # zero possible in rational case
        if cf.a == 0 and cf.b == 0:
            return cf_list

        cf.div_replace(Surd(1, 0, 1, cf.r), cf)
        cnt = cnt + 1
    return cf_list


####----- end function -----


def cntd_frac(x):
    """Continued fraction expansion of quadratic Surd.

    We take a Surd argument and compute its continued fraction
    expansion.  We test for a zero remainder as the input could
    be a rational with a small expansion.  The other case is a
    repeating continued fraction expansion.  We find and
    delineate the first repeating section of the continued
    fraction for this Surd.
    """

    cf_list = []
    cf = Surd(x.a, x.b, x.d, x.r)

    # case of rational number
    if x.b == 0:
        while 1:
            # take floor, then reciprocal
            # __floordiv__ is `//'
            cf_floor = cf.a // cf.d
            # int converts e.g. 5.0 to 5
            cf_list.append(int(cf_floor))
            # get into range [0,1)
            cf.a = cf.a - cf_floor * cf.d
            # zero means end of contd frac
            if cf.a == 0:
                return cf_list

            cf.div_replace(Surd(1, 0, 1, cf.r), cf)

            # irrational quadratic surd
    else:
        cf_partial_list = []
        cf_partial_list.append(str(cf))

        while 1:
            # take floor, then reciprocal
            # r**0.5 won't work for large r
            cf_top_approx = cf.a + cf.b * cf.r**0.5
            # __floordiv__ is `//'
            cf_floor = cf_top_approx // cf.d
            cf_floor = int(cf_floor)
            cf_list.append(cf_floor)
            # get into range [0,1)
            cf.a = cf.a - cf_floor * cf.d

            cf.div_replace(Surd(1, 0, 1, cf.r), cf)
            idx = -1
            #                       -----Debugging--------
            #  print " ### cf_partial_list: ", cf_partial_list
            #
            for string in cf_partial_list:
                idx = idx + 1
                if string == str(cf):
                    cf_list_repeat = cf_list[:idx]
                    cf_list_repeat.append(cf_list[idx:])
                    return cf_list_repeat
                    # not yet in the list
            cf_partial_list.append(str(cf))


####----- end function -----


def cflist_to_rtnl(cflist, r):
    """Convert finite continued fraction to rational.

    We take a list, interpreted as a finite continued fraction
    expansion, and return the corresponding rational number.
    The output is a Surd, however.
    Variables -- r is a dummy for the quadratic term.
    """

    # initialize with bottom int
    # and work your way to the top.
    # This accounts for cflist[-1]
    # and the indices 2, etc.
    rtnl = Surd(cflist[-1], 0, 1, r)
    for i in range(2, len(cflist) + 1):
        rtnl.div_replace(Surd(1, 0, 1, r), rtnl)
        rtnl.add(Surd(cflist[-i], 0, 1, r))

    return rtnl


####----- end function -----


def cflist_pureperiod_to_surd(cflist):
    """Convert purely periodic continued fraction to surd.

    We take a list, interpreted as a purely periodic continued
    fraction expansion, and return the corresponding surd.
    The strategy is to work with what are essentially linear
    fractional transformations, and find the fixed points of
    them.
            Variables -- none.
    """

    # initialize at bottom to x
    temp = linfractrans.LFT(1, 0, 0, 1)
    # start with bottom int and
    # work to the top.  This accounts
    # for cflist[-i].
    for i in range(1, len(cflist) + 1):
        temp.reciprocal()
        temp.add_const(cflist[-i])

        # get the 2 fixed points
    pp_surd = temp.fixed_pts()
    # return larger of fixed points
    return pp_surd[0]


####----- end function -----


def cflist_to_surd(cflist):
    """Convert finite continued fraction to surd or rational.

    We take a list, interpreted as a finite continued fraction
    expansion, and return the corresponding value.
            Variables -- r is a dummy for the quadratic term.
    """

    # Check last element of list.
    # if integer then its a finite
    # continued fraction,. If a list
    # then it is a surd.
    # Section 5.15 library ref, types
    # an int is an instance of 'int'.
    if isinstance(cflist[-1], int):
        return cflist_to_rtnl(cflist, 2)

        ## purely periodic part at bottom
    elif isinstance(cflist[-1], list):
        cf = cflist_pureperiod_to_surd(cflist[-1])
        cf_r = cf.r
        # and work your way to the top.
        # This accounts for cflist[-1]
        # and the indices 2, etc.
        for i in range(2, len(cflist) + 1):
            cf.div_replace(Surd(1, 0, 1, cf_r), cf)
            cf.add(Surd(cflist[-i], 0, 1, cf_r))

        return cf

    else:
        print("### Type error")
        return


####----- end function -----
