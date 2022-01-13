#-------------------------------------------------------------
# Module for Linear Fractional Transformation over $\Bbb Z$.
#
# Copyright 2007 Jesse I. Deutsch
#
#   $Id:$
#-------------------------------------------------------------



from elem_nt import euclid_alg, square_part
import surd 


class LFT:
		""" Class for linear fractional transformations over $\Bbb Z$

		We are dealing with expressions of the form
						a*z + b
						------- .
						c*z + d
		We implement addition of constants, taking reciprocals,
		and should do composition of these functions.  Also we
		want to find the fixed points of these functions.
		"""

		def __init__ (self, a, b, c, d):
			self.a = a
			self.b = b
			self.c = c
			self.d = d
			if c == 0 and d == 0:
				raise ZeroDivisionError
			LFT.normalize (self)
             
										# string version of object
										# useful for print
		def __str__(self):
			sgn_b = '+'
			self.b_abs = self.b
			if self.b < 0:
				sgn_b = '-'
				self.b_abs = -1 * self.b

			sgn_d = '+'
			self.d_abs = self.d
			if self.d < 0:
				sgn_d = '-'
				self.d_abs = -1 * self.d

										# case where denom is not 1
			return "(" + repr (self.a) + "z" + sgn_b + repr(self.b_abs) \
				+ ") / " \
				+ "(" + repr (self.c) + "z" + sgn_d + repr(self.d_abs) + ")"


		def normalize (self):
			if self.c < 0:
				self.d = -1 * self.d
				self.a = -1 * self.a
				self.b = -1 * self.b
				self.c = -1 * self.c

			if self.c == 0 and self.d < 0:
				self.d = -1 * self.d
				self.a = -1 * self.a
				self.b = -1 * self.b
				self.c = -1 * self.c

										# now must clean out gcd
			self.gcd = euclid_alg (self.a, self.b)
			self.gcd = euclid_alg (self.gcd, self.d)
			self.gcd = euclid_alg (self.gcd, self.c)
			if self.gcd == 0:
				return
			self.a   = self.a // self.gcd
			self.b   = self.b // self.gcd
			self.c   = self.c // self.gcd
			self.d   = self.d // self.gcd
			return										


		def add_const (self, y):
			self.temp_a = self.a + self.c * y
			self.temp_b = self.b + self.d * y
			self.a = self.temp_a
			self.b = self.temp_b
			LFT.normalize (self)

		def reciprocal (self):
			self.temp_a = self.c
			self.temp_b = self.d
			self.c = self.a
			self.d = self.b
			self.a = self.temp_a
			self.b = self.temp_b
			LFT.normalize (self)

		def fixed_pts (self):
			""" Returns a list of the fixed points.
			
			There is a pair of fixed points, which satisfy the
			equation
						cz^2 + (d-a)z - b = 0.
			The discriminant is disc = (d-a)^2 + 4 bc, and the
			roots are ( -(d-a) \pm sqrt(disc)) / (2 * c).
			"""
			self.disc = pow(self.d - self.a, 2) + 4 * self.b * self.c			
										# debugging
			print ("### last linfractrans, a, b, c, d, and discriminant")
			print ("###", self.a, self.b, self.c, self.d, self.disc)
										# end debugging

			self.gcd1 = euclid_alg (self.d - self.a, self.c)
			self.gcd  = euclid_alg (self.gcd1, self.b)
										# debugging
			#	print ("### gcd: ", self.gcd)
										# end debugging
			if self.gcd != 0:
				self.d_minus_a = (self.d - self.a) // self.gcd
				self.b = self.b // self.gcd
				self.c = self.c // self.gcd
				self.disc = self.disc // (self.gcd * self.gcd)


													# break into square x sq_free
			self.decomp = square_part (self.disc)
			self.coeff_radix = self.decomp[0]
			self.radix = self.decomp[1]
			self.fixed_pt1 = surd.Surd (-(self.d_minus_a), self.coeff_radix, \
										 	2 * self.c, self.radix)
			self.fixed_pt2 = surd.Surd (-(self.d_minus_a), -self.coeff_radix,\
											2 * self.c, self.radix)
			return [self.fixed_pt1, self.fixed_pt2]


#----- end of class -------


