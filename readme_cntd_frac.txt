README for Continued Fraction Implementation in Python

Copyright 2018 Jesse I. Deutsch


  We wish to find the continued fraction expansion of
quadratic expressions over the rationals.  The typical
expression will be a surd, that is 

 a+ b*sqrt(r)
-------------- .
      d

Here a, b, d are rational integers and r is positive integer
greater than one that is not a perfect square.

  A typical session would look like this.  Invoke python 3,
then write

>>> from cntd_frac import *
>>> y = Surd(1, 2, 5, 3)
>>> print (y)
[1+2*rt(3)] / 5
>>> cntd_frac(y)
[0, 1, [8, 3, 34, 3]]
>>> z = Surd(1, 1, 2, 5)
>>> print (z)
[1+1*rt(5)] / 2
>>> cntd_frac(z)
[[1]]

