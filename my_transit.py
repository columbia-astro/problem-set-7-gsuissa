"""
   Some functions of use for Calculating Transit light curves
"""

import math

#
# This implements the planet obscuration function for a uniform
#  (stellar) disk described in Mandel & Agol (2002).
#  http://dx.doi.org/10.1086/345520
# 

def kappa0(p, z):
    return acos((p**2 + z**2 - 1)/(2*p*z))

def kappa1(p, z):
    return acos((1 - p**2 + z**2)/(2*z))

def lambd(p, z):
    if 1 + p < z:
        return 0
    if z <= 1 - p:
        return p**2
    if z <= p - 1:
        return 1
    arg = (4*z**2 - (1 + z**2 - p**2)**2)/4
    if arg >= 0:
        return (kappa0(p, z)*p**2 + kappa1(p, z) - math.sqrt(arg))/pi
    else:
        return 0

def FluxRatio(p, z):
    """
    Compute the ratio of obscured/unobscured flux for a planet transit.
    
    Arguments:
       p - ratio of planet radius to stellar radius
       z - distance between star and planet divided by stellar radius
    Returns: 
       FluxRatio - ratio of obscured to unobscured stellar flux

    """
    return 1 - lambd(p, abs(z))


# The following functions are used in computing the transit light curve
#

def beta(p, r, z):
    return 2.0 * math.acos((z**2 - p**2 + r**2)/(2*z*r))

def delta(p, r, z):
    """
    Compute the ratio of obscured/unobscured flux for a planet transit
       for a given radial ring.
    
    Arguments:
       p - ratio of planet radius to stellar radius
       r - radius of ring (divided by stellar radius)
       z - distance between star and planet divided by stellar radius
    Returns: 
       delta - ratio of obscured to unobscured stellar flux for that radius
    """
    if r >= z + p or r <= z - p:
        return 0
    if r + z <= p:
        return 1
    b = beta(p, r, z)
    return b/(2*math.pi)

