# -*- coding: utf-8 -*-
"""
@author: A28
"""
# Imports
from math import *


# normal stress due to bending
# returns the normal stress sigma_x
def norm_strs(M_z, I_z_z, y):
    sigma_x = M_z / I_z_z * y
    return sigma_x


# deflection function
# returns v prime prime
def defl(M_x, I_z_z, E):
    v = -1 * M_x / (E * I_z_z)
    return v

#TODO numerical implementation of a contour integration
# rate of twist
# returns d(theta)/dz
def ROT(A, G, t, q_s):
    theta = 1 / (2 * A)
    return dTheta

#maximum deflection
def maxdefl(dTheta, v):
    leading = dTheta
    trailing = dTheta + v
    return leading, trailing


