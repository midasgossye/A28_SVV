# -*- coding: utf-8 -*-
"""
@author: A28
"""
# Imports
from math import *
import scipy.integrate as intg


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


# TODO numerical implementation of a contour integration(done with numpy.quad): check with test data
# rate of twist
# returns d(theta)/dz
def ROT(A, G, t, q_s):
    s = 11
    dTheta = 1 / (2 * A) * intg.quad(q_s / G / t, 0, s)
    return dTheta


# maximum deflection
def maxdefl(dTheta, v, LEdist, TEdist, z):
    leading = dTheta * z * LEdist - v
    trailing = dTheta * z * TEdist + v
    return leading, trailing
