# -*- coding: utf-8 -*-
"""
@author: A28
"""
# Imports
from math import *


# normal stress due to bending
def norm_strs(M_z, I_z_z, y):
    sigma_x = M_z / I_z_z * y
    return sigma_x


# deflection function
# returns v prime prime
def defl(M_x, I_z_z, E):
    v = -1 * M_x / (E * I_z_z)
    return v

