# -*- coding: utf-8 -*-
"""
@author: A28
"""
# Imports
from math import *


# TODO: blind built without specifications from other modules, change if inertia function is finalized
# normal stress due to bending
def norm_strs(M_z, I_z_z, y):
    sigma_x = M_z / I_z_z * y
    return sigma_x
