# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 14:07:00 2018

@author: A28
"""
# Imports
from math import *

C_a = 0.515  # Chord length aileron [m]


# functions
def Inertia(None):
    return None


# calculating the cross section
def cross_section(ha, ca, tskin, tspar):
    # C shape
    cshape = 0.5 * math.pi * ((ha / 2) ** 2) - 0.5 * math.pi * ((ha - (2 * tskin) / 2) ** 2)
    # spar
    spar = tspar * (ha - (2 * tskin))
    # triangle
    triangle = 0.5 * ha * (ca - 0.5 * ha) - 0.5 * (ha - 2 * tskin) * (ca - 0.5 * ha - tskin)
    return cshape + spar + triangle
