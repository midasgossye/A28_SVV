# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 14:07:00 2018

@author: A28
"""
# Imports
from math import *

# Global variables
C_a = 0.515         # Chord length aileron [m]
l_a = 2.691         # Span of the aileron [m]
x_1 = 0.174         # x-location of hinge 1 [m]
x_2 = 1.051         # x-location of hinge 2 [m]
x_3 = 2.512         # x-location of hinge 3 [m]
x_a = 0.300         # Distance between actuator 1 and 2 [m]
h   = 0.248         # Aileron height [m]
t_sk = 1.1*10**(-3) # Skin thickness [m]
t_sp = 2.2*10**(-3) # Spar thickness [m]
t_st = 1.2*10**(-3) # Thickness of stiffener [m]
h_st = 1.5*10**(-2) # Height of stiffener [m]
w_st = 3.0*10**(-2) # Width of stiffener [m]
n_st = 11           # Number of stiffeners [-]
d_1 = 10.34*10**(-2)# Vertical displacement hinge 1 [m]
d_3 = 20.66*10**(-2)# Vertical displacement hinge 3 [m]
theta = 25          # Maximum upward deflection [deg]
P = 20.6*10**3      # Load in actuator 2 [N]
q = 1.00*10**3      # Net aerodynamic load [N/m]
        
# functions
def Inertia():
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
    
ail_1 = Aileron()


