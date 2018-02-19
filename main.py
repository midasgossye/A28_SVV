# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 14:07:00 2018

@author: A28
"""
# Imports
from math import *

# Global variables
C_a = 0.515                 # Chord length aileron [m]
l_a = 2.691                 # Span of the aileron [m]
x_1 = 0.174                 # x-location of hinge 1 [m]
x_2 = 1.051                 # x-location of hinge 2 [m]
x_3 = 2.512                 # x-location of hinge 3 [m]
x_a = 0.300                 # Distance between actuator 1 and 2 [m]
h = 0.248                   # Aileron height [m]
t_sk = 1.1 * 10 ** (-3)     # Skin thickness [m]
t_sp = 2.2 * 10 ** (-3)     # Spar thickness [m]
t_st = 1.2 * 10 ** (-3)     # Thickness of stiffener [m]
h_st = 1.5 * 10 ** (-2)     # Height of stiffener [m]
w_st = 3.0 * 10 ** (-2)     # Width of stiffener [m]
n_st = 11                   # Number of stiffeners [-]
d_1 = 10.34 * 10 ** (-2)    # Vertical displacement hinge 1 [m]
d_3 = 20.66 * 10 ** (-2)    # Vertical displacement hinge 3 [m]
theta = 25                  # Maximum upward deflection [deg]
P = 20.6 * 10 ** 3          # Load in actuator 2 [N]
q = 1.00 * 10 ** 3          # Net aerodynamic load [N/m]





# calculating the cross section of components of the aileron
def cross_section(ha, ca, tskin, tspar, stiffener_amount, w_stiffener, t_stiffener, h_stiffener):
    # C shape
    cshape = 0.5 * pi * ((ha / 2) ** 2) - 0.5 * pi * ((ha - (2 * tskin)) / 2) ** 2
    # spar
    spar = tspar * (ha - (2 * tskin))
    # triangle
    triangle = 0.5 * ha * (ca - 0.5 * ha) - 0.5 * (ha - 2 * tskin) * (ca - 0.5 * ha - tskin)
    # stiffeners
    stiffeners = stiffener_amount * (w_stiffener * t_stiffener + (h_stiffener - t_stiffener) * t_stiffener)
    return cshape, spar, triangle, stiffeners  # unit: m^2


# calculating the enclosed cross sectional area
def enc_area(ha, ca, tskin, tspar):
    A_1 = 0.5 * pi * ((ha - (1 * tskin)) / 2) ** 2  # Circular section enclosed area
    A_2 = 0.5 * (ha - 1 * tskin) * (ca - 0.5 * ha - tskin)  # Triangular section enclosed area
    return A_1, A_2


# functions
def Inertia(h, t_sk, n_st):
    circle_perim = 0.5 * pi * (0.5 * h - t_sk)
    slope = sqrt((0.5 * h - t_sk) ** 2 + (C_a - 0.5 * h - t_sk) ** 2)
    total_perimeter = circle_perim + slope
    slope_angle = atan(0.5 * h / (C_a - 0.5 * h)) - radians(180)
    emptyspacetrailingedge = (h_st * cos(slope_angle)) / (sin(slope_angle))
    print total_perimeter
    print emptyspacetrailingedge
    spacing = (total_perimeter - emptyspacetrailingedge) / ((n_st + 1) / 2)

    for i in xrange(7):
        local_spacing = i * spacing
        if local_spacing < circle_perim:
            angle = (local_spacing / circle_perim) * radians(90)
            z_coordinate = -0.5 * h - t_sk + cos(angle) * (0.5 * h - t_sk)
            y_coordinate = sin(angle) * (0.5 * h - t_sk)
            rot_angle = angle + radians(90)

        else:
            rot_angle = atan(0.5 * h / (C_a - 0.5 * h)) - radians(180)
            z_coordinate = -((h - t_sk) / 2) - (local_spacing - circle_perim) * cos(atan(0.5 * h / (C_a - 0.5 * h)))
            y_coordinate = (slope - (local_spacing - circle_perim)) * sin(atan(0.5 * h / (C_a - 0.5 * h)))

        print "Stif.", i, "\t z:", z_coordinate, "\t y:", y_coordinate, "\t angle:", degrees(rot_angle)
    print "the max height of the possible stiffener location is:", slope * sin(atan(0.5 * h / (C_a - 0.5 * h)))
    print "max possible height should be:", h / 2 - t_sk


# testing
Inertia(h, t_sk, n_st)
