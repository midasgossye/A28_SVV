# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 14:07:00 2018

@author: A28
"""
# Imports
from math import *
import unittest
import scipy.integrate as integrate
import numpy as np
from matplotlib import pyplot as plt

from int_stress_and_defl import *
import internal_shear_and_moment as intsm
import total_shear_calc as totshear

# Global variables
C_a = 0.515  # Chord length aileron [m]
l_a = 2.691  # Span of the aileron [m]
x_1 = 0.174  # x-location of hinge 1 [m]
x_2 = 1.051  # x-location of hinge 2 [m]
x_3 = 2.512  # x-location of hinge 3 [m]
x_a = 0.300  # Distance between actuator 1 and 2 [m]
h = 0.248  # Aileron height [m]
t_sk = 1.1 * 10 ** (-3)  # Skin thickness [m]
t_sp = 2.2 * 10 ** (-3)  # Spar thickness [m]
t_st = 1.2 * 10 ** (-3)  # Thickness of stiffener [m]
h_st = 1.5 * 10 ** (-2)  # Height of stiffener [m]
w_st = 3.0 * 10 ** (-2)  # Width of stiffener [m]
n_st = 11  # Number of stiffeners [-]
d_1 = 10.34 * 10 ** (-2)  # Vertical displacement hinge 1 [m]
d_3 = 20.66 * 10 ** (-2)  # Vertical displacement hinge 3 [m]
theta = 25  # Maximum upward deflection [deg]
P = 20.6 * 10 ** 3  # Load in actuator 2 [N]
q = 1.00 * 10 ** 3  # Net aerodynamic load [N/m]
G = 28 * 10 ** 9  # Shear modulus in Pa (28 GPa, source: http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=ma2024t3)
n = 270  # sections to be analysed

datax = [0.175, 0.902, 1.052, 1.202, 2.513]
dataptx = []
for i in xrange(len(datax)):
    dataptx.append(int(datax[i] / l_a * n))
print dataptx, "<--- positions are from 1000 sections"

E = 71 * 10 ** 9


# functions

# calculating the cross section of components of the aileron
# input  height aileron ha, chord length aileron ca, skin thickness tskin, spar thickness tspar,
# stiffener_amount, width stiffener w_stiffener, thickness stiffener t_stiffener, height stiffener h_stiffener
# return cshape, spar, triangle, stiffeners  # unit: m^2
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
# input height aileron ha, chord length aileron ca, thickness tskin
# returns Circular section enclosed area A_1, Triangular section enclosed area A_2 #areas m^2
def enc_area(ha, ca, tskin):
    A_1 = 0.5 * pi * ((ha - (1 * tskin)) / 2) ** 2  # Circular section enclosed area
    A_2 = 0.5 * (ha - 1 * tskin) * (ca - 0.5 * ha - tskin)  # Triangular section enclosed area
    return A_1, A_2


# inertia

# returns stiffener z,y locations and rotation
# return z_y_angle_coords  # [(stringer0 z,y,rot),(stringer1 z,y,rot)] m,m,rad
def stif_loc(h, t_sk, n_st):
    circle_perim = 0.5 * pi * (0.5 * h - t_sk)
    total_perimeter = circle_perim + sqrt((0.5 * h - t_sk) ** 2 + (C_a - 0.5 * h - t_sk) ** 2)  # m

    spacing = total_perimeter / ((n_st + 1) / 2)
    z_y_angle_coords = []
    for i in xrange(6):
        local_spacing = i * spacing
        if local_spacing < circle_perim:
            angle = (local_spacing / circle_perim) * radians(90)
            z_coordinate = -1 * (0.5 * h - (0.5 * h - t_sk + cos(angle) * (0.5 * h - t_sk)))
            y_coordinate = sin(angle) * (0.5 * h - t_sk)
            rot_angle = angle + radians(90)

        else:
            rot_angle = atan(0.5 * h / (C_a - 0.5 * h)) - radians(180)
            z_coordinate = (-1) * (local_spacing - circle_perim) * cos(atan(0.5 * h / (C_a - 0.5 * h)))
            y_coordinate = h / 2 - (local_spacing - circle_perim) * sin(atan(0.5 * h / (C_a - 0.5 * h)))

        apnd_itm = (z_coordinate, y_coordinate, rot_angle)
        z_y_angle_coords.append(apnd_itm)
        if i > 0:
            apnd_itm = (z_coordinate, -y_coordinate, -rot_angle)
            z_y_angle_coords.append(apnd_itm)

        # print "Stif.", i, "\t z:", z_coordinate, "\t y:", y_coordinate, "\t angle:", degrees(rot_angle)

    return z_y_angle_coords  # [(stringer0 z,y,rot),(stringer1 z,y,rot), ...]


# function to calculate torsional constant
# input height h, thickness skin t_sk, chord length aileron C_a
# return J  # torsional constant
def torsional_constant(h, t_sk, C_a):
    midcircle_perim = pi * (0.5 * h - 0.5 * t_sk)  # wall mid line perimeter circular
    midtriangle_perim = 2 * (
            sqrt((0.5 * h - t_sk) ** 2 + (C_a - 0.5 * h - t_sk) ** 2) - 0.5 * t_sk)  # wall mid line perimeter triangle
    p = midcircle_perim + midtriangle_perim  # wall mid line perimeter
    AeC, AeT = enc_area(h, C_a, t_sk)  # enclosed area of circular part and triangle part
    Ae = AeC + AeT  # total enclosed area
    J = (4 * Ae ** 2 * t_sk) / p

    return J  # torsional constant


# function for transforming axis to a rotated version
# input MMOI I_zz, I_yy, I_zy, and rotation angle rot_angle
# outputs new rotated I_uu, I_vv, I_uv
def axis_transformation(I_zz, I_yy, I_zy, rot_angle):
    # Axis transformation for rotated axis system used for Inertia calculations
    I_uu = (I_zz + I_yy) * 0.5 + (I_zz - I_yy) * 0.5 * cos(2 * rot_angle) - I_zy * sin(2 * rot_angle)
    I_vv = (I_zz + I_yy) * 0.5 - (I_zz - I_yy) * 0.5 * cos(2 * rot_angle) + I_zy * sin(2 * rot_angle)
    I_uv = (I_zz - I_yy) * 0.5 * sin(2 * rot_angle) + I_zy * cos(2 * rot_angle)
    return I_uu, I_vv, I_uv


# function to calculate MMOI
# input stiffener data as z_y_angle_coords, thickness skin t_st, height stiffeners h_st, width stiffener w_st,
#  thickness spar t_sp,aileron height  h, maximum upward deflection theta
# returns IZZ in body ref:  TOT_I_zz_br, IYY body ref:  TOT_I_yy_br, IZZ: TOT_I_zz, IYY: TOT_I_yy, IZY: TOT_I_zy
def moment_of_inertia(z_y_angle_coords, t_st, h_st, w_st, t_sp, h, theta):
    #
    # Calculate Inertias for simple beam axis system
    #   |        
    #   |        ^ (y)
    # -------  <--| (z)

    # === Determine base and height values of inv-T beam rectangles
    b_1 = w_st
    h_1 = t_st
    b_2 = t_st
    h_2 = h_st - t_st
    # ===

    # === Calculate individual I_zz and I_yy and sum steiner term
    I_zz_1 = (b_1 * (h_1 ** 3)) / 12 + b_1 * h_1 * ((t_st * 0.5) ** 2)
    I_yy_1 = ((b_1 ** 3) * h_1) / 12

    I_yy_2 = ((b_2 ** 3) * h_2) / 12
    I_zz_2 = (b_2 * (h_2 ** 3) / 12) + b_2 * h_2 * ((h_2 * 0.5 + t_st) ** 2)
    # ===

    # === BASE INERTIAS AND AREA FOR INVERSE-T BEAM
    I_zz = I_zz_1 + I_zz_2
    I_yy = I_yy_1 + I_yy_2
    I_zy = 0
    A_st = w_st * t_st + t_st * (h_st - t_st)
    # ===

    TOT_I_zz_br = 0
    TOT_I_yy_br = 0
    TOT_I_zy_br = 0
    for coords in z_y_angle_coords:
        z_coord, y_coord, rot_angle = coords  # Get z,y and rotation angle for each stiffener
        stiff_I_zz, stiff_I_yy, stiff_I_zy = axis_transformation(I_zz, I_yy, I_zy,
                                                                 rot_angle)  # perform inertia axis angle transformation
        I_zz_body_ref = stiff_I_zz + A_st * (y_coord ** 2)  # Apply parallel axis theorem
        I_yy_body_ref = stiff_I_yy + A_st * (z_coord ** 2)  # Apply parallel axis theorem
        I_zy_body_ref = stiff_I_zy + A_st * y_coord * z_coord  # Apply parallel axis theorem

        # === SUM ALL STIFFENER MOMENTS OF INERTIA's W.R.T. BODY REFERENCE SYSTEM
        # NOTE: TOTAL I_zy inertia should be zero, because total cross-section has an axis of symmetry
        #       If calculated TOTAL I_zy is NOT equal to zero, there is an error in the computation
        TOT_I_zz_br += I_zz_body_ref
        TOT_I_yy_br += I_yy_body_ref
        TOT_I_zy_br += I_zy_body_ref  # Should be zero, if not => check values!

    # === Semi_circle Moment of inertia:

    I_zz_s_circ = integrate.quad(lambda x: t_sk * ((0.5 * h * sin(x)) ** 2) * 0.5 * h, -pi / 2, pi / 2)[0]

    I_yy_s_circ = I_zz_s_circ
    TOT_I_zz_br += I_zz_s_circ
    TOT_I_yy_br += I_yy_s_circ

    # ===

    # === Triangle skin moment of inertia
    a = sqrt((0.5 * h - t_sk) ** 2 + (C_a - 0.5 * h - t_sk) ** 2)
    angle = atan(0.5 * h / (C_a - 0.5 * h))

    I_zz_t = ((a ** 3 * t_sk * (sin(angle)) ** 2) / 12 + a * t_sk * (0.25 * (h - t_sk)) ** 2) * 2
    # print angle, I_zz_t
    I_yy_t = 2 * ((a ** 3 * t_sk * (cos(angle)) ** 2) / 12) + 2 * a * t_sk * ((C_a - 0.5 * h - t_sk) * 0.5) ** 2

    TOT_I_zz_br += I_zz_t
    TOT_I_yy_br += I_yy_t
    # ===

    # === Spar Moment of Inertia
    I_zz_spar = (t_sp * (h - 2 * t_sk) ** 3) / 12
    # I_yy of spar is negligible since you have a t^3 term if using the thin walled approx.
    # NOTE: t/h << 1

    TOT_I_zz_br += I_zz_spar
    # ===

    # === Transform Inertias from Body Reference system to Main Reference system
    TOT_I_zz, TOT_I_yy, TOT_I_zy = axis_transformation(TOT_I_zz_br, TOT_I_yy_br, TOT_I_zy_br, theta)

    # ===

    # Returns I_zz and I_yy in our OWN DEFINED BODY REFERENCE SYSTEM, followed by the I_zz, I_yy and I_zy in the main reference system
    # NOTE: All reported values are in m^4
    return TOT_I_zz_br, TOT_I_yy_br, TOT_I_zz, TOT_I_yy, TOT_I_zy


# calculates boom area
# input stiffener loaction, thickness stiffeners, height stiffeners,width stiffeners, thickness spar, height
# outputs boom_arear of booms in a list b_i_arr, b_i_spar
def boom_area_calc(stif_loc, t_st, h_st, w_st, t_sp, h):
    A_st = w_st * t_st + (h_st - t_st) * t_st

    circle_perim = 0.5 * pi * (0.5 * h - t_sk)
    total_perimeter = circle_perim + sqrt((0.5 * h - t_sk) ** 2 + (C_a - 0.5 * h - t_sk) ** 2)  # m

    spacing = total_perimeter / ((n_st + 1) / 2)
    B_i_arr = []
    for i in xrange(n_st):
        if i == 0:
            sigma_ratio = (stif_loc[i][1]) / (stif_loc[i + 1][1])
            B_i = A_st + ((t_sk * spacing) / 6) * (2 + sigma_ratio)
            B_i_arr.append(B_i)
        elif i == (n_st - 2):
            sigma_ratio = (stif_loc[i - 2][1]) / (stif_loc[i][1])
            B_i = A_st + ((t_sk * spacing) / 6) * (2 + sigma_ratio)
            B_i_arr.append(B_i)
            B_i_arr.append(B_i)
        elif i < (n_st - 2):
            sigma_ratio = (stif_loc[i][1]) / (stif_loc[i + 2][1])
            B_i = A_st + ((t_sk * spacing) / 6) * (2 + sigma_ratio)
            B_i_arr.append(B_i)

    sigma_ratio = -1
    B_spar_end = ((t_sp * (h - 2 * t_sk)) / 6) * (2 + sigma_ratio)

    # returns an array with all stiffener boom area's and the value of the spar end_cap area
    # 0-th boom is boom at LE, 1st 
    return B_i_arr, B_spar_end


def plot_numerical_bending(Inertia_bend, loc_data, moment_data, E):
    # This function can calculate the numerical bending defelction both in y and z directions
    # Parameters: Inertia refernced to bending direction, location data in x-direction, Moments [Nm], Young's Modulus [Pa]

    least_error = 100
    best_value = (0.0, 0.0)

    # The following for-loops will try to find the best fitting integration constants (c_1 and c_2) to match the deflection to the boundary conditions
    for c_1 in np.arange(-0.1, -0.1, 0.001):
        for c_2 in np.arange(-0.1, -0.1, 0.001):
            # === Initialise first integration
            Int = 0
            Int_arr = np.array([])
            # === 
            for i in xrange(270):
                x_1 = loc_data[i]  #
                x_2 = loc_data[i + 1]
                m_x_1 = (-1 * moment_data[i]) / (E * Inertia_bend)
                m_x_2 = (-1 * moment_data[i + 1]) / (E * Inertia_bend)

                dx = x_2 - x_1

                Int += dx * (m_x_1 + m_x_2) / 2 + c_1
                Int_arr = np.append(Int_arr, Int)

            Int_2 = 0
            Int_2_arr = np.array([])
            for i in xrange(269):
                x_1 = loc_data[i]
                x_2 = loc_data[i + 1]
                y_1 = Int_arr[i]
                y_2 = Int_arr[i + 1]

                dx = x_2 - x_1

                Int_2 += dx * (y_1 + y_2) / 2 + c_2
                Int_2_arr = np.append(Int_2_arr, Int_2)

            print "c1:", c_1, "\tc2:", c_2
            print Int_2_arr[16], "\t", Int_2_arr[105]
            max_defl = max(Int_2_arr)

            error = abs(Int_2_arr[16] + max_defl + (0.1034 / 2.54) / 2) + abs(
                Int_2_arr[251] + max_defl + (0.2066 / 2.54)) / 2  # + abs(Int_2_arr[105])/3
            # + abs(Int_2_arr[105])/3
            if error < least_error:
                best_value = (c_1, c_2)
                least_error = error

    print best_value
    c_1 = best_value[0]

    c_2 = best_value[1]
    c_1 = 0.002
    c_2 = -0.003
    # c_1 = 0.0011
    # c_2 = -0.001
    Int = 0
    Int_arr = np.array([])
    for i in xrange(270):
        x_1 = loc_data[i]
        x_2 = loc_data[i + 1]
        m_x_1 = (-1 * moment_data[i]) / (float(E * I_bend))
        m_x_2 = (-1 * moment_data[i + 1]) / (float(E * I_bend))
        print I_bend
        dx = x_2 - x_1

        Int += dx * (m_x_1 + m_x_2) / 2 + c_1
        Int_arr = np.append(Int_arr, Int)

    Int_2 = 0
    Int_2_arr = np.array([])
    for i in xrange(269):
        x_1 = x_coor[i]
        x_2 = x_coor[i + 1]
        y_1 = Int_arr[i]
        y_2 = Int_arr[i + 1]

        dx = x_2 - x_1

        Int_2 += dx * (y_1 + y_2) / 2 + c_2
        Int_2_arr = np.append(Int_2_arr, Int_2)
    max_defl = max(Int_2_arr)
    print "c1:", c_1, "\tc2:", c_2
    print Int_2_arr[16], "\t", Int_2_arr[105]
    plt.grid()
    Int_2_arr_inv = np.array([])
    for i in xrange(len(Int_2_arr) - 1, -1, -1):
        Int_2_arr_inv = np.append(Int_2_arr_inv, Int_2_arr[i])
    print max_defl
    plot_arr = Int_2_arr_inv + max_defl + 0.088962 - 0.00077
    plt.plot(x_coor[:269], -Int_2_arr_inv + max_defl + 0.22479 - 0.359795)  # +1*max_defl+0.098)
    plt.xlabel("x-coordinate [m]")
    plt.ylabel("z-deflection [m]")
    plt.show()

    return x_coor[:269], plot_arr


# B_i_arr[0]
print "Moments: (I_z'z', I_y'y', I_zz, I_yy, I_zy) All in m^4"
print moment_of_inertia(stif_loc(h, t_sk, n_st), t_st, h_st, w_st, t_sp, h, theta)

Moment_data = np.genfromtxt("M_y.txt")
x_coor = Moment_data[:, 0]
I_bend = 6.385385647322895e-05
M_x = Moment_data[:, 1]

#TODO: uncomment at end
# x_coor, plot_arr = plot_numerical_bending(I_bend, x_coor, M_x, E)
# file_n = open("defl_data_z.txt", "w")
#
# for i in xrange(len(plot_arr)):
#     file_n.write(str(plot_arr[i]))
#     file_n.write("\n")
# file_n.close()

# print "Moments: (I_z'z', I_y'y', I_zz, I_yy, I_zy) All in m^4"
# print moment_of_inertia(stif_loc(h, t_sk, n_st), t_st, h_st, w_st, t_sp, h, theta)


# main
stif_data = stif_loc(h, t_sk, n_st)  # initialize stiffener properties
# Boom idealization
b_r = []  # list for the resultant boom areas

# calculating the stiffeners' total boom area
b_sp = []

b_r, b_sp = boom_area_calc(stif_data, t_st, h_st, w_st, t_sp,
                           h)  # b_r is the list off stiffener boom areas, b_sp for spar(single value in m^2)

J = torsional_constant(h, t_sk, C_a)
# crosssection
A = sum(cross_section(h, C_a, t_sk, t_sp, n_st, w_st, t_st, h_st))  # A is sum of cross section

I_zz_br, I_yy_br, I_zz, I_yy, I_zy = moment_of_inertia(stif_data, t_st, h_st, w_st, t_sp, h,
                                                       theta)  # values for the MMOI

enclosed = sum(enc_area(h, C_a, t_sk))  # enclosed area size

model = []  # whole model
section_length = l_a / n
verifdata = []
qribdata = []


def iteration(section_number):
    x_start = section_number * section_length
    mid = x_start + section_length / 2
    M, V_y, V_z, V_ypr, V_zpr = intsm.internal(mid, I_zz)

    stif_data = stif_loc(h, t_sk, n_st)
    bir, bisp = boom_area_calc(stif_data, t_st, h_st, w_st, t_sp, h)
    totshearvalue, qrib = totshear.totalshear(stif_data, V_zpr, V_ypr, bir, bisp, I_zz_br, I_yy_br)
    verifdata.append(totshearvalue[7])
    verifdata.append(totshearvalue[9])
    verifdata.append(totshearvalue[0])
    verifdata.append(totshearvalue[5])
    # qribdata.append(qrib[0])
    # qribdata.append(qrib[1])
    # print "section: ",section_number,"at x: ", mid
    # print totshearvalue
    # print qrib
    return


for i in xrange(len(dataptx)):
    iteration(dataptx[i])
verifdata += qribdata
print verifdata


# print I_zz_br, I_yy_br

# x_start = 0 * section_length
# mid = x_start + section_length / 2
# M, V_y, V_z, V_ypr, V_zpr = intsm.internal(mid, I_zz)
#
# stif_data = stif_loc(h, t_sk, n_st)
# bir, bisp = boom_area_calc(stif_data, t_st, h_st, w_st, t_sp, h)
# totshearvalue = totshear.totalshear(stif_data, V_zpr, V_ypr, bir, bisp, I_zz_br, I_yy_br)
# print totshearvalue

x_start = n / 2 * section_length
mid = x_start + section_length / 2
M, V_y, V_z, V_ypr, V_zpr = intsm.internal(mid, I_zz)

stif_data = stif_loc(h, t_sk, n_st)
bir, bisp = boom_area_calc(stif_data, t_st, h_st, w_st, t_sp, h)
totshearvalue = totshear.totalshear(stif_data, V_zpr, V_ypr, bir, bisp, I_zz_br, I_yy_br)
# print totshearvalue

# for y in xrange(n):
#     model.append(iteration(y))


# internal stress and deflection


# test
# print "stiff location print:", stif_loc(h, t_sk, n_st)
# print "torsional constant", torsional_constant(h, t_sk, C_a)
# testunits for unittests
