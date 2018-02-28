import numpy as np
from math import *
from scipy.linalg import solve
import sympy as sp

# ----------------------------------------------------------------
# Shear calculations


# ----------------------------------------------------------------
# ----DEFINE CONSTANTS--------------------------------------------
# ----------------------------------------------------------------

# ----Given constants---------------------------------------------

P = 20.6 * 10 ** (3)  # N
q = 1.00 * 10 ** (3)  # N/m
C_a = 0.515  # m
l_a = 2.691  # m
x_1 = 0.174  # m
x_2 = 1.051  # m
x_3 = 2.512  # m
x_a = 0.3  # m
h = 0.248  # m
theta = 25 * pi / 180  # radians
f_s = 0.4  # assumed form factor for cross section
A = 0.0021108225743942223  # cross sectional area [m2]
G = 28 * 10 ** 9  # shear modulus aluminum 2014-T3 [Pa]
d_1 = 0.134  # m
d_3 = 0.2066  # m
E = 73.1 * 10 ** 9  # Pa
I_zprimezprime = 1.437615789242078e-05  # m**4
I_yprimeyprime = 6.473600927242316e-05  # m**4
S_zprime =
S_yprime =

# ----Calculating distances---------------------------------------

y_hingeline = h / 2 * sin((pi / 2) - theta) - h / 2 * sin(theta)  # m
z_hingeline = (0.25 * C_a - h / 2) * cos(theta)  # m
z_star = h / 2 * cos(theta)  # m

# ----Calculating basic shear flow q_b----------------------------

# ----Top Circular part of cross-section------------------------------


# Starting at the cut between boom 0 and 1 and moving along the top of the cross section in the direction of spar 1
q_b1 = [0]
q_open_b, 1 = -S_zprime / I_zprimezprime * sum(B_i_arr[1] * z_y_angle_coords[1][0]) - S_yprime / I_yprimeyprime * sum(
    B_i_arr[1] * z_y_angle_coords[1][1])
q_open_b, 1 = q_b1[0] + q_open_b, 1
q_b1.append(q_open_b, 1)
q_open_spar1_left = -S_zprime / I_zprimezprime * sum(
    B_spar_end * z_y_angle_coords[1][0]) - S_yprime / I_yprimeyprime * sum(B_spar_end * z_y_angle_coords[1][1])
# coordinaten hierboven kloppen niet.
q_open_spar1_left = q_b1[1] + q_open_spar1_left
q_b1.append(q_open_spar1_left)

# Starting at the trailing edge of the cross-section and going over the top of the cross-section
q_b2 = [0]
q_open_b, 9 = -S_zprime / I_zprimezprime * sum(B_i_arr[9] * z_y_angle_coords[9][0]) - S_yprime / I_yprimeyprime * sum(
    B_i_arr[9] * z_y_angle_coords[9][1])
q_open_b, 9 = q_b2[0] + q_open_b, 9
q_b2.append(q_open_b, 9)
q_open_b, 7 = -S_zprime / I_zprimezprime * sum(B_i_arr[7] * z_y_angle_coords[7][0]) - S_yprime / I_yprimeyprime * sum(
    B_i_arr[7] * z_y_angle_coords[7][1])
q_open_b, 7 = q_b2[1] + q_open_b, 7
q_b2.append(q_open_b, 7)
q_open_b, 5 = -S_zprime / I_zprimezprime * sum(B_i_arr[5] * z_y_angle_coords[5][0]) - S_yprime / I_yprimeyprime * sum(
    B_i_arr[5] * z_y_angle_coords[5][1])
q_open_b, 5 = q_b2[2] + q_open_b, 5
q_b2.append(q_open_b, 5)
q_open_b, 3 = -S_zprime / I_zprimezprime * sum(B_i_arr[3] * z_y_angle_coords[3][0]) - S_yprime / I_yprimeyprime * sum(
    B_i_arr[3] * z_y_angle_coords[3][1])
q_open_b, 3 = q_b2[3] + q_open_b, 3
q_b2.append(q_open_b, 3)

# total shear flow at spar 1
q_b_atspar1 = q_b1[2] + q_b2[4]

# shear flow in the rib
q_b3 = q_b_atspar1 - S_zprime / I_zprimezprime * sum(
    B_spar_end * z_y_angle_coords[1][0]) - S_yprime / I_yprimeyprime * sum(B_spar_end * z_y_angle_coords[1][1])

# Starting at the trailing edge of the cross-section and going along the bottom side of the cross-section
q_b4 = [0]
q_open_b, 10 = -S_zprime / I_zprimezprime * sum(
    B_i_arr[10] * z_y_angle_coords[10][0]) - S_yprime / I_yprimeyprime * sum(B_i_arr[10] * z_y_angle_coords[10][1])
q_open_b, 10 = q_b4[0] + q_open_b, 10
q_b2.append(q_open_b, 10)
q_open_b, 8 = -S_zprime / I_zprimezprime * sum(B_i_arr[8] * z_y_angle_coords[8][0]) - S_yprime / I_yprimeyprime * sum(
    B_i_arr[8] * z_y_angle_coords[8][1])
q_open_b, 8 = q_b4[1] + q_open_b, 8
q_b2.append(q_open_b, 8)
q_open_b, 6 = -S_zprime / I_zprimezprime * sum(B_i_arr[6] * z_y_angle_coords[6][0]) - S_yprime / I_yprimeyprime * sum(
    B_i_arr[6] * z_y_angle_coords[6][1])
q_open_b, 6 = q_b4[2] + q_open_b, 6
q_b2.append(q_open_b, 6)
q_open_b, 4 = -S_zprime / I_zprimezprime * sum(B_i_arr[4] * z_y_angle_coords[4][0]) - S_yprime / I_yprimeyprime * sum(
    B_i_arr[4] * z_y_angle_coords[4][1])
q_open_b, 4 = q_b4[3] + q_open_b, 4
q_b2.append(q_open_b, 4)

# total shear flow at spar 2
q_b_atspar2 = q_b3 + q_b4[4]

# starting at spar 2 and moving along the bottom circular part of the cross-section
q_b5 = []
q_open_afterspar2 = -S_zprime / I_zprimezprime * sum(
    B_spar_end * z_y_angle_coords[4][0]) - S_yprime / I_yprimeyprime * sum(B_spar_end * z_y_angle_coords[4][1])
q_open_afterspar2 = q_b_atspar2 + q_open_afterspar2
q_b5.append(q_open_afterspar2)
q_open_b, 2 = -S_zprime / I_zprimezprime * sum(B_i_arr[2] * z_y_angle_coords[2][0]) - S_yprime / I_yprimeyprime * sum(
    B_i_arr[2] * z_y_angle_coords[2][1])
q_open_b, 2 = q_b5[0] + q_open_b, 2
q_b5.append(q_open_b, 2)

# make a list of all the shearflows
q_basicshear = [q_b1[0], q_b1[1], q_b1[2], q_b2[0], q_b2[1], q_b2[2], q_b2[3], q_b2[4], q_b3, q_b4[0], q_b4[1], q_b4[2],
                q_b4[3], q_b4[4], q_b5[0], q_b5[1]]
