import numpy as np
from math import *
from scipy.linalg import solve
import sympy as sp
import Manouschka as Mnk

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

# De volgende constanten staan denk ik ook ergens in het grote main bestand?
A_1 = pi * (h / 2) ** 2 / 2  # m^2
A_2 = h * (l_a - (h / 2))  # m^2
t_skin = 0.0011  # m
t_spar = 0.0022  # m
phi = 0.3071  # radians, the angle between the symmetry axis and the top or bottom of the triangle of the cross-section


# --THE FOLLOWING VALUES WILL FOLLOW FROM THE NUMERICAL METHOD
# --S_zprime and S_yprime are the values of the resultant forces along the body axes.
# I_zprimezprime = 1.437615789242078e-05  # m**4
# I_yprimeyprime = 6.473600927242316e-05  # m**4
# S_zprime = 714.647480604  # N
# S_yprime = -1532.56646788  # N
# B_i_arr = np.ones(15) * 0.00125  # m**2


# print totalshear


def totalshear(z_y_angle_coords, S_zprime, S_yprime, B_i_arr, B_spar_end, I_zprimezprime, I_yprimeyprime):
    S_zprime = S_zprime[0]
    S_yprime = S_yprime[0]
    lengths = Mnk.lengths()
    # ----Calculating distances---------------------------------------

    y_hingeline = h / 2 * sin((pi / 2) - theta) - h / 2 * sin(theta)  # m
    z_hingeline = (0.25 * C_a - h / 2) * cos(theta)  # m
    z_star = h / 2 * cos(theta)  # m

    # ----------------------------------------------------------------
    # ----CALCULATING BASIC SHEAR FLOW Q_B----------------------------
    # ----------------------------------------------------------------

    # ----Top Circular part of cross-section------------------------------

    # Starting at the cut between boom 0 and 1 and moving along the top of the cross section in the direction of spar 1
    q_b1 = [0]
    q_open_b_1 = -S_zprime / I_zprimezprime * B_i_arr[1] * z_y_angle_coords[1][0] - S_yprime / I_yprimeyprime * B_i_arr[
        1] * z_y_angle_coords[1][1]
    q_open_b_1 = q_b1[0] + q_open_b_1
    q_b1.append(q_open_b_1)

    # Starting at the trailing edge of the cross-section and going over the top of the cross-section
    q_b2 = [0]
    q_open_b_9 = -S_zprime / I_zprimezprime * B_i_arr[9] * z_y_angle_coords[9][0] - S_yprime / I_yprimeyprime * B_i_arr[
        9] * z_y_angle_coords[9][1]
    q_open_b_9 = q_b2[0] + q_open_b_9
    q_b2.append(q_open_b_9)
    q_open_b_7 = -S_zprime / I_zprimezprime * B_i_arr[7] * z_y_angle_coords[7][0] - S_yprime / I_yprimeyprime * B_i_arr[
        7] * z_y_angle_coords[7][1]
    q_open_b_7 = q_b2[1] + q_open_b_7
    q_b2.append(q_open_b_7)
    q_open_b_5 = -S_zprime / I_zprimezprime * B_i_arr[5] * z_y_angle_coords[5][0] - S_yprime / I_yprimeyprime * B_i_arr[
        5] * z_y_angle_coords[5][1]
    q_open_b_5 = q_b2[2] + q_open_b_5
    q_b2.append(q_open_b_5)
    q_open_b_3 = -S_zprime / I_zprimezprime * B_i_arr[3] * z_y_angle_coords[3][0] - S_yprime / I_yprimeyprime * B_i_arr[
        3] * z_y_angle_coords[3][1]
    q_open_b_3 = q_b2[3] + q_open_b_3
    q_b2.append(q_open_b_3)

    # total shear flow at spar 1
    q_b_atspar1 = q_b1[1] + q_b2[4]

    # shear flow in the rib
    q_b3 = q_b_atspar1 - S_zprime / I_zprimezprime * B_spar_end * z_y_angle_coords[1][
        0] - S_yprime / I_yprimeyprime * B_spar_end * z_y_angle_coords[1][1]

    # Starting at the trailing edge of the cross-section and going along the bottom side of the cross-section
    q_b4 = [0]
    q_open_b_10 = -S_zprime / I_zprimezprime * B_i_arr[10] * z_y_angle_coords[10][0] - S_yprime / I_yprimeyprime * \
                  B_i_arr[10] * z_y_angle_coords[10][1]
    q_open_b_10 = q_b4[0] + q_open_b_10
    q_b4.append(q_open_b_10)
    q_open_b_8 = -S_zprime / I_zprimezprime * B_i_arr[8] * z_y_angle_coords[8][0] - S_yprime / I_yprimeyprime * B_i_arr[
        8] * z_y_angle_coords[8][1]
    q_open_b_8 = q_b4[1] + q_open_b_8
    q_b4.append(q_open_b_8)
    q_open_b_6 = -S_zprime / I_zprimezprime * B_i_arr[6] * z_y_angle_coords[6][0] - S_yprime / I_yprimeyprime * B_i_arr[
        6] * z_y_angle_coords[6][1]
    q_open_b_6 = q_b4[2] + q_open_b_6
    q_b4.append(q_open_b_6)
    q_open_b_4 = -S_zprime / I_zprimezprime * B_i_arr[4] * z_y_angle_coords[4][0] - S_yprime / I_yprimeyprime * B_i_arr[
        4] * z_y_angle_coords[4][1]
    q_open_b_4 = q_b4[3] + q_open_b_4
    q_b4.append(q_open_b_4)

    # total shear flow at spar 2
    q_b_atspar2 = q_b3 + q_b4[4]

    # starting at spar 2 and moving along the bottom circular part of the cross-section
    q_b5 = []
    q_open_afterspar2 = -S_zprime / I_zprimezprime * B_spar_end * z_y_angle_coords[4][
        0] - S_yprime / I_yprimeyprime * B_spar_end * z_y_angle_coords[4][1]
    q_open_afterspar2 = q_b_atspar2 + q_open_afterspar2
    q_b5.append(q_open_afterspar2)
    q_open_b_2 = -S_zprime / I_zprimezprime * B_i_arr[2] * z_y_angle_coords[2][0] - S_yprime / I_yprimeyprime * B_i_arr[
        2] * z_y_angle_coords[2][1]
    q_open_b_2 = q_b5[0] + q_open_b_2
    q_b5.append(q_open_b_2)

    # make a list of all the shearflows
    q_basicshear = [q_b1[0], q_b1[1], q_b2[0], q_b2[1], q_b2[2], q_b2[3], q_b2[4], q_b3, q_b4[0], q_b4[1], q_b4[2],
                    q_b4[3],
                    q_b4[4], q_b5[0], q_b5[1]]

    # ----------------------------------------------------------------
    # ----CALCULATING CORRECTIONAL SHEAR FLOW Q_S,0-------------------
    # ----------------------------------------------------------------

    ##q_s_1 = sp.Symbol('q_s_1')
    ##q_s_2 = sp.Symbol('q_s_2')

    # Compute rate of twist for both sections
    # Assume there is a list called lengths, with all the lengths of the sections between the booms of the cross-section, starting from the right along the bottom of the cross-section with the spar as the last input
    # twist1 = 1/(2*A_1*G)* (1/t_skin*((q_basicshear[0]+q_s_1)* lengths[7]+(q_basicshear[1]+q_s_1)* lengths[8]+(q_basicshear[13]+q_s_1)* lengths[5]+(q_basicshear[14]+q_s_1)* lengths[6])+1/t_spar*(q_basicshear[8]+q_s_1-q_s_2)*lengths[13])
    # Hoe zit het hier met de orientatie van de krachten? je moet ze denk ik van elkaar afhalen aangezien de correctional shear flow de andere kant op staat.
    # twist3 = 1/(2*A_2*G)* (1/t_skin*((-q_basicshear[5]+q_s_2)* lengths[12]+ (-q_basicshear[4]+q_s_2)* lengths[11] + (-q_basicshear[3]+q_s_2)* lengths[10] + (-q_basicshear[2]+q_s_2)* lengths[9] + (q_basicshear[12]+q_s_2)*lengths[4] + (q_basicshear[11]+q_s_2)*lengths[3] + (q_basicshear[10]+q_s_2)*lengths[2] + (q_basicshear[9]+q_s_2)*lengths[1]) + (1/t_spar*(-q_basicshear[7]+q_s_2)* lengths[14]) - (1/t_spar*(-q_s_1)* lengths[14]))

    # take moment around midpoint in the spar, clockwise positive
    # Assume there is a list called lengths, with all the lengths of the sectionsbetween the booms on the bottom of the cross-section, starting from the right with lengths[0]
    # All the forces act from the boom they depart from and the shear flow is decomposed in q_z and q_y due to the angle phi.
    # - 2*A_1*q_s_1 - 2*A_2*q_s_2 = q_basicshear[9]*lengths[1]*(cos(phi)*-z_y_angle_coords[10][2]+sin(phi)*-z_y_angle_coords[10][1]) + q_basicshear[10]*lengths[2]*(cos(phi)*-z_y_angle_coords[8][2]+sin(phi)*-z_y_angle_coords[8][1])+q_basicshear[11]*lengths[3]*(cos(phi)*-z_y_angle_coords[6][2]+sin(phi)*-z_y_angle_coords[6][1])+q_basicshear[12]*lengths[4]*(cos(phi)*-z_y_angle_coords[4][2]+sin(phi)*-z_y_angle_coords[4][1]) + q_basicshear[13]*lengths[5]*(h/2) + q_basicshear[14]*lengths[6]*(h/2) + q_basicshear[0]*lengths[7]*(h/2) + q_basicshear[1]*lengths[8]*(h/2) +  q_basicshear[2]*lengths[9]*(-cos(phi)*z_y_angle_coords[3][2]+sin(phi)*z_y_angle_coords[3][1]) + q_basicshear[3]*lengths[10]*(-cos(phi)*z_y_angle_coords[5][2]+sin(phi)*z_y_angle_coords[5][1]) + q_basicshear[4]*lengths[11]*(-cos(phi)*z_y_angle_coords[7][2]+sin(phi)*z_y_angle_coords[7][1]) + q_basicshear[5]*lengths[12]*(-cos(phi)*z_y_angle_coords[9][2]+sin(phi)*z_y_angle_coords[9][1])

    ##Dit is alleen om het in de matrix te zetten, dit zijn de twist equations omgescreven naar matrix formaat
    ## q_s_1*(1/(2*A_1*G)* (1/t_skin*(lengths[7]+ lengths[8]+ lengths[5]+ lengths[6])+ 1/t_spar*lengths[13])- 1/(2*A_2*G)*  (1/t_spar* lengths[14])
    ##- q_s_2*(1/(2*A_1*G)* 1/t_spar*(lengths[13])+ (1/(2*A_2*G)* 1/t_skin*((lengths[12]+  lengths[11] +  lengths[10] +  lengths[9] + lengths[4] + lengths[3] + lengths[2] + lengths[1]) + (1/t_spar*lengths[14]))))
    ##= 1/(2*A_2*G)*
    ##(1/t_skin*(-q_basicshear[5]* lengths[12]+ -q_basicshear[4]* lengths[11] + -q_basicshear[3]* lengths[10] + -q_basicshear[2]* lengths[9] + q_basicshear[12]*lengths[4] + q_basicshear[11]*lengths[3] + q_basicshear[10]*lengths[2] + q_basicshear[9]*lengths[1]) + (1/t_spar*-q_basicshear[7]* lengths[14]))
    ##-1/(2*A_1*G)* ((1/t_skin* (q_basicshear[0]* lengths[7]+q_basicshear[1]* lengths[8]+q_basicshear[13]* lengths[5]+q_basicshear[14]* lengths[6])) +1/t_spar*(q_basicshear[8]*lengths[13]))

    # Simplify input matrix
    correctional1 = (1 / (2 * A_1 * G) * (
            1 / t_skin * (lengths[7] + lengths[8] + lengths[5] + lengths[6]) + 1 / t_spar * lengths[13]) - 1 / (
                             2 * A_2 * G) * (1 / t_spar * lengths[14]))
    correctional2 = 1 / (2 * A_1 * G) * 1 / t_spar * (lengths[13]) + (1 / (2 * A_2 * G) * 1 / t_skin * ((lengths[12] +
                                                                                                         lengths[11] +
                                                                                                         lengths[10] +
                                                                                                         lengths[9] +
                                                                                                         lengths[4] +
                                                                                                         lengths[3] +
                                                                                                         lengths[2] +
                                                                                                         lengths[1]) + (
                                                                                                                1 / t_spar *
                                                                                                                lengths[
                                                                                                                    14])))
    twistrightside = 1 / (2 * A_2 * G) * (1 / t_skin * (
            -q_basicshear[5] * lengths[12] + -q_basicshear[4] * lengths[11] + -q_basicshear[3] * lengths[10] + -
    q_basicshear[2] * lengths[9] + q_basicshear[12] * lengths[4] + q_basicshear[11] * lengths[3] + q_basicshear[10] *
            lengths[2] + q_basicshear[9] * lengths[1]) + (1 / t_spar * -q_basicshear[7] * lengths[14]))
    -1 / (2 * A_1 * G) * ((1 / t_skin * (
            q_basicshear[0] * lengths[7] + q_basicshear[1] * lengths[8] + q_basicshear[13] * lengths[5] + q_basicshear[
        14] * lengths[6])) + 1 / t_spar * (q_basicshear[8] * lengths[13]))
    momentrightside = q_basicshear[9] * lengths[1] * (
            cos(phi) * -z_y_angle_coords[10][2] + sin(phi) * -z_y_angle_coords[10][1]) + q_basicshear[10] * lengths[
                          2] * (cos(phi) * -z_y_angle_coords[8][2] + sin(phi) * -z_y_angle_coords[8][1]) + q_basicshear[
                          11] * lengths[3] * (cos(phi) * -z_y_angle_coords[6][2] + sin(phi) * -z_y_angle_coords[6][1]) + \
                      q_basicshear[12] * lengths[4] * (
                              cos(phi) * -z_y_angle_coords[4][2] + sin(phi) * -z_y_angle_coords[4][1]) + q_basicshear[
                          13] * lengths[5] * (h / 2) + q_basicshear[14] * lengths[6] * (h / 2) + q_basicshear[0] * \
                      lengths[
                          7] * (h / 2) + q_basicshear[1] * lengths[8] * (h / 2) + q_basicshear[2] * lengths[9] * (
                              -cos(phi) * z_y_angle_coords[3][2] + sin(phi) * z_y_angle_coords[3][1]) + q_basicshear[
                          3] * lengths[10] * (-cos(phi) * z_y_angle_coords[5][2] + sin(phi) * z_y_angle_coords[5][1]) + \
                      q_basicshear[4] * lengths[11] * (
                              -cos(phi) * z_y_angle_coords[7][2] + sin(phi) * z_y_angle_coords[7][1]) + q_basicshear[
                          5] * lengths[12] * (-cos(phi) * z_y_angle_coords[9][2] + sin(phi) * z_y_angle_coords[9][1])

    # The equilibrium equations are implemented in 2 matrices, A and B
    # Matrix A consists of q_s_1  q_s_2
    # Matrix B consists of the known values on the right side of the twist and moment equations
    A = np.matrix([[correctional1, correctional2],  # Matrix including the twist and moment equations
                   [- 2 * A_1, - 2 * A_2]])

    B = np.matrix([[twistrightside],  # Matrix inlcuding the right side of the twist and moment equations
                   [momentrightside]])

    correctionalshearflows = solve(A, B)
    # print correctionalshearflows

    # ------------------------------------------------------------
    # ----NOW ALL TOTAL SHEAR FLOWS CAN BE CALCULATED-------------
    # ------------------------------------------------------------

    # Clockwise is positive
    # totalshear is a list,
    totalshear = []
    for n in [8, 9, 10, 11, 12]:
        shear = q_basicshear[n] + correctionalshearflows[1]
        totalshear.append(shear[0])
    for n in [13, 14, 0, 1]:
        shear = q_basicshear[n] + correctionalshearflows[0]
        totalshear.append(shear[0])
    for n in [3, 4, 5, 6, 7]:
        shear = q_basicshear[n] - correctionalshearflows[1]
        totalshear.append(shear[0])
    for n in [7]:
        shear = q_basicshear[n] + correctionalshearflows[0] + correctionalshearflows[1]
        totalshear.append(shear[0])

    return totalshear
