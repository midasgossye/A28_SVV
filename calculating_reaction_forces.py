import numpy as np
from math import *
from scipy.linalg import solve
import sympy as sp

# ----------------------------------------------------------------
# Program will return the 6 unknown reaction forces.


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
I = 1.5258310691614982e-05  # m**4  (Izz) # input MMOI

# ----Calculating distances---------------------------------------

y_hingeline = h / 2 * sin((pi / 2) - theta) - h / 2 * sin(theta)  # m
z_hingeline = (0.25 * C_a - h / 2) * cos(theta)  # m
z_star = h / 2 * cos(theta)  # m

# print y_hingeline
# print z_hingeline
# print z_star

# ----------------------------------------------------------------
# ----GOVERNING EQUATIONS-----------------------------------------
# ----------------------------------------------------------------

# Castigliano's theorem: Internal complimentary energy for bending only:
# ---------------------------------------------------------------

# Define every unknown constant as an algebraic symbol:
H_1_y = sp.Symbol('H_1_y')
H_2_y = sp.Symbol('H_2_y')
H_3_y = sp.Symbol('H_3_y')
# q = sp.Symbol('q')
x = sp.Symbol('x')


# I = sp.Symbol('I')
# E = sp.Symbol('E')

# #Squared bending moment equations for different sections of the beam
# m_1 = (-q*x**2/2)**2                            #for 0<x<x_1                                #for 0<x<x_1
# m_2 = (H_1_y*x-q*x**2/2)**2                     #for x_1<x<x_2
# m_3 = (H_2_y*x-q*x**2/2)**2                     #for x_2<x<x_3
# m_4 = (H_3_y*x-q*x**2/2)**2                     #for x_3<x<l_a
#
# #shear equations for different sections of the beam
# FR = (f_s)/(2*A*G)
# integral_1 = sp.integrate(q*x,(x,0,x_1))
# integral_2 = sp.integrate(q*x-H_1_y,(x,x_1,x_2-x_1))
# integral_3 = sp.integrate(q*x-H_1_y-H_2_y,(x,x_2-x_1,x_3-x_2))
# integral_4 = sp.integrate(q*x-H_1_y-H_2_y-H_3_y,(x,x_3-x_2,l_a))
#
# #Take integrals over dedicated distances
# integral1 = sp.integrate(m_1, (x, 0, x_1))
# integral2 = sp.integrate(m_2, (x, x_1, x_2))
# integral3 = sp.integrate(m_3, (x, x_2, x_3))
# integral4 = sp.integrate(m_4, (x, x_3, l_a))
#
# #Total complimentary energy due to bending
# CE_bending = (1/(2*E*I))*(integral1 + integral2 + integral3 + integral4)
#
# #total complimentary energy due to shear
# CE_shear = FR*(integral_1+integral_2+integral_3+integral_4)
#
# #Print the solution
# #print "Complementary energy due to bending:", CE_bending
#
# #Add complimentary energies of bending and shear
# CE_total = CE_bending + CE_shear
# #print "total complimentary energy", CE_total
#
# #Use castigliano's theorem
# #With respect to H_1_y
# delta_1 = sp.diff(CE_total, H_1_y)
# #print "deflection at hinge 1 equals", delta_1
# #With respect to H_3_y
# delta_3 = sp.diff(CE_total, H_3_y)
# #print "deflection at hinge 3 equals", delta_3
#
# #coefficienten op een rijtje
# a = sp.Poly(delta_1)
# #print a.coeffs()
#
# #The equilibrium equations are implemented in 2 matrices, A and B
# #Matrix A consists of H_1,y    H_1,z   H_2,y     H_2,z    H_3,y    A_1,z #output
# #Matrix B consists of the known values on the right side of the governing equations
# #sixth equation is castigliano's
# A = np.matrix([[1, 0, 1, 0, 1, 0],                          #Matrix including the equilibruim equations with unknown reaction forces
#                [0, 1, 0, 1, 0, 1],
#                [0, 0, 0, 0, 0, y_hingeline],
#                [0, -x_1, 0, -x_2, 0, -(x_2-x_a/2)],
#                [x_1, 0, x_2, 0, x_3, 0],
#                [1.9695*10**-6, 0, 2.87622*10**-6, 0, 4.38789*10**-7, 0]])
#
# B = np.matrix([[q*l_a],                                       #Matrix inlcuding the right side of the equilibrium equations, all the known forces
#                [P],
#                [q*l_a*z_hingeline+P*y_hingeline],
#                [-P*(x_2+x_a/2)],
#                [q*l_a**(2)/2],
#                 [0.13843]])
# # outputs 6 variables H_1,y    H_1,z   H_2,y     H_2,z    H_3,y    A_1,z
# forces = solve(A, B)
#
# # print forces
# # print 'These forces are not correct'

def calc_reac_f(Izz):
    I = Izz
    # Define every unknown constant as an algebraic symbol:
    H_1_y = sp.Symbol('H_1_y')
    H_2_y = sp.Symbol('H_2_y')
    H_3_y = sp.Symbol('H_3_y')
    # q = sp.Symbol('q')
    x = sp.Symbol('x')
    # I = sp.Symbol('I')
    # E = sp.Symbol('E')

    # Squared bending moment equations for different sections of the beam
    m_1 = (-q * x ** 2 / 2) ** 2  # for 0<x<x_1                                #for 0<x<x_1
    m_2 = (H_1_y * x - q * x ** 2 / 2) ** 2  # for x_1<x<x_2
    m_3 = (H_2_y * x - q * x ** 2 / 2) ** 2  # for x_2<x<x_3
    m_4 = (H_3_y * x - q * x ** 2 / 2) ** 2  # for x_3<x<l_a

    # shear equations for different sections of the beam
    FR = (f_s) / (2 * A * G)
    integral_1 = sp.integrate(q * x, (x, 0, x_1))
    integral_2 = sp.integrate(q * x - H_1_y, (x, x_1, x_2 - x_1))
    integral_3 = sp.integrate(q * x - H_1_y - H_2_y, (x, x_2 - x_1, x_3 - x_2))
    integral_4 = sp.integrate(q * x - H_1_y - H_2_y - H_3_y, (x, x_3 - x_2, l_a))

    # Take integrals over dedicated distances
    integral1 = sp.integrate(m_1, (x, 0, x_1))
    integral2 = sp.integrate(m_2, (x, x_1, x_2))
    integral3 = sp.integrate(m_3, (x, x_2, x_3))
    integral4 = sp.integrate(m_4, (x, x_3, l_a))

    # Total complimentary energy due to bending
    CE_bending = (1 / (2 * E * I)) * (integral1 + integral2 + integral3 + integral4)

    # total complimentary energy due to shear
    CE_shear = FR * (integral_1 + integral_2 + integral_3 + integral_4)

    # Print the solution
    # print "Complementary energy due to bending:", CE_bending

    # Add complimentary energies of bending and shear
    CE_total = CE_bending + CE_shear
    # print "total complimentary energy", CE_total

    # Use castigliano's theorem
    # With respect to H_1_y
    delta_1 = sp.diff(CE_total, H_1_y)
    # print "deflection at hinge 1 equals", delta_1
    # With respect to H_3_y
    delta_3 = sp.diff(CE_total, H_3_y)
    # print "deflection at hinge 3 equals", delta_3

    # coefficienten op een rijtje
    a = sp.Poly(delta_1)
    # print a.coeffs()

    # The equilibrium equations are implemented in 2 matrices, A and B
    # Matrix A consists of H_1,y    H_1,z   H_2,y     H_2,z    H_3,y    A_1,z #output
    # Matrix B consists of the known values on the right side of the governing equations
    # sixth equation is castigliano's
    matA = np.matrix([[1, 0, 1, 0, 1, 0],  # Matrix including the equilibruim equations with unknown reaction forces
                   [0, 1, 0, 1, 0, 1],
                   [0, 0, 0, 0, 0, y_hingeline],
                   [0, -x_1, 0, -x_2, 0, -(x_2 - x_a / 2)],
                   [x_1, 0, x_2, 0, x_3, 0],
                   [1.9695 * 10 ** -6, 0, 2.87622 * 10 ** -6, 0, 4.38789 * 10 ** -7, 0]])

    matB = np.matrix([[q * l_a],  # Matrix inlcuding the right side of the equilibrium equations, all the known forces
                   [P],
                   [q * l_a * z_hingeline + P * y_hingeline],
                   [-P * (x_2 + x_a / 2)],
                   [q * l_a ** (2) / 2],
                   [0.13843]])
    # outputs 6 variables H_1_y, H_1_z, H_2_y, H_2_z, H_3_y, A_1_z
    forces = solve(matA, matB)
    return forces
