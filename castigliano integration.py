#######################
## CASTIGLIANO SHEAR ##
#######################

##------------------------------------ IMPORTS ---------------------------------------##

import sympy as sp
from math import *

##------------------------ PARAMETERS/VARIABLES/CONSTANTS-----------------------------##


#--UNKNOWN CONSTANTS AS SYMBOL--#

H_1_y = sp.Symbol('H_1_y')
H_2_y = sp.Symbol('H_2_y')
H_3_y = sp.Symbol('H_3_y')
q = sp.Symbol('q')
x = sp.Symbol('x')

#--LENGTHS--#

l_a = 2.691                                                       # span of the aileron
z_hingeline = (0.25*0.515-12.4)*cos(25*pi/180)                     # distance from body-fixed reference frame to hingeline of the aileron along z-dir. [m]
y_hingeline = 12.4*sin((90-25)*pi/180)-12.4*sin(25*pi/180)        # distance from body-fixed reference frame to hingeline of the aileron along y-dir. [m]
z_star = 12.4*cos(25*pi/180)                                      # relevant dimension for moment calculation [m]
x_1 = 0.174                                                       # location hinge 1 [m]
x_2 = 1.051                                                       # location hinge 2 [m]
x_3 = 2.512                                                       # location hinge 3 [m]
x_a = 0.3                                                         # distance between actuators I & II [m]
f_s = 0.4                                                         # assumed form factor for cross section
A  = 0.0021108225743942223                                        # cross sectional area [m2]
G  = 27*10**9                                                     # shear modulus aluminum 2014-T3 [Pa]

FR = (f_s)/(2*A*G)
integral_1 = sp.integrate(q*x,(x,0,x_1))
integral_2 = sp.integrate(q*x-H_1_y,(x,x_1,x_2-x_1))
integral_3 = sp.integrate(q*x-H_1_y-H_2_y,(x,x_2-x_1,x_3-x_2))
integral_4 = sp.integrate(q*x-H_1_y-H_2_y-H_3_y,(x,x_3-x_2,l_a))


#INTEGRATE equation ANALYTICALLY w.r.t. x:
C_i = FR*(integral_1+integral_2+integral_3+integral_4)

# Print the analytical solution:
print "C_i:", C_i
