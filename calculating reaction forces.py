import numpy as np
from math import *
from scipy.linalg import solve

#----------------------------------------------------------------
#----DEFINE CONSTANTS--------------------------------------------
#----------------------------------------------------------------

#----Given constants---------------------------------------------

P = 20.6*10**(3)    #N
q = 1.00*10**(3)    #N/m
C_a = 0.515         #m
l_a = 2.691         #m
x_1 = 0.174         #m
x_2 = 1.051         #m
x_3 = 2.512         #m
x_a = 0.3           #m
h = 0.248           #m
theta = 25*pi/180   #radians

#----Calculating distances---------------------------------------

y_hingeline = h/2*sin((pi/2)-theta) - h/2*sin(theta)        #m
z_hingeline = (0.25*C_a-h/2)*cos(theta)                     #m
z_star = h/2*cos(theta)                                     #m

#print y_hingeline
#print z_hingeline
#print z_star

#----------------------------------------------------------------
#----GOVERNING EQUATIONS-----------------------------------------
#----------------------------------------------------------------

#The equilibrium equations are implemented in 2 matrices, A and B
#Matrix A consists of H_1,y    H_1,z   H_2,y     H_2,z    H_3,y    A_1,z
#Matrix B consists of the known values on the right side of the governing equations
#Last equation will be castigliano's
A = np.matrix([[1, 0, 1, 0, 1, 0],                          #Matrix including the equilibruim equations with unknown reaction forces 
               [0, 1, 0, 1, 0, 1],
               [0, 0, 0, 0, 0, y_hingeline],
               [0, -x_1, 0, -x_2, 0, -(x_2-x_a/2)],
               [x_1, 0, x_2, 0, x_3, 0]])

B = np.matrix([[q*l],                                       #Matrix inlcuding the right side of the equilibrium equations, all the known forces
               [P],
               [q*l*z_hingeline+P*y_hingeline],
               [-P*(x_2+X_a/2)],
               [q*l**(2)/2]])



