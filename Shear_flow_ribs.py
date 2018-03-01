import numpy as np
from math import *
from scipy.linalg import solve
import sympy as sp


#----------------------------------------------------------------
#Shear calculations in ribs


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
f_s = 0.4           # assumed form factor for cross section
A  = 0.0021108225743942223                                        # cross sectional area [m2]
G = 28*10**9        #shear modulus aluminum 2014-T3 [Pa]
d_1 = 0.134         #m
d_3 = 0.2066        #m
E = 73.1*10**9      #Pa

#De volgende constanten staan denk ik ook ergens in het grote main bestand?
phi = 0.3071        #radians, the angle between the symmetry axis and the top or bottom of the triangle of the cross-section

#----------------------------------------------------------------
#----CALCULATING SHEAR FLOW IN RIBS------------------------------
#----------------------------------------------------------------

#Consider every cell seperately, first the left cell, only in vertical direction
#beta 1 is half the angle between spar and boom 2, beta 2 is half the angle between boom 2 and boom 0
beta1 = (lengths[5]/(h/2))/2
beta2 = (lengths[6]/(h/2))/2

#calculate the vertical force due to the total shearflows in the left cell
V_yprime_1 = totalshear[5]*lengths[5]*sin(beta1/2) + totalshear[6]*lengths[6]*sin(beta2/2) + totalshear[7]*lengths[7]*sin(beta2/2) + totalshear[8]*lengths[9]*sin(beta1/2)

#Calculate the vertical due to the total shearflows in the right cell
V_yprime2 = 0
for n in [0, 4]:
    Verticalforce = totalshear[n] * lengths[n] * sin(phi)
    V_yprime2 = V_yprime2 + Verticalforce
for n in [9, 13]:
    Verticalforce = totalshear[n] * lengths[n] * sin(phi)
    V_yprime2 = V_yprime2 + Verticalforce

#Calculate the force in the ribs in the left cell
F_ribs_left = V_prime_1 - S_yprime      #N

#Calculate the forces carried by the flanges in cell 2, by taking moment around point 1, clockwise positive
#Start with calculating the force per shear section
forcesperboom = []
for n in [0, 4]:
    #vertical direction
    verticalforce = totalshear[n]*lengths[n]*sin(phi)
    #horizontal direction
    horizontalforce = totalshear[n]*lengths[n]*cos(phi)
    #append to list with [horizontalforce, verticalforce] per boom
    forcesperboom.append(horizontalforce, verticalforce)


#Compute moment of first shear flow by hand, because it is not in the list with coordinates
momentspersection = [forcesperboom[0][1]*(l_a-(h/2)+forcesperboom[0][0]* (h/2)]
#Compute the rest of the moments and append to list momentspersection
for n in [1, 4]:
    momentofshearflow = verticalforce * -z_y_angle_coords[n][0] + horizontalforce* (-z_y_angle_coords[n][1]+h/2)
    momentspersection.append(momentofshearflow)
#Compute total moment due to shear flows around point 1
totalmoment_1 = sum(momentspersection)
#Do moment equilibrium and rewrite to K_d, the horizontal component of K_d
K_d = totalmoment_1/h
#Compute the vertical component of K_d
K_d_vertical = K_d*tan(phi)

#Calculate the vertical component of K_b, by taking equilibrium in the z' direction of the external forces only
K_b_vertical = -K_d_vertical

#Compute force carried by rib
Rib_force = V_yprime2 - K_d_vertical - K_b_vertical

#Finally, compute the shear flow in the ribs
q_ribs = Rib_force/h    #N/m

