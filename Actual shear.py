from math import *
import numpy as np
import sympy as sp

#aileron dimensions
C_a = 0.515
l_a = 2.691
h_a = 0.248
r_a = h_a/2
alpha = atan(r_a/(C_a-r_a))         #angle at TE of aileron
l_tr = (C_a-r_a)/(np.cos(alpha))       #length straight part triangle
#hinge and actuator locations
x1 = 0.174
x2 = 1.051
x3 = 2.512
x_a = 0.3
#skin properties
t_sk = 0.0011
t_sp = 0.0022
E = 71*10**9
G = 27*10**9
#stringers
t_st = 0.0012
h_st = 0.15
w_st = 0.03
n_st = 11
#load case
d1 = 0.1034/2.54
d3 = 0.2066/2.54
theta = 25*np.pi/180
P = -20600
q = -1000

#areas
A_st = t_st*w_st + (h_st-t_st)*t_sk     #stringer area
A_tr = t_sk*l_tr                        #area straight part skin 
A_sp = t_sp*h_a                         #area spar
A_cir = np.pi*r_a*t_sk                  #area semi circular skin
A_triangle = h_a*(C_a-r_a)/2            #enclosed area triangular part
A_cir_encl = 0.5*np.pi*r_a*r_a          #enclosed area semi circular part 

#stringer locations: y distance measured from symmetry line, z distance measured from TE & from spar 
d_st = (0.5*np.pi*r_a+l_tr)/6
y_st2 = r_a*np.sin(d_st/r_a)
y_st3 = (l_tr-(2*d_st-0.5*np.pi*r_a))*np.sin(alpha)
y_st4 = (l_tr-(3*d_st-0.5*np.pi*r_a))*np.sin(alpha)
y_st5 = (l_tr-(4*d_st-0.5*np.pi*r_a))*np.sin(alpha)
y_st6 = (l_tr-(5*d_st-0.5*np.pi*r_a))*np.sin(alpha)
z_st1 = C_a
z_st2 = (C_a-r_a) + r_a*np.cos(d_st/r_a)
z_st3 = (l_tr - 2*d_st + 0.5*np.pi*r_a)*np.cos(alpha)
z_st4 = (l_tr - 3*d_st + 0.5*np.pi*r_a)*np.cos(alpha)
z_st5 = (l_tr - 4*d_st + 0.5*np.pi*r_a)*np.cos(alpha)
z_st6 = (l_tr - 5*d_st + 0.5*np.pi*r_a)*np.cos(alpha)
z_str1 = z_st1 - (C_a - r_a)
z_str2 = z_st2 - (C_a - r_a)
z_str3 = -(C_a - r_a) + z_st3
z_str4 = -(C_a - r_a) + z_st4 
z_str5 = -(C_a - r_a) + z_st5 
z_str6 = -(C_a - r_a) + z_st6

#centroid z-location, measured from TE
centroidz = (2*A_tr*0.5*l_tr*np.cos(alpha) + A_cir*(C_a-r_a+2*r_a/np.pi) + A_sp*(C_a-r_a) + 2*A_st*z_st6 + 2*A_st*z_st5 + 2*A_st*z_st4 + 2*A_st*z_st3 + 2*A_st*z_st2 + A_st*z_st1)/(11*A_st + 2*A_tr + A_cir + A_sp)
centrzspar = centroidz - (C_a - r_a)

#moments of Inertia
I_zzbody = 1.46173*10**(-5)
I_yybody = 6.5726*10**(-5)

#forces at shear center
S_yb = 71308.4508
S_zb = 37311.32143
#T_x = -P*(abs(k)*sin(theta)-0.059977501)
k = 0.05984369622874337

#shear
x = sp.Symbol('x')
y = sp.Symbol('y')
q_bAB = -(S_zb*t_sk/I_yybody) * sp.integrate((abs(centrzspar)+r_a*sp.cos(y))*r_a, (y, 0, y)) - (S_yb*t_sk/I_zzbody) * sp.integrate(r_a*r_a*sp.sin(y), (y, 0, y))
q_bB1 = -(S_zb*t_sk/I_yybody) * sp.integrate((abs(centrzspar)+r_a*sp.cos(y))*r_a, (y, 0, 0.5*sp.pi)) - (S_yb*t_sk/I_zzbody) * sp.integrate(r_a*r_a*sp.sin(y), (y, 0, 0.5*sp.pi))
q_bCB = -(S_zb*t_sk/I_yybody) * sp.integrate(-(C_a-r_a-abs(centrzspar)-x*sp.cos(alpha)), (x,0,x)) - (S_yb*t_sk/I_zzbody) * sp.integrate(x*sp.sin(alpha), (x,0,x))
q_bB2 = -(S_zb*t_sk/I_yybody) * sp.integrate(-(C_a-r_a-abs(centrzspar)-x*sp.cos(alpha)), (x, 0, l_tr)) - (S_yb*t_sk/I_zzbody) * sp.integrate(x*sp.sin(alpha), (x, 0, l_tr))
q_bBD = q_bB1 + q_bB2 - (S_zb*t_sp/I_yybody) * sp.integrate(abs(centrzspar), (x, 0, x)) - (S_yb*t_sp/I_zzbody) * sp.integrate(r_a - x, (x, 0, x))
q_bD1 = q_bB1 + q_bB2 -(S_zb*t_sp/I_yybody) * sp.integrate(abs(centrzspar), (x, 0, h_a)) - (S_yb*t_sp/I_zzbody) * sp.integrate(r_a - x, (x, 0, h_a))
q_bCD = -(S_zb*t_sk/I_yybody) * sp.integrate(-(C_a-r_a-abs(centrzspar)-x*sp.cos(alpha)), (x, 0, x)) - (S_yb*t_sk/I_zzbody) * sp.integrate(-x*sp.sin(alpha), (x, 0, x))
q_bD2 = -(S_zb*t_sk/I_yybody) * sp.integrate(-(C_a-r_a-abs(centrzspar)-x*sp.cos(alpha)), (x, 0, l_tr)) - (S_yb*t_sk/I_zzbody) * sp.integrate(-x*sp.sin(alpha), (x, 0, l_tr))
q_bDA = q_bD1 + q_bD2 -(S_zb*t_sk/I_yybody) * sp.integrate(r_a*(abs(centrzspar)+r_a*sp.sin(y)), (y, 0, y)) - (S_yb*t_sk/I_zzbody) * sp.integrate(-r_a*r_a*sp.cos(y), (y, 0, y))
q_bA = q_bD1 + q_bD2 -(S_zb*t_sk/I_yybody) * sp.integrate(r_a*(abs(centrzspar)+r_a*sp.sin(y)), (y, 0, 0.5*sp.pi)) - (S_yb*t_sk/I_zzbody) * sp.integrate(-r_a*r_a*sp.cos(y), (y, 0, 0.5*sp.pi))
#q_bAD = -S_zb*t_sk/I_yybody * sp.integrate((abs(centrzspar)+r_a*sp.cos(y))*r_a, y) - S_yb*t_sk/I_zzbody * sp.integrate(-r_a*r_a*sp.sin(y), y)
#q_bD3 = -S_zb*t_sk/I_yybody * sp.integrate((abs(centrzspar)+r_a*sp.cos(y))*r_a, (y, 0, 0.5*sp.pi)) - S_yb*t_sk/I_zzbody * sp.integrate(-r_a*r_a*sp.sin(y), (y, 0, 0.5*sp.pi))

#solving for q_s01 and q_s02
coeff1q_s01 = A_triangle*(np.pi*r_a/t_sk + h_a/t_sp) + A_cir_encl*h_a/t_sp          #rate of twist
coeff1q_s02 = -A_triangle*h_a/t_sp - A_cir_encl*(h_a/t_sp + 2*l_tr/t_sk)            #rate of twist
coeff2q_s01 = 2*A_cir_encl                  #moment around middle point spar
coeff2q_s02 = 2*A_triangle                  #moment around middle point spar
constant1 = -(A_cir_encl/t_sk) * sp.integrate(q_bCB, (x, 0, l_tr)) + (A_cir_encl/t_sk) * sp.integrate(q_bCD, (x, 0, l_tr)) - (A_cir_encl/t_sp) * sp.integrate(q_bBD, (x, 0, h_a)) - (A_triangle*r_a/t_sk) * sp.integrate(q_bAB, (y, 0, 0.5*sp.pi)) - (A_triangle/t_sp) * sp.integrate(q_bBD, (x, 0, h_a)) - (A_triangle*r_a/t_sk) * sp.integrate(q_bDA, (y, 0, 0.5*sp.pi))
constant2 = -S_yb*k + sp.integrate(q_bCB*sp.sin(alpha)*(C_a-r_a), (x, 0, l_tr)) - sp.integrate(q_bCD*sp.sin(alpha)*(C_a-r_a), (x, 0, l_tr)) - sp.integrate(q_bAB*r_a*r_a, (y, 0, 0.5*sp.pi)) - sp.integrate(q_bDA*r_a*r_a, (y, 0, 0.5*sp.pi))
constant1 = float(constant1)
constant2 = float(constant2)
matA = np.matrix([[coeff1q_s01, coeff1q_s02], [coeff2q_s01, coeff2q_s02]])
matB = np.matrix([[constant1], [constant2]])
q_s0 = np.linalg.solve(matA, matB)
q_s01 = q_s0[0]
q_s02 = q_s0[1]

#total shear flow in skin
q_AB = q_bAB + q_s01
q_CB = q_bCB - q_s02
q_CD = q_bCD + q_s02
q_DA = q_bDA + q_s01
q_BD = q_bBD + q_s01 - q_s02

#shear flow in ribs
V_1y = sp.integrate((q_AB*sp.sin(0.25*sp.pi))*r_a, (y, 0, 0.5*sp.pi)) + sp.integrate((q_DA*sp.sin(0.25*sp.pi))*r_a, (y, 0, 0.5*sp.pi))
R_1y = float(V_1y - S_yb)
q_rib1 = R_1y/h_a

V_2y = sp.integrate(q_CB*sp.sin(alpha), (x, 0, l_tr)) - sp.integrate(q_CD*sp.sin(alpha), (x, 0, l_tr))
K_D2z = -(sp.integrate(q_CB*sp.cos(alpha)*r_a, (x, 0, l_tr)) - sp.integrate(q_CB*sp.sin(alpha)*(C_a-r_a), (x, 0, l_tr)) + sp.integrate(q_CD*sp.cos(alpha)*h_a, (x, 0, l_tr)) + sp.integrate(q_CD*sp.sin(alpha)*(C_a-r_a), (x, 0, l_tr)))/h_a
K_B2z = K_D2z
K_D2y = K_D2z*np.tan(alpha)
K_B2y = K_B2z*np.tan(alpha)
R_2y = V_2y - K_D2y - K_B2y
q_rib2 = R_2y/h_a
#print(q_s0)
#print(q_rib1)
#print(q_rib2)

#torsion
#coeff1 = A_cir_encl*h_a/t_sp + A_triangle*(pi*r_a/t_sk+h_a/t_sp)
#coeff2 = -(A_cir_encl*(2*l_tr/t_sk)+A_triangle*h_a/t_sp)
#matrixA = np.matrix([[coeff1, coeff2], [2*A_cir_encl, 2*A_triangle]])
#matrixB = np.matrix([[0], [T_x]])
#q_sT = np.linalg.solve(matrixA, matrixB)
#q_sT1 = q_sT[0]
#q_sT2 = q_sT[1]
#dthetadx = (q_sT1*(pi*r_a/t_sk + h_a/t_sp) - q_sT2*h_a/t_sp)/(2*A_cir_encl*G)
#dthetadx2 = (-q_sT1*h_a/t_sp + q_sT2*(2*l_tr/t_sk + h_a/t_sp))/(2*A_triangle*G)

#shear
#q_b12 = 0
#q_b23 = -S_zb/I_yybody * A_st * (z_st2 - centroidz) - S_yb/I_zzbody * A_st * y_st2
#q_b34 = q_b23 -S_zb/I_yybody * A_st * (z_st3 - centroidz) - S_yb/I_zzbody * A_st * y_st3
#q_b45 = q_b34 -S_zb/I_yybody * A_st * (z_st4 - centroidz) - S_yb/I_zzbody * A_st * y_st4
#q_b56 = q_b45 -S_zb/I_yybody * A_st * (z_st5 - centroidz) - S_yb/I_zzbody * A_st * y_st5
#q_b67 = q_b56 -S_zb/I_yybody * A_st * (z_st6 - centroidz) - S_yb/I_zzbody * A_st * y_st6
#q_b78 = q_b67 -S_zb/I_yybody * A_st * (z_st6 - centroidz) + S_yb/I_zzbody * A_st * y_st6
#q_b89 = q_b78 -S_zb/I_yybody * A_st * (z_st5 - centroidz) + S_yb/I_zzbody * A_st * y_st5
#q_b910 = q_b89 -S_zb/I_yybody * A_st * (z_st4 - centroidz) + S_yb/I_zzbody * A_st * y_st4
#q_b1011 = q_b910 -S_zb/I_yybody * A_st * (z_st3 - centroidz) + S_yb/I_zzbody * A_st * y_st3
#q_b111 = q_b1011 -S_zb/I_yybody * A_st * (z_st2 - centroidz) + S_yb/I_zzbody * A_st * y_st2