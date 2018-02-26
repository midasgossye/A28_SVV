##################################################
######### INTERNAL SHEAR CALCULATION #############
##################################################


##-------------------------------------------------IMPORTS--------------------------------------------##


from math import *



##------------------------------------------------BEGIN TESTING---------------------------------------##



value = input('ENTER X-VALUE ON BEAM: ')


##------------------------------------------------CONSTANTS-------------------------------------------##


l_a = 2.691                                                     #[m]
x_1 = 0.174                                                     #[m]
x_2 = 1.051                                                     #[m]
x_3 = 2.512                                                     #[m]

q   = 1                                                         #[kN/m]
H_1_y = -55704                                                  #[N]
H_2_y = 91291                                                   #[N]
H_3_y = -32895                                                  #[N]


##---------------------------------------------SHEAR COMPUTATION--------------------------------------##

for value in range(0,int(x_1)):

    v = q*value

for value in range(int(x_1),int(x_2)):

    v = q*value - H_1_y

for value in range(int(x_2),int(x_3)):

    v = q*value - H_1_y - H_2_y

for value in range(int(x_3),int(l_a)):

    v = q*value - H_1_y - H_2_y - H_3_y


print 'shear at x-value',value,'is equal to',v
