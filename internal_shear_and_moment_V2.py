# ------------------------------------------
# -------INTERNAL-LOADS---------------------
# ------------------------------------------
# Everything is in the x,y-plane, except for z-dir. shear

import pylab as pl

# -------Defining constants-----------------

l_a = 2.691                         # m
x_1 = 0.174                         # m
x_2 = 1.051                         # m
x_3 = 2.512                         # m
x_a = 0.3                           # m
q = 1.00 * 10 ** 3                  # N/m
P = -20.6 * 10 ** 3                 # N

step = 0.00001                      # m

H_1_y = -55704.98965647             # N
H_1_z = -7079.78619831              # N
H_2_y = 91291.79830036              # N
H_2_z = 6886.63622554               # N
H_3_y = -32895.80864389             # N
A_1_z = 20793.14997276              # N


var = input('Please enter an x coordinate between 0 and 2.691 m:  ')

if var in list(pl.frange(0.0,x_1,step)):
    m = (-q * var ** 2 / 2)
    v_y = q * var
    v_z = 0
    print 'The moment at point x=', var, 'equals', m, 'N/m'
    print 'The shear_y at point x=', var, 'equals', v_y, 'N'
    print 'The shear_z at point x=', var, 'equals', v_z, 'N'
elif var in list(pl.frange(x_1,x_2-x_a/2,step)):
    m = (H_1_y * var - q * var ** 2 / 2)
    v_y = q * var - H_1_y
    v_z = H_1_z
    print 'The moment at point x=', var, 'equals', m, 'N/m'
    print 'The shear_y at point x=', var, 'equals', v_y, 'N'
    print 'The shear_z at point x=', var, 'equals', v_z, 'N'
elif var in list(pl.frange(x_2-x_a/2,x_2,step)):
    m = (H_1_y * var - q * var ** 2 / 2)
    v_y = q * var - H_1_y
    v_z = H_1_z + A_1_z
    print 'The moment at point x=', var, 'equals', m, 'N/m'
    print 'The shear_y at point x=', var, 'equals', v_y, 'N'
    print 'The shear_z at point x=', var, 'equals', v_z, 'N'
elif var in list(pl.frange(x_2,x_2+x_a/2,step)):
    m = (H_2_y * var - q * var ** 2 / 2)
    v_y = q * var - H_1_y - H_2_y
    v_z = H_1_z + A_1_z + H_2_z
    print 'The moment at point x=', var, 'equals', m, 'N/m'
    print 'The shear_y at point x=', var, 'equals', v_y, 'N'
    print 'The shear_z at point x=', var, 'equals', v_z, 'N'
elif var in list(pl.frange(x_2+x_a/2,x_3,step)):
    m = (H_2_y * var - q * var ** 2 / 2)
    v_y = q * var - H_1_y - H_2_y
    v_z = 0
    print 'The moment at point x=', var, 'equals', m, 'N/m'
    print 'The shear_y at point x=', var, 'equals', v_y, 'N'
    print 'The shear_z at point x=', var, 'equals', v_z, 'N'
else:
    m = (H_3_y * var - q * var ** 2 / 2)
    v_y = q * var - H_1_y - H_2_y - H_3_y
    v_z = 0
    print 'The moment at point x=', var, 'equals', m, 'N/m'
    print 'The shear_y at point x=', var, 'equals', v_y, 'N'
    print 'The shear_z at point x=', var, 'equals', v_z, 'N'



    

def interval(var):

    if int(var)>2.691 or int(var)<0:
        raise ValueError('internal shear and moment module error, x coordinates is out of bounds')

    
    if var in list(pl.frange(0.0,x_1,step)):
        m = (-q * var ** 2 / 2)
        v_y = q * var
        v_z = 0
    elif var in list(pl.frange(x_1,x_2-x_a/2,step)):
        m = (H_1_y * var - q * var ** 2 / 2)
        v_y = q * var - H_1_y
        v_z = H_1_z
    elif var in list(pl.frange(x_2-x_a/2,x_2,step)):
        m = (H_1_y * var - q * var ** 2 / 2)
        v_y = q * var - H_1_y
        v_z = H_1_z + A_1_z
    elif var in list(pl.frange(x_2,x_2+x_a/2,step)):
        m = (H_2_y * var - q * var ** 2 / 2)
        v_y = q * var - H_1_y - H_2_y
        v_z = H_1_z + A_1_z + H_2_z
    elif var in list(pl.frange(x_2+x_a/2,x_3,step)):
        m = (H_2_y * var - q * var ** 2 / 2)
        v_y = q * var - H_1_y - H_2_y
        v_z = 0
    else:
        m = (H_3_y * var - q * var ** 2 / 2)
        v_y = q * var - H_1_y - H_2_y - H_3_y
        v_z = 0 

        
    return m, v_y, v_z
