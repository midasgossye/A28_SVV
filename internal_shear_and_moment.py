import calculating_reaction_forces as reac

# ------------------------------------------
# -------INTERNAL-LOADS---------------------
# ------------------------------------------
# Everything is in the x,y-plane

# -------Defining constants-----------------
l_a = 2.691  # m
x_1 = 0.174  # m
x_2 = 1.051  # m
x_3 = 2.512  # m
x_a = 0.3  # m
q = 1.00 * 10 ** 3  # N/m
P = P = 20.6 * 10 ** 3  # N


# var = input('Please enter an x coordinate between 0 and 2.691 m:')
#
# if int(var) in range(0, int(x_1)):
#     m = (-q * var ** 2 / 2)
#     v = q * var
#     print 'The moment at point x=', var, 'equals', m, 'N/m'
#     print 'The shear at point x=', var, 'equals', v, 'N'
# elif int(var) in range(int(x_1), int(x_2)):
#     m = (H_1_y * var - q * var ** 2 / 2)
#     v = q * var - H_1_y
#     print 'The moment at point x=', var, 'equals', m, 'N/m'
#     print 'The shear at point x=', var, 'equals', v, 'N'
# elif int(var) in range(int(x_2), int(x_3)):
#     m = (H_2_y * var - q * var ** 2 / 2)
#     v = q * var - H_1_y - H_2_y
#     print 'The moment at point x=', var, 'equals', m, 'N/m'
#     print 'The shear at point x=', var, 'equals', v, 'N'
# else:
#     m = (H_3_y * var - q * var ** 2 / 2)
#     v = q * var - H_1_y - H_2_y - H_3_y
#     print 'The moment at point x=', var, 'equals', m, 'N/m'
#     print 'The shear at point x=', var, 'equals', v, 'N'


def internal(var, Izz):
    H_1_y, H_1_z, H_2_y, H_2_z, H_3_y, A_1_z = reac.calc_reac_f(Izz)

    # H_1_y = -55704  # N
    # H_2_y = 91291  # N
    # H_3_y = -32895  # N
    if int(var) > 2.691 or int(var) < 0:
        raise ValueError('internal shear and moment module error, x coordinates is out of bounds')
    if int(var) in range(0, int(x_1)):
        m = (-q * var ** 2 / 2)
        v = q * var
    elif int(var) in range(int(x_1), int(x_2)):
        m = (H_1_y * var - q * var ** 2 / 2)
        v = q * var - H_1_y
    elif int(var) in range(int(x_2), int(x_3)):
        m = (H_2_y * var - q * var ** 2 / 2)
        v = q * var - H_1_y - H_2_y
    else:
        m = (H_3_y * var - q * var ** 2 / 2)
        v = q * var - H_1_y - H_2_y - H_3_y
    return m, v
