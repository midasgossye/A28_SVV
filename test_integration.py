# -*- coding: utf-8 -*-
"""
Created on Sun Feb 25 13:32:32 2018

@author: A28
"""

# TEST PROGRAM FOR ANALYTICAL INTEGRATION IN PYTHON


import sympy as sp
#from math import *
# Equation to solve:
#(-H_1_y/2)*x^2 + (q/2)*x^3 - (H_1_y/2)*x^2

# Define every unknown constant as an algebraic symbol:
H_1_y = sp.Symbol('H_1_y')
H_2_y = sp.Symbol('H_2_y')
H_3_y = sp.Symbol('H_3_y')
q = sp.Symbol('q')
x = sp.Symbol('x')
#INTEGRATE equation ANALYTICALLY w.r.t. x:
solution = sp.integrate(x*(sp.sin(x))**2, (x, 0, q))

# Print the analytical solution:
print "Solution:", solution

