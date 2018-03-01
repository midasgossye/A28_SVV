# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 09:28:16 2018

@author: Kaat
"""
from math import *
import numpy as np

h = 0.248
t_sk = 0.0011
C_a = 0.515
n_st = 11

circle_perim = 0.5 * pi * (0.5 * h - t_sk)
total_perimeter = circle_perim + sqrt((0.5 * h - t_sk) ** 2 + (C_a - 0.5 * h - t_sk) ** 2)  # m

spacing = total_perimeter / ((n_st + 1) / 2)

lengths = [spacing, spacing, spacing, spacing, total_perimeter - circle_perim - 4*spacing, 5*spacing - (total_perimeter - circle_perim), spacing, spacing, 5*spacing - (total_perimeter - circle_perim), total_perimeter - circle_perim - 4*spacing, spacing, spacing , spacing, spacing, h] 