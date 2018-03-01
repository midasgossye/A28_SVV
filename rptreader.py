# imports
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from math import *

# functions
# main
DO228N_ULC1 = np.genfromtxt("Do228n_ULC1.rpt", skip_header=19, max_rows=6156)
DO228N_ULC1 = np.transpose(DO228N_ULC1)
DO228N_ULC1ft = np.genfromtxt("Do228n_ULC1.rpt", skip_header=6191, max_rows=16)
DO228N_ULC1ft = np.transpose(DO228N_ULC1ft)
DO228N_SLC1_1 = np.genfromtxt("Do228n_SLC1.rpt", skip_header=23, max_rows=6848)
DO228N_SLC1_1 = np.transpose(DO228N_SLC1_1)
DO228N_SLC1_2 = np.genfromtxt("Do228n_SLC1.rpt", skip_header=6892, max_rows=4280)
DO228N_SLC1_2 = np.transpose(DO228N_SLC1_2)
DO228N_SLC1_3 = np.genfromtxt("Do228n_SLC1.rpt", skip_header=11193, max_rows=6848)
DO228N_SLC1_3 = np.transpose(DO228N_SLC1_3)
DO228N_SLC1_4 = np.genfromtxt("Do228n_SLC1.rpt", skip_header=18062, max_rows=6848)
DO228N_SLC1_4 = np.transpose(DO228N_SLC1_4)

inp = np.genfromtxt("Do228n.inp", delimiter=",", comments="*", max_rows=6156)  # nodenumber, x, y, z
inp = np.transpose(inp)
# print max(inp[1])
# print min(inp[1])

# for i in xrange(len(inp) / 4):
#     sortedinp[0].append(inp[i * 4])
#     sortedinp[1].append(inp[i * 4 + 1])
#     sortedinp[2].append(inp[i * 4 + 2])
#     sortedinp[3].append(inp[i * 4 + 3])
# print sortedinp
# plotting room
# ax = Axes3D(plt.gcf())
# ax.scatter(sortedinp[3], sortedinp[1], sortedinp[2])
#
# plt.show()
