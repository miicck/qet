#
# This module provides routines for the optimization
# of lennard-jones systems
#
import numpy as np
from scipy.optimize import minimize

LJ_UNIT_CELL_CUTOFF = 3

# Calculate the lennard-jones potential for a
# lattice, fractional coordinates (x0), minima
# distances (rms) and depths (epsilons)
def lj_potential(lattice, x0, rms, epsilons):

    pot = 0.0
    
    # Index of first atom
    for i in range(0, len(x0)):

        # Position of first atom
        x1  = x0[i][0] * lattice[0]
        x1 += x0[i][1] * lattice[1]
        x1 += x0[i][2] * lattice[2]

        # Index of second atom
        for j in range(0, i):

            # Loop over periodic images
            for nx in range(-LJ_UNIT_CELL_CUTOFF, LJ_UNIT_CELL_CUTOFF+1):
                for ny in range(-LJ_UNIT_CELL_CUTOFF, LJ_UNIT_CELL_CUTOFF+1):
                    for nz in range(-LJ_UNIT_CELL_CUTOFF, LJ_UNIT_CELL_CUTOFF+1):

                        # Position of second atom
                        x2  = (x0[j][0] + nx) * lattice[0]
                        x2 += (x0[j][1] + ny) * lattice[1]
                        x2 += (x0[j][2] + nz) * lattice[2] 

                        # Use the average values of epsilon/rm
                        # for these two atoms
                        rm = 0.5*(rms[i] + rms[j])
                        ep = 0.5*(epsilons[i] + epsilons[j])
                        r  = np.linalg.norm(x1-x2)
                        pot += rm/r
                        continue

                        # Calculate the LJ pairwise potential
                        r6 = (rm/r)**6.0
                        pot += ep*(r6**2.0 - 2*r6)
    return pot

# Unflatten the vector of fractional coords
def unflatten_x(x):
    for i in range(0, len(x)): x[i] -= np.floor(x[i])
    n1 = int(len(x)/3)
    return x.reshape((n1,3))

# Optimize a lennard-jones crystal with the given
# lattice, fractional coordinates (x0), minima
# distances (rms) and depths (epsilons)
# if epsilons = None, we set espilon_i = 1 for all i
def optimize_lj_crystal(lattice, x0, rms, epsilons=None):

    if epsilons is None:
        epsilons = [1.0 for x in x0]

    lattice = np.array(lattice)
    x0      = np.array(x0)

    to_min  = lambda x : lj_potential(lattice, unflatten_x(x), rms, epsilons)
    min_res = minimize(to_min, x0.flatten())
    print(min_res)
    xf = unflatten_x(min_res.x)
    return [lattice, xf]

def plot_crystal(ax, lattice, x0):
    ax.axvline(0.0)
    ax.axvline(1.0)
    ax.axhline(0.0)
    ax.axhline(1.0)
    for x in x0:
        ax.plot([x[0]], [x[1]], marker="+", markersize=20-10*x[2], color="red")
    ax.set_xlim([-0.1,1.1])
    ax.set_ylim([-0.1,1.1])

# Test the above
def test():

    import matplotlib.pyplot as plt
    
    # Put some atoms in a cubic box
    lattice = np.identity(3)*2
    x0      = [[0.4,0.4,0.4], [0.5,0.5,0.5], [0.6,0.6,0.6]]
    rms     = [1.0 for x in x0]
    eps     = [1.0 for x in x0]

    plot_crystal(plt.subplot(221), lattice, x0)

    lattice, x = optimize_lj_crystal(lattice, x0, rms, eps)

    plot_crystal(plt.subplot(222), lattice, x)

    plt.show()

test()
