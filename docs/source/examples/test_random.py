"""
Test random structure generation
"""

import numpy as np
import time

import ligning.monomer as mono
import ligning.polymer as poly
import ligning.utils as ut

# Set the random seed
np.random.seed(1)

# Starting from a random monomer
P0 = mono.Monomer('S')
polymer = poly.Polymer(P0)

start = time.time()

# Set the number of polymerization iteration
n_iter = 12

# list for graphs
P_Gs = []

# list for mols
P_mols = []

# Polymerization loop
for i in range(n_iter):
    polymer.add_random_monomer()
    Px = polymer.G.copy()
    P_Gs.append(Px)
    Px_mol = ut.graph_to_mol(Px)
    P_mols.append(Px_mol)

end = time.time()

n_monomer = n_iter + 2
run_time = end - start 

print("Lignin polymerization: %2d monomers takes %5.5f s" %(n_monomer, run_time))


