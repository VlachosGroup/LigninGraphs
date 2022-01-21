"""
Test for visualizing atomic and big graphs
"""

import os
import numpy as np
import time

import ligning.monomer as mono
import ligning.polymer as poly
import ligning.utils as ut


import networkx as nx
import time 
from rdkit import Chem

output_path = os.path.join(os.getcwd(), 'results', 'molecules')

'''
Starting from a H monomer
'''
P0 = mono.Monomer("H")
polymer = poly.Polymer(P0)
G1 = polymer.G.copy()
ut.draw_atomic_graph(G1)

'''
Add a H monomer with beta-O-4 bond
'''
polymer.add_specific_linkage('beta-O-4', monomer_type = 'H')


#%%
'''
#Add n random monomers with beta-O-4 bonds and time the speed
'''
# 100 monomer takes ~0.9s
# 1000 monomer takes ~80s
start = time.time()
n_iter = 10
for i in range(n_iter):
    if i%2 ==0:
        polymer.add_specific_linkage('beta-O-4', monomer_type =  'G')
    else:
        polymer.add_specific_linkage('beta-O-4', monomer_type =  'S')
end = time.time()

n_monomer = n_iter + 2
run_time = end - start 

print("Lignin polymerization: %2d monomers takes %5.2f s" %(n_monomer, run_time))

# Plot the atomic and big graph
ut.draw_atomic_graph(polymer.G)
ut.draw_big_graph(polymer.bigG)

# Convert graph to smiles
Pn_smiles = ut.graph_to_smile(polymer.G)

# Convert smiles to mol and save the mol structure to an image file
Pn_mol = ut.graph_to_mol(polymer.G, save_mol=True, name='test_graphs_viz', save_path=output_path)
