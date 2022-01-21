"""
Test for ring formation inside a lignin polymer
"""
import ligning.characterization as ch
import ligning.monomer as mono
import ligning.polymer as poly
import ligning.utils as ut

import time 
from rdkit import Chem

#%% Constructing a polymer for testing
'''
Starting from a H monomer
'''
P0 = mono.Monomer("S")
P0_G = P0.create()
ut.draw_graph(P0_G)

polymer = poly.Polymer(P0)

'''
Add a G monomer with a linkage
'''
polymer.add_specific_linkage(linkage_type = '4-O-5', monomer_type = 'G')
P1_G = polymer.G
ut.draw_graph(P1_G)

'''
#Add n random monomers with beta-O-4 bonds and time the speed
'''
start = time.time()
n_iter = 10
for i in range(n_iter):
    polymer.add_specific_linkage(linkage_type ='beta-O-4', monomer_type =  'G')

n_monomer = n_iter + 2


#%% Add the rings
# Generate hypothetical structures with many rings
n_inter = 50 # try 50 iterations
for i in range(n_inter):
    polymer.add_random_ring()

Pn_G = polymer.G
Pn_bigG = polymer.bigG

ut.draw_graph(Pn_G)
ut.draw_graph(Pn_bigG)

# Convert graph to smiles
Pn_smiles = ut.graph_to_smile(Pn_G)
print(Pn_smiles)

# Convert smiles to mol and save the mol structure to an image file
Pn_mol = ut.graph_to_mol(Pn_G)
end = time.time()
run_time = end - start 
print("Lignin polymerization: %2d monomers takes %5.5f s" %(n_monomer, run_time))


#%% Characterize the polymer
cPn = ch.Characterize(polymer)

# Calculate and print the polymer properties
print("\nCharaterizing the polymer...\n")
cPn.cal_all()
