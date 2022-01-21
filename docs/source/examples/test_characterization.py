"""
Test for characterization the lignin structures
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
n_iter = 5
for i in range(n_iter):
    polymer.add_specific_linkage(linkage_type ='beta-O-4', monomer_type =  'G')
end = time.time()

n_monomer = n_iter + 2
run_time = end - start 

print("Lignin polymerization: %2d monomers takes %5.5f s" %(n_monomer, run_time))
Pn_G = polymer.G
ut.draw_graph(Pn_G)

# Convert graph to smiles
Pn_smiles = ut.graph_to_smile(Pn_G)
print(Pn_smiles)

# Convert smiles to mol and save the mol structure to an image file
Pn_mol = ut.graph_to_mol(Pn_G)

#%% Characterize the polymer

# Calculate and print the polymer properties using the graphs
cPn_G = ch.CharacterizeGraph(Pn_G)
print("\nCharaterizing the polymer...\n")
cPn_G.cal_all()


# Calculate and print the polymer properties using the graphs and big graphs
cPn = ch.Characterize(polymer)
print("\nCharaterizing the polymer...\n")
cPn.cal_all()