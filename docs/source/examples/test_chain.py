# -*- coding: utf-8 -*-
"""
test chain formation
"""
import ligning.characterization as ch
import ligning.monomer as mono
import ligning.polymer as poly
import ligning.utils as ut
from ligning.rules import linkage_ring
import time 
from rdkit import Chem

#%% Constructing a polymer for testing
'''
Starting from a H monomer
'''
P0 = mono.Monomer("G")
P0_G = P0.create()
ut.draw_graph(P0_G)

polymer = poly.Polymer(P0)

'''
Add a G monomer with a linkage
'''
polymer.add_specific_linkage(linkage_type = 'beta-O-4', monomer_type = 'G')
P1_G = polymer.G
ut.draw_graph(P1_G)

'''
#Add n random monomers with beta-O-4 bonds and time the speed
'''
start = time.time()
n_iter = 100
for i in range(n_iter):
    #polymer.add_specific_linkage(linkage_type ='beta-O-4', monomer_type = 'G')
    polymer.add_specific_linkage(linkage_type ='beta-O-4', monomer_type = 'G')

n_iter = 10
for i in range(n_iter):
    #polymer.add_specific_linkage(linkage_type ='beta-O-4', monomer_type = 'G')
    polymer.add_specific_linkage(linkage_type ='beta-beta', monomer_type = 'G')
    
    
n_monomer = n_iter + 2
P_mol = ut.graph_to_mol(polymer.G)

#%% Characterize the polymer

# Calculate and print the polymer properties using the graphs
ch1 = ch.CharacterizeGraph(polymer.G)
print("\nCharaterizing the polymer...\n")
ch1.cal_all()

# Calculate and print the polymer properties using the graphs and big graphs
ch2 = ch.Characterize(polymer)
print("\nCharaterizing the polymer...\n")
ch2.cal_all()