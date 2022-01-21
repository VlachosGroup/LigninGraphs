"""
Test for ring formation inside a lignin polymer
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
    #polymer.add_specific_linkage(linkage_type ='beta-O-4', monomer_type = 'G')
    polymer.add_specific_linkage(linkage_type ='5-5', monomer_type = 'G')

n_monomer = n_iter + 2

 
for linkage_i in linkage_ring:
    print(linkage_i)
    polymer.add_specific_ring(linkage_i)