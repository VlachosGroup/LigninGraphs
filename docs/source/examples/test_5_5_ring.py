
# set ligning path (optional if installed via pip)
import sys, os
project_path = os.path.abspath(os.path.join(os.getcwd(), '../../..'))
print(project_path)
sys.path.insert(0, project_path)

import os
import ligning.monomer as mono
import ligning.polymer as poly
import ligning.utils as ut
from ligning.rules import monomer1_select_monomer2_linkage_name
import time 
from rdkit.Chem.Descriptors import ExactMolWt

output_path = os.path.join(os.getcwd(), 'demo_results', 'molecules')

H_G = mono.monomer_H()
G_G = mono.monomer_G()
S_G = mono.monomer_S()

# Initialize monomer objects
G = mono.Monomer('G')
G_G = G.create()

# Draw the connected graph
ut.draw_atomic_graph(G_G)
G_mol = ut.graph_to_mol(S_G, save_mol=True, name='G')
# initialize the polymer object with a monomer
polymer_h = poly.Polymer(G)

is_formed = polymer_h.add_specific_monomer('G', 'beta-O-4')
ut.draw_atomic_graph(polymer_h.G)
print(polymer_h.find_available_C1_in_polymer())
        
is_formed = polymer_h.add_specific_monomer('G', '5-5')
ut.draw_atomic_graph(polymer_h.G)
print(polymer_h.find_available_C1_in_polymer())

dimer_mol = ut.graph_to_mol(polymer_h.G)

edges_5_5 = [ei for ei in polymer_h.G.edges if polymer_h.G.edges[ei]['btype'] == '5-5']

# Add 5-5 ring with an extra monomer, the branding state needs to be True   
print(polymer_h.O_pair_indices_in_polymer)     
is_formed = polymer_h.add_specific_monomer('S', '5-5-ring', True)
ut.draw_atomic_graph(polymer_h.G)
ut.draw_big_graph(polymer_h.bigG) 
trimer_mol = ut.graph_to_mol(polymer_h.G)
trimer_G = polymer_h.G.copy()
C1_available = polymer_h.find_available_C1_in_polymer()
print(C1_available)
# [0, 3, 6, 7, 11, 14, 17, 18]

# add one more beta-O-4
is_formed = polymer_h.add_specific_monomer('H', '4-O-5', True)
ut.draw_atomic_graph(polymer_h.G)
ut.draw_big_graph(polymer_h.bigG) 
tetramer_mol = ut.graph_to_mol(polymer_h.G)
tetramer_G = polymer_h.G.copy()
C1_available = polymer_h.find_available_C1_in_polymer()
print(C1_available)

# add one more beta-O-4
is_formed = polymer_h.add_specific_monomer('H', 'beta-O-4', True)
ut.draw_atomic_graph(polymer_h.G)
ut.draw_big_graph(polymer_h.bigG) 
x_mol = ut.graph_to_mol(polymer_h.G)

