
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
H = mono.Monomer('H')
H_G = H.create()

# Draw the connected graph
ut.draw_atomic_graph(H_G)
H_mol = ut.graph_to_mol(H_G, save_mol=True, name='H')
# initialize the polymer object with a monomer
polymer_h = poly.Polymer(H)
        
is_formed = polymer_h.add_specific_monomer('G', '5-5')
ut.draw_atomic_graph(polymer_h.G)

dimer_mol = ut.graph_to_mol(polymer_h.G)

# Add 5-5 ring with an extra monomer, the branding state needs to be True        
is_formed = polymer_h.add_specific_monomer('S', '5-5', True)
ut.draw_atomic_graph(polymer_h.G)
ut.draw_big_graph(polymer_h.bigG) 
trimer_mol = ut.graph_to_mol(polymer_h.G)

# add one more beta-O-4
is_formed = polymer_h.add_specific_monomer('H', 'beta-O-4', True)
ut.draw_atomic_graph(polymer_h.G)
ut.draw_big_graph(polymer_h.bigG) 
tetramer_mol = ut.graph_to_mol(polymer_h.G)