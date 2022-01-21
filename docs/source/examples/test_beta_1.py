"""
Test beta-1
"""


import ligning.characterization as ch
import ligning.monomer as mono
import ligning.polymer as poly
import ligning.utils as ut
from ligning.rules import default_color
import time 
from rdkit import Chem


#%% Step by step testing
linkage_new_name = 'beta-1'
linkage_new_types = ('H', 'H')

# try 8, 1 first, how about 1, 8???
linkage_new = (1, 8) # (8, 1)

H_G_1 = mono.Monomer(linkage_new_types[0]).create()
H_G_2 = mono.Monomer(linkage_new_types[1]).create()

dimmer_G = ut.join_two(H_G_1, H_G_2)

C1_index_in_polymer = linkage_new[0] - 1
C2_index_in_polymer = len(H_G_1) + linkage_new[1] - 1 
O_index_in_polymer = len(dimmer_G)

if linkage_new[0] == 8:
    alpha1_index_in_polymer = C1_index_in_polymer-1
    alpha2_index_in_polymer = C2_index_in_polymer-1
else: 
    alpha1_index_in_polymer = C2_index_in_polymer-1
    alpha2_index_in_polymer = C1_index_in_polymer-1

# Remove nodes 
atom_indices_to_delete = [7, 8, 9, 11]
atom_indices_to_delete_in_polymer = [ni + alpha2_index_in_polymer for ni in atom_indices_to_delete]
dimmer_G.remove_nodes_from(atom_indices_to_delete_in_polymer)
#ut.draw_graph(dimmer_G)


# Add the new -OH group
O_3_index = 16
# the properties are based on alpha1
monomer_index_1 = dimmer_G.nodes[alpha1_index_in_polymer]['mi']
monomer_type_1 = dimmer_G.nodes[alpha1_index_in_polymer]['mtype']
color = dimmer_G.nodes[alpha1_index_in_polymer]['color']

dimmer_G.add_node(O_index_in_polymer, element = 'O', aromatic = False, group = '7OH', index = O_3_index, mtype = monomer_type_1, bonding = False, color = color, mi=monomer_index_1)

# the old alpha beta bond
linkage_old_index = (alpha1_index_in_polymer, alpha1_index_in_polymer+1)

# set the alpha-beta bond order back to 1
dimmer_G.edges[linkage_old_index]['order'] = 1

# add the new linkage
linkage_new_index_list = [(C1_index_in_polymer, C2_index_in_polymer), (alpha1_index_in_polymer, O_index_in_polymer)]
dimmer_G.add_edges_from(linkage_new_index_list, order= 1, index = linkage_new, mtype = linkage_new_types, btype = linkage_new_name)
ut.draw_graph(dimmer_G)

dimmer_mol = ut.graph_to_mol(dimmer_G)

#%% Use the polymer object
H = mono.Monomer('H')
H_G_1 = H.create()

# initialize the polymer object with a monomer
polymer_h = poly.Polymer(H)

# add a bond monomer combination
monomer_i = 'H'
linkage_i = 'beta-1'
polymer_h.add_specific_monomer(monomer_i, linkage_i)

ut.draw_graph(polymer_h.G)
dimmer_mol = ut.graph_to_mol(polymer_h.G)

# note that the beta-1 is still not accurate
# there should have a -OOH group at the end
# in addition, we need to product a dimmer and attach lignin chain to the second part