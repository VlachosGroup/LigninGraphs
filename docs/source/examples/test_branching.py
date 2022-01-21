"""
Test on branching propensity
"""

import ligning.characterization as ch
import ligning.monomer as mono
import ligning.polymer as poly
import ligning.utils as ut
from ligning.rules import linkage_ring
import time 
from rdkit import Chem

#%% Constructing a dimmer for testing
m0 = mono.Monomer("S")
polymer = poly.Polymer(m0)
polymer.add_specific_linkage(linkage_type = '4-O-5', monomer_type = 'G')
#%% Test on branched polymer
polymer_0 = poly.Polymer(polymer)
start = time.time()
n_iter = 10
for i in range(n_iter):
    flag = polymer_0.add_specific_linkage(linkage_type ='beta-O-4', monomer_type = 'G', branching_state=True)
    #polymer.add_specific_linkage(linkage_type ='5-5', monomer_type = 'G')

P_mol_0 = ut.graph_to_mol(polymer_0.G)
# Characterize the polymer
ch_p = ch.Characterize(polymer_0)
ch_p.cal_branching()
branching_coeff_0 = ch_p.branching_coeff
#assert(branching_coeff_0 == 0.25)

#%% Test on non-branched polymer
polymer_1 = poly.Polymer(polymer)
start = time.time()
n_iter = 10
for i in range(n_iter):
    polymer_1.add_specific_linkage(linkage_type ='beta-O-4', monomer_type = 'G', branching_state=False)
    #polymer.add_specific_linkage(linkage_type ='5-5', monomer_type = 'G')

P_mol_1 = ut.graph_to_mol(polymer_1.G)
# Characterize the polymer
ch_p = ch.Characterize(polymer_1)
ch_p.cal_branching()
branching_coeff_1 = ch_p.branching_coeff
#assert(branching_coeff_1 == 0)

#%% Test on drawing functions
ut.draw_graph(polymer_1.bigG)
ut.draw_graph(polymer_1.G)
ut.draw_atomic_graph(polymer_1.G)
ut.draw_big_graph(polymer_1.bigG)

