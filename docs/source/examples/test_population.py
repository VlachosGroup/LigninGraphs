
# set ligning path (optional if installed via pip)
import sys, os
project_path = os.path.abspath(os.path.join(os.getcwd(), '../../..'))
print(project_path)
sys.path.insert(0, project_path)

import os
import ligning.monomer as mono
import ligning.polymer as poly
import ligning.utils as ut
import ligning.optimization as opt
import ligning.characterization as ch

from ligning.rules import monomer1_select_monomer2_linkage_name
import time 
from rdkit.Chem.Descriptors import ExactMolWt

output_path = os.path.join(os.getcwd(), 'demo_results', 'molecules')

ResultsName='demo_results'
library_name = 'new_branch'
trial_index=0

# 4-O-5, alpha-O-4, beta-O-4, 5-5, beta-5, beta-beta, beta-1
# lift the 5-5 percentage 
linkage_distribution_input =[7, 0, 75, 20, 7, 7, 0] #[7, 0, 75, 3, 7, 7, 0]


# H, G, S
monomer_distribution_input = [0, 30, 70]


verbose = True
additional_metrics = 0.0

branching_propensity = 0.5 # lift the branching propensity

population_metrics = [2540, 4170]


expected_size = 2540
max_size = 12000

# size distribution scaling factor
distribution_scaling = 0.15

# size in MW
size_in_MW = True

Tmetro = 5
Tmetro_out = 5


i_max = 1000
i_max_out = 1000
i_max_ring = 500

n_population = 10

seed_init = 1

trial_index = 0

metrics_weights = [1, 1, 1, 1, 10, 10, 1, 1, 1, 1, 1, 10, 10]

sim = opt.Simulation(linkage_distribution_input=linkage_distribution_input,
                     monomer_distribution_input=monomer_distribution_input,
                     use_beta_distribution=False,
                     use_uniform_distribution=False,
                     expected_size=expected_size,
                     max_size=max_size,
                     distribution_scaling=distribution_scaling,
                     Tmetro=Tmetro,
                     Tmetro_out=Tmetro_out,
                     seed_init=seed_init,
                     ResultsName=ResultsName,
                     library_name=library_name,
                     trial_index=trial_index,
                     n_population=n_population,
                     i_max=i_max,
                     i_max_out=i_max_out,
                     i_max_ring=i_max_ring,
                     additional_metrics=additional_metrics,
                     population_metrics=population_metrics,
                     size_in_MW=size_in_MW,
                     metrics_weights=metrics_weights,
                     branching_propensity=branching_propensity,
                     form_rings=False,
                     show_plots=False)

sim.run()

P_population = sim.P_population
population = ch.Population(P_population, name=library_name, ResultsName=ResultsName, TrialIndex=str(trial_index))
population.analyze()
#%%
population_MWs = population.data['MW']

#%%
# see ut.graph_to_mol(P_population[5].G) for the 5-5
x_mol = ut.graph_to_mol(P_population[5].G)
x_ch = ch.Characterize(P_population[5]).count_linkages()
ut.draw_big_graph(P_population[5].bigG)
