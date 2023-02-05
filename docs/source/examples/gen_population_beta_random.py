
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
library_name = 'gen_beta_random' #Update it here
trial_index=0

# 4-O-5, alpha-O-4, beta-O-4, 5-5, beta-5, beta-beta, beta-1, Update it here
# lift the 5-5 percentage 
linkage_distribution_input = [7, 0, 75, 3, 7, 7, 0]


# H, G, S
monomer_distribution_input = [0, 30, 70]


verbose = True
additional_metrics = 0.0

branching_propensity = 0.5 # lift the branching propensity

# size in MW
size_in_MW = True

Tmetro = 5
Tmetro_out = 5


i_max = 2000
i_max_out = 5000
i_max_ring = 500

# total population size, Update it here
n_population = 10

seed_init = 1

trial_index = 0

# only 11 metrics in this case
metrics_weights = [1, 1, 1, 1, 10, 10, 1, 1, 1, 1, 1]
MW_upper = 25100
MW_lower = 300
population_metrics =  [MW_upper, MW_upper]


sim = opt.Simulation(linkage_distribution_input=linkage_distribution_input,
                     monomer_distribution_input=monomer_distribution_input,
                     use_beta_distribution=True,
                     max_size=MW_upper,
                     min_size=None,
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
                     population_metrics=None,
                     size_in_MW=size_in_MW,
                     metrics_weights=metrics_weights,
                     branching_propensity=branching_propensity,
                     form_rings=False,
                     show_plots=False)

sim.run()

P_population = sim.P_population
population = ch.Population(P_population, name=library_name, ResultsName=ResultsName, TrialIndex=str(trial_index))
population.analyze()

population_MWs = population.data['MW']

#%% plot the distributions 
from scipy.stats import beta
import matplotlib.pyplot as plt

fig, ax = plt.subplots(1, 1)

a, b = 2, 7.5
loc = 147
scale = 24907

r =beta.rvs(a, b, loc = loc, scale = scale, size=1000)
ax.hist(r, bins=20, alpha=0.5, color='blue', label='reference',density=True)
ax.hist(population_MWs, bins=20, alpha=0.5, color='red', label='actual',density=True)
ax.legend(loc='best', frameon=False)
plt.show()
