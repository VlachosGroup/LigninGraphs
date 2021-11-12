"""
Miscanthus simulation
"""

import ligning.optimization as opt
#%% Set the optimization target
library_name = 'miscanthus_weights'

# linkage and monomer distribution
# ['4-O-5', 'alpha-O-4', 'beta-O-4', '5-5', 'beta-5', 'beta-beta', 'beta-1']
linkage_distribution_input =  [0, 0, 68, 0, 15, 17, 0] 

# ['H', 'G', 'S']
monomer_distribution_input = [4, 46, 50]

# additional metrics: branching coeffcient
additional_metrics = 0.0

# branching propensity
branching_propensity = 0.0

# population metrics: number average MW, weight average MW
population_metrics = [1240, 2310]

# size limit, monomer count or number average MW
expected_size = 1240
max_size = 10000

# size distribution scaling factor
distribution_scaling = 0.1

# size in MW
size_in_MW = True

# Metropolis temperature
Tmetro = 5
Tmetro_out = 5

# number of iterations
i_max = 1000
i_max_out = 1000
i_max_ring = 500

# max polymer in population
n_population = 100

# random seed init
seed_init = 1

# weights for 13 metrics in total 
# the last two correspond to number average MW, weight average MW
# we assign them more weights
metrics_weights = [1, 1, 1, 1, 10, 10, 1, 1, 1, 1, 1, 10, 10]

#%% Simulation main body
sim = opt.Simulation(linkage_distribution_input=linkage_distribution_input,
                     monomer_distribution_input=monomer_distribution_input,
                     expected_size=expected_size,
                     max_size=max_size,
                     distribution_scaling=distribution_scaling,
                     Tmetro=Tmetro,
                     Tmetro_out=Tmetro_out,
                     seed_init=seed_init,
                     library_name=library_name,
                     n_population=n_population,
                     i_max=i_max,
                     i_max_out=i_max_out,
                     i_max_ring=i_max_ring,
                     additional_metrics=additional_metrics,
                     population_metrics=population_metrics,
                     size_in_MW=size_in_MW, 
                     metrics_weights=metrics_weights,
                     branching_propensity=branching_propensity,
                     show_plots=False)

sim.run()
