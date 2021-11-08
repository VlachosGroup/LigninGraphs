"""
Popular simulation
"""

import ligning.optimization as opt
#%% Set the optimization target
library_name = 'poplar_weights'

# linkage and monomer distribution
# ['4-O-5', 'alpha-O-4', 'beta-O-4', '5-5', 'beta-5', 'beta-beta', 'beta-1']
linkage_distribution_input = [0, 0, 78, 0.1, 7, 15, 0] 

# ['H', 'G', 'S']
monomer_distribution_input = [0, 37, 63]

# additional output during simulation
verbose = True

# additional metrics: branching coeffcient
additional_metrics = True

# population metrics: number average MW, weight average MW
population_metrics = None

# size limit, monomer count or number average MW
expected_size = 50
max_size = 250

# size distribution scaling factor
distribution_scaling = .1

# size in MW
size_in_MW = False

# branching propensity
branching_propensity = 0.1

# Metropolis temperature
Tmetro = 10
Tmetro_out = 10

# number of iterations
i_max = 5000
i_max_out = 1000
i_max_ring = 500

# max polymer in population
n_population = 100

# random seed init
seed_init = 1

# weights for 11 metrics in total 
metrics_weights = [1, 1, 1, 1, 1, 100, 1, 1, 1, 1, 1]

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
                     branching_propensity=branching_propensity,
                     metrics_weights = metrics_weights,
                     verbose=verbose)

sim.run()
