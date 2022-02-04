from networkx.utils.decorators import random_state
import numpy as np
import time
from numpy.core.fromnumeric import size
import scipy.stats as stats
import os
import shutil

try:
    from nextorch import plotting, bo, doe
except ImportError:
    pass

import copy
from typing import Optional, Tuple, TypeVar
import pprint
pp = pprint.PrettyPrinter(indent=4)

from ligning.rules import monomer_types, linkage_names
import ligning.utils as ut
import ligning.characterization as ch
from ligning.monomer import Monomer
from ligning.polymer import Polymer
from ligning.characterization import Characterize
from ligning.utils import nxgraph, molecule, nparray


# Boltzmann constant
kb = 8.617333262145e-05

# Distribution factor from Boltzmann distributin
w = 0

# Temperature of the simulation in K
Tmax = 2000 #K the max temperature, starting temperature


#%% MCMC simulation classes 
class Trajectory():
    """ Single MCMC trajectory """ 
    def __init__(
        self, 
        linkage_distribution_input: list, 
        monomer_distribution_input: list, 
        Tmetro: int,
        expected_size: Tuple[int, float], 
        max_size: Tuple[int, float],
        distribution_scaling: Optional[float] = 1.0,
        additional_metrics: Optional[Tuple[list, float]] = None,
        size_in_MW: Optional[bool] = False,
        branching_propensity: Optional[float] = None,
        metrics_weights: Optional[Tuple[list, nparray]] = None,
        verbose: Optional[bool] = True,
        file: Optional[str] = None):
        """Initialize the target metrics

        Parameters
        ----------
        linkage_distribution_input : list
            target linkage distribution
        monomer_distribution_input : list
            target monomer distribution
        Tmetro : int
            temperature for metropolis monte carlo
        expected_size : Tuple[int, float]
            target polymer size (mean)
        max_size : Tuple[int, float]
            maximum polymer size
        distribution_scaling : Optional[float], optional
            std of the size distribution, by default 1.0
        additional_metrics : Optional[Tuple[list, float]], optional
            additional metrics such as branching coefficient, by default None
        size_in_MW : Optional[bool], optional
            flag indicating whether size is in MW, by default False
        branching_propensity : Optional[float], optional
            branching propensity, by default None
        metrics_weights : Optional[Tuple[list, nparray]], optional
            assigned weights for each metrics, by default None
        verbose : Optional[bool], optional
            flag to control print outputs, by default True
        file: Optional[str], optional
            output file path, by default None
        """        
   
        self.linkage_distribution = np.array(linkage_distribution_input)/np.sum(linkage_distribution_input)
        self.monomer_distribution = np.array(monomer_distribution_input)/np.sum(monomer_distribution_input)
        self.expected_size = expected_size
        self.max_size = max_size
        self.distribution_scaling = distribution_scaling
        self.additional_metrics = None
        self.additional = False # flag for addditional metrics
        self.metrics_weights = metrics_weights # an array of metrics weights

        metrics_target = np.concatenate((self.monomer_distribution, self.linkage_distribution), axis= 0)

        # Add additional metrics if it is provided
        if additional_metrics is not None:
            if isinstance(additional_metrics, float) or isinstance(additional_metrics, int) :
                additional_metrics = [additional_metrics]

            additional_metrics =  np.array(additional_metrics)
            
            metrics_target = np.concatenate((metrics_target, additional_metrics), axis= 0)
            self.additional = True

        self.metrics_target = metrics_target
        self.additional_metrics = additional_metrics
        
        # metropolis temperature
        self.Tmetro = Tmetro

        # flag for whether moleculcar weight is include in the metrics
        self.size_in_MW = size_in_MW
        self.cal_MW = True

        # the likelihood to form a branch
        self.branching_propensity = branching_propensity

        # print flag
        self.verbose = verbose

        # hard monomer count limit
        self.max_monomer_count = 100

        # if output file path is provided
        self.file = file



    def run_MCMC(
        self, 
        rseed: int, 
        i_max: Optional[int] = 500
    ) -> Tuple[Polymer, dict]:
        """run MCMC simulations

        Parameters
        ----------
        rseed : int
            random seed
        i_max : Optional[bool], optional
            maximum number of iterations in the inner loop, 
            by default 500

        Returns
        -------
        polymer_final : Polymer
            final polymer object
        MC_dict : dict
            MC trajectory
        
        """        
        '''
        Starting from a random monomer
        '''
        random_state = ut.set_random_state(rseed)
        # if self.verbose:
        #     print("\tCreating a new polymer:")
        # if self.file is not None:
        #     self.file.write("\tCreating a new polymer:\n")

        # Select one random size as the stopping criteria 
        self.stop_size = ut.generate_random_size_from_distribution(self.expected_size, 
                                                                    self.max_size, 
                                                                    distribution_scaling=self.distribution_scaling,
                                                                    size_in_MW=self.size_in_MW)

        # set the branching state
        branching_state = None # no input, no restrictions

        monomer_random = ut.generate_random_monomer(self.monomer_distribution, random_state)
        M_init = Monomer(monomer_random)
        polymer = Polymer(M_init, False) #, verbose= self.verbose)
        
        metrics_P, monomer_count_P, MW_P = ch.get_metrics_polymer(polymer, 
                                                                       additional=self.additional, 
                                                                       cal_MW=self.cal_MW)
        distance_init = ut.cal_distance(self.metrics_target, metrics_P, self.metrics_weights)

        '''
        Form a dimmer
        '''
        linkage_new = ut.generate_random_linkage(self.linkage_distribution, random_state)
        monomer_new = ut.generate_random_monomer(self.monomer_distribution, random_state)
        new_linkage_flag = polymer.add_specific_monomer(monomer_type = monomer_new, 
                                        linkage_type = linkage_new)
        metrics_P, monomer_count_P, MW_P = ch.get_metrics_polymer(polymer, 
                                                                  additional=self.additional, 
                                                                  cal_MW=self.cal_MW)
        distance_init = ut.cal_distance(self.metrics_target, metrics_P, self.metrics_weights)

        # update the current size
        if self.size_in_MW: current_size = MW_P
        else: current_size = monomer_count_P
        
        '''
        Simulation main 
        '''
        # Trajectory arrays
        d = distance_init
        acceptance_count = 0
        distance = [d]
        start = time.time()
        i_step = 0 # event indices

        while (current_size <= self.stop_size) and (i_step <= i_max):

            Tmetro = self.Tmetro # Annealing version: (1 - monomer_count_P/self.stop_size) * self.Tmetro
            
            polymer_i = Polymer(polymer, verbose=False)# self.verbose)
            # update the branching state based on the branching propensity
            if self.branching_propensity is not None:
                if self.branching_propensity > 0.0: # randomly generate a branching state based on branching propensity
                    branching_state = ut.generate_random_branching_state(self.branching_propensity)
                if self.branching_propensity == 0.0: # no branching
                    branching_state = False 

            linkage_new = ut.generate_random_linkage(self.linkage_distribution, random_state)
            monomer_new = ut.generate_random_monomer(self.monomer_distribution, random_state)

            new_linkage_flag = polymer_i.add_specific_monomer(monomer_type = monomer_new, 
                                                             linkage_type = linkage_new,
                                                             branching_state=branching_state)

            if new_linkage_flag:

                metrics_P, monomer_count_P, MW_P = ch.get_metrics_polymer(polymer, 
                                                                  additional=self.additional, 
                                                                  cal_MW=self.cal_MW)
                # check if the molecule is valid
                if MW_P < 100: continue

                d_new = ut.cal_distance(self.metrics_target, metrics_P, self.metrics_weights)
                energy_flag = False
                delta_d = d_new - d

                # accept the change if energy going downhill
                if delta_d <= 0 or Tmetro == np.inf :
                    energy_flag = True
                # test using Boltzmann distribution
                else:
                    if Tmetro > 0:
                        w = np.exp(-delta_d / kb /Tmetro)
                        if np.random.rand() <= w:
                            energy_flag = True

                if energy_flag: #and distance_flag and NN_flag):

                    d = d_new
                    polymer = polymer_i        
                    acceptance_count += 1

                    # update the current size
                    if self.size_in_MW: current_size = MW_P
                    else: current_size = monomer_count_P
                else:
                    pass
                
                # if self.verbose:
                #     print('\t\t{} monomers after MC step {}'.format(monomer_count_P, i_step+1))
                if self.file is not None:
                    self.file.write('\t\t{} monomers after MC step {}\n'.format(monomer_count_P, i_step+1))
            
            i_step += 1
            distance.append(d)
        
        end = time.time()

        # if self.verbose:
        #     print("\t\tPolymerization takes {0:.1f} min".format((end-start)/60))
        if self.file is not None:
            self.file.write("\t\tPolymerization takes {0:.1f} min\n".format((end-start)/60))
        
        # Extract the polymer object
        polymer_final = polymer

        return polymer_final, distance, monomer_count_P, i_step


    def run_MCMC_ring(
        self, 
        polymer_init: Polymer, 
        rseed: int, 
        i_max_ring: Optional[int] = 100):
        """Run MCMC simulations to add the ring

        Parameters
        ----------
        polymer_init : Polymer
            initial polymer object
        rseed : int
            random seed
        i_max_ring : Optional[int], optional
            maximum number of iterations in the inner loop, 
            by default 500

        Returns
        -------
        polymer_final : Polymer
            final polymer object
        MC_dict : dict
            MC trajectory
        acceptance_count: int
            accepted ring addition
        
        """        
        # Seed the random state
        random_state = ut.set_random_state(rseed)
        # np.random.seed(rseed)

        polymer = Polymer(polymer_init, verbose= self.verbose)
        metrics_P, monomer_count_P, MW_P = ch.get_metrics_polymer(polymer, 
                                                                       additional=self.additional, 
                                                                       cal_MW=self.cal_MW)

        distance_init = ut.cal_distance(self.metrics_target, metrics_P, self.metrics_weights)

        '''
        Simulation main 
        '''
        # Trajectory arrays
        d = distance_init
        distance = [d]
        acceptance_count = 0

        start = time.time()
        i_step = 0 # event indices

        while i_step <= i_max_ring:

            Tmetro = self.Tmetro # Annealing version: (1 - i_step/i_max) * Tmax #
            
            polymer_i = Polymer(polymer, verbose=False)# self.verbose)
            new_linkage_flag = polymer_i.add_specific_ring(linkage_type = ut.generate_random_linkage(self.linkage_distribution, random_state))
            
            
            if new_linkage_flag:
                #P_new = polymer.P
                metrics_P, monomer_count_P, MW_P =   ch.get_metrics_polymer(polymer, 
                                                                  additional=self.additional, 
                                                                  cal_MW=self.cal_MW)
                # check if the molecule is valid
                if MW_P < 100: continue

                d_new = ut.cal_distance(self.metrics_target, metrics_P, self.metrics_weights)
                energy_flag = False
                delta_d = d_new - d

                # accept the change if energy going downhill
                if delta_d <= 0 or Tmetro == np.inf :
                    energy_flag = True
                # test using Boltzmann distribution
                else:
                    if Tmetro > 0:
                        w = np.exp(-delta_d / kb /Tmetro)
                        if np.random.rand() <= w:
                            energy_flag = True

                if energy_flag: #and distance_flag and NN_flag):

                    d = d_new
                    polymer = polymer_i
                    acceptance_count += 1
                else:
                    pass

                # if self.verbose:
                #     print('\t\t{} ring(s) added after {} MC steps'.format(acceptance_count, i_step + 1))
                if self.file is not None:
                    self.file.write('\t\t{} ring(s) added after {} MC steps\n'.format(acceptance_count, i_step + 1))

            distance.append(d)
            i_step += 1
        
        
        end = time.time()
        # if self.verbose:
        #     print("\t\tRing addition takes {0:.1f} min".format((end-start)/60))
        if self.file is not None:
            self.file.write("\t\tRing addition takes {0:.1f} min\n".format((end-start)/60))
        
        # Extract the final structure
        polymer_final = polymer #accepted[-1]

        return polymer_final, distance, acceptance_count, i_step




class Simulation(Trajectory):
    """Full MCMC simulation 
    to generate a population of lignin structures
    """
    def __init__(self, 
        linkage_distribution_input: list, 
        monomer_distribution_input: list, 
        expected_size: Tuple[int, float], 
        max_size: Tuple[int, float], 
        distribution_scaling: float,
        Tmetro: int,
        Tmetro_out: int,
        seed_init: Optional[int] = 1,
        library_name: Optional[str] = 'lignin_x',
        ResultsName: Optional[str] = 'results',
        trial_index: Optional[int] = None,
        n_population: Optional[int] = 100,
        i_max: Optional[int] = 1000,
        i_max_out: Optional[int] = 1000,
        i_max_ring: Optional[int] = 500,
        additional_metrics: Optional[Tuple[list, float]] = None,
        population_metrics: Optional[list] = None,
        size_in_MW: Optional[bool] = False,
        branching_propensity: Optional[float] = None,
        metrics_weights: Optional[Tuple[list, nparray]] = None,
        verbose: Optional[bool] = True,
        show_plots: Optional[bool] = True,
        save_path: Optional[str] = os.getcwd()):
        """Initialize the simulation parameters

        Parameters
        ----------
        linkage_distribution_input : list
            target linkage distribution
        monomer_distribution_input : list
            target monomer distribution
        expected_size : Tuple[int, float]
            target polymer size (mean)
        max_size : Tuple[int, float]
            maximum polymer size
        distribution_scaling : Optional[float], optional
            std of the size distribution, by default 1.0
        Tmetro : int
            temperature for metropolis monte carlo in the inner loop
        Tmetro_out : int
            temperature for metropolis monte carlo in the outer loop
        seed_init : Optional[int], optional
            random seed, by default 1
        library_name : Optional[str], optional
            name of structure library or 
            type of lignin, by default 'lignin_x'
        trial_index : Optional[int], optional
            trial no., by default None
        n_population : Optional[int], optional
            population size, by default 100
        i_max : Optional[int], optional
            maximum number of iterations in the inner loop, by default 1000
        i_max_out : Optional[int], optional
            maximum number of iterations in the outter loop, by default 1000
        i_max_ring : Optional[int], optional
            maximum number of iterations in the ring loop, by default 500
        additional_metrics : Optional[Tuple[list, float]], optional
            additional metrics such as branching coefficient, by default None
        population_metrics : Optional[list], optional
            list of population metrics such as average MW, by default None
        size_in_MW : Optional[bool], optional
            flag indicating whether size is in MW, by default False
        branching_propensity : Optional[float], optional
            branching propensity, by default None
        metrics_weights : Optional[Tuple[list, nparray]], optional
            assigned weights for each metrics, by default None
        verbose : Optional[bool], optional
            flag to control printing outputs, by default True
        show_plots : Optional[bool], optional
            flag to control showing plots, by default True
        save_path : Optional[str], optional
            output file locations, by default os.getcwd()
        """        
        # inherit from Trajectory
        super().__init__(linkage_distribution_input, 
                        monomer_distribution_input,
                        Tmetro,
                        expected_size,
                        max_size,
                        distribution_scaling,
                        additional_metrics,
                        size_in_MW,
                        branching_propensity,
                        metrics_weights, 
                        verbose)
        
        # Set up io path
        ResultsPathParent = os.path.join(save_path, ResultsName, library_name)
        if not os.path.exists(ResultsPathParent): os.mkdir(ResultsPathParent)

        # Determine previous trials in library
        if trial_index==None:
            trial_count = 0
            trial_list=os.scandir(ResultsPathParent)
            for trial_folder in trial_list:
                # print(trial_folder)
                if trial_folder.is_dir():
                    trial_count=trial_count+1
            trial_index=trial_count
                    
        # create a new folder for the new trial            
        if trial_index is None:
            ResultsPath = ResultsPathParent
        else:
            ResultsPath = os.path.join(ResultsPathParent, 'i'+str(trial_index))
        
        print("Starting a new trial, No.{}:\n". format(trial_index))
            
        if os.path.exists(ResultsPath): shutil.rmtree(ResultsPath)
        os.mkdir(ResultsPath)

        # assign to self
        self.linkage_distribution_input = linkage_distribution_input
        self.monomer_distribution_input = monomer_distribution_input
        self.library_name = library_name
        self.n_population = n_population
        self.Tmetro_out = Tmetro_out
        self.seed_init = seed_init
        self.i_max = i_max
        self.i_max_out = i_max_out
        self.i_max_ring = i_max_ring
        self.population_metrics = population_metrics
        self.ResultsName=ResultsName
        self.ResultsPath = ResultsPath
        self.P_population = None
        self.metrics_current_dict= None
        self.metrics_target_dict = None
        self.metrics_names_to_plot = None
        self.metrics_population_to_plot=None
        self.trial_index=trial_index
        self.show_plots = show_plots

        # estimate the maximum MW for scaling purposes 
        if self.size_in_MW:
            self.max_MW = self.max_size
        else:
            self.max_MW = 150 * self.max_size
        self.cal_MW = True

    
    def run(self):
        """
        Main body
        run the MCMC simulations
        """

        # Set the seed for this entire simulation
        np.random.seed(self.seed_init)

        FilePath = os.path.join(self.ResultsPath, self.library_name + '.out')
        file = open(FilePath,"w")
        file.write('==========Generating {} libraray =========\n'.format(self.library_name))
        file.write('==========Current trial No. {} =========\n'.format(self.trial_index))

        # Initialize the similuation trajectory 
        n_polymers = 0
         # Set the distance tolerance
        d_average = np.inf # set to a very large value
        # Starting seed for each trajectory
        rseed = 0 

        # assign individual trajectory weights
        # exclude the last two population metrics from the input metrics
        metrics_weights_individual = None
        if self.metrics_weights is not None:
            if self.population_metrics is not None:
                metrics_weights_individual = self.metrics_weights[:-2]
            else:
                metrics_weights_individual = self.metrics_weights

        traj = Trajectory(self.linkage_distribution_input, 
                        self.monomer_distribution_input, 
                        self.Tmetro, 
                        expected_size = self.expected_size, 
                        max_size = self.max_size,
                        distribution_scaling=self.distribution_scaling,
                        size_in_MW = self.size_in_MW, 
                        additional_metrics = self.additional_metrics,
                        branching_propensity = self.branching_propensity,
                        metrics_weights=metrics_weights_individual,
                        verbose = self.verbose,
                        file=file)

        # Set the optimization target
        metrics_target = traj.metrics_target
        metrics_names_to_match = monomer_types + linkage_names
        if self.additional:
            metrics_names_to_match += ['branching_coeff']
        if self.population_metrics is not None:
            metrics_names_to_match += ['MW', 'MW_weighted']

        # Add population metrics to match, use normalized value for population metrics (MWs)
        metrics_target_to_match = metrics_target.copy()
        metrics_target_to_match_original = metrics_target.copy()
        if self.population_metrics is not None:
            population_metrics_normalized = [pi/self.max_MW for pi in self.population_metrics]
            metrics_target_to_match = np.append(metrics_target_to_match, population_metrics_normalized)
            metrics_target_to_match_original = np.append(metrics_target_to_match_original, self.population_metrics)
            
        self.metrics_target_dict = ut.metrics_array_to_dict(metrics_target_to_match_original, metrics_names_to_match)

        # Simulation main body
        P_population = []
        counts_population = []
        monomer_count_population = []
        MW_population = []
        rseeds = []
        metrics_population = []
        distance = []
        distance_polymer = []
        i_step = 0
        start = time.time()
        monomer_accepted=0
        monomer_iterations=0

        while n_polymers <= self.n_population and i_step <= self.i_max_out:
                
            P_i, distance_i, monomer_accepted_count, monomer_iteration_count = traj.run_MCMC(rseed, self.i_max)
            counts_P, monomer_count_P, MW_P = ch.get_counts_polymer(P_i, \
                additional=self.additional, cal_MW=self.cal_MW)
            # check if the molecule is valid
            if MW_P < 100: continue

            metrics_P = ut.counts_to_metrics(counts_P, additional=self.additional)

            # make a copy
            counts_population_copy = counts_population.copy()
            monomer_count_population_copy = monomer_count_population.copy()
            MW_population_copy = MW_population.copy()
            metrics_population_copy = metrics_population.copy()
            
            # append the new 
            counts_population_copy.append(counts_P)
            monomer_count_population_copy.append(monomer_count_P)
            MW_population_copy.append(MW_P)
            metrics_population_copy.append(metrics_P)
            
            # compute the metrics from the sum of counts
            counts_sum = np.sum(np.array(counts_population_copy), axis = 0)
            metrics_average_to_match = ut.counts_to_metrics(counts_sum, additional=self.additional)              
            
            # add population metrics to match 
            if self.population_metrics is not None:
                MW_average = ut.MW_array_to_number_average(np.array(MW_population_copy))
                MW_weight_average = ut.MW_array_to_weight_average(np.array(MW_population_copy))

                metrics_average_to_match = np.append(metrics_average_to_match, [MW_average/self.max_MW, MW_weight_average/self.max_MW])
            
            # Compute the new distance
            d_average_new = ut.cal_distance(metrics_target_to_match, metrics_average_to_match, self.metrics_weights)
            
            energy_flag = False
            delta_d = d_average_new - d_average

            # accept the change if energy going downhill
            if delta_d <= 0 or self.Tmetro_out == np.inf :
                energy_flag = True
            # test using Boltzmann distribution
            else:
                if self.Tmetro_out > 0:
                    w = np.exp(-delta_d / kb /self.Tmetro_out)
                    if np.random.rand() <= w:
                        energy_flag = True

            if energy_flag: #and distance_flag and NN_flag):
                d_average = d_average_new
                counts_population = counts_population_copy
                monomer_count_population = monomer_count_population_copy
                MW_population = MW_population_copy
                metrics_population = metrics_population_copy
            
                P_population.append(P_i)
                n_polymers += 1
                
                file.write('\tn_polymer {} added on iteration no {} \n'.format(n_polymers-1, i_step))
                if self.verbose:
                    print('\tn_polymer {} added on iteration no {}'.format(n_polymers-1, i_step))
            else:
                file.write('\tPolymer addition rejected\n')
                if self.verbose:
                    print('\tPolymer addition rejected')
            
            distance_polymer.append(distance_i)
            distance.append(d_average)
                
            rseed += 1
            rseeds.append(rseed)
            
            monomer_accepted += monomer_accepted_count
            monomer_iterations += monomer_iteration_count
            
            i_step += 1
            

        distance_final = d_average
        end = time.time()
        file.write('Runtime for creating all polymers : {:.2f} minutes \n'.format((end-start)/60))
        print('Runtime for creating all polymers : {:.2f} minutes \n'.format((end-start)/60))

        # Add rings
        distance_polymer_w_ring = []
        distance_w_ring = np.Inf
        if (self.branching_propensity is None) or (self.branching_propensity > 0.0):
            counts_population_w_ring = []
            monomer_count_population_w_ring = []
            MW_population_w_ring = []
            P_population_w_ring = []
            
            metrics_population_w_ring = []
            
            start = time.time()
            polymer_no=1
            ring_count=0
            ring_iterations=0
            for Pi, ri in zip(P_population, rseeds):
                P_i, distance_i, acceptance_count, ring_its = traj.run_MCMC_ring(Pi, ri, self.i_max_ring)

                counts_w_ring_P, monomer_count_w_ring_P, MW_w_ring_P = ch.get_counts_polymer(P_i, \
                    additional=self.additional, cal_MW=self.cal_MW)
                # check if the molecule is valid
                if MW_w_ring_P < 100: continue

                metrics_w_ring_P = ut.counts_to_metrics(counts_w_ring_P, additional=self.additional)
                file.write('\t{} ring(s) added to polymer {}\n'.format(acceptance_count, polymer_no))
                if(acceptance_count>0) and self.verbose:
                    print('\t{} ring(s) added to polymer {}'.format(acceptance_count, polymer_no))
                
                P_population_w_ring.append(P_i)
                counts_population_w_ring.append(counts_w_ring_P)
                monomer_count_population_w_ring.append(monomer_count_w_ring_P)
                MW_population_w_ring.append(MW_w_ring_P)
                metrics_population_w_ring.append(metrics_w_ring_P)
                distance_polymer_w_ring.append(distance_i)
                
                ring_count=ring_count+acceptance_count
                ring_iterations=ring_iterations+ring_its
                polymer_no=polymer_no+1
            
            # record the simulation time
            end = time.time()
            
            file.write('Runtime for adding the rings : {:.2f} minutes \n'.format((end-start)/60))
            print('Runtime for adding the rings : {:.2f} minutes \n'.format((end-start)/60))
            
            # compute the metrics from the sum of counts
            counts_w_ring_sum = np.sum(np.array(counts_population_w_ring), axis = 0)
            metrics_w_ring_to_match = ut.counts_to_metrics(counts_w_ring_sum, additional=self.additional)    

            # add population metrics to match 
            if self.population_metrics is not None:
                MW_w_ring_average = ut.MW_array_to_number_average(np.array(MW_population_w_ring))
                MW_w_ring_weight_average = ut.MW_array_to_weight_average(np.array(MW_population_w_ring))
                metrics_w_ring_to_match = np.append(metrics_w_ring_to_match, [MW_w_ring_average/self.max_MW, MW_w_ring_weight_average/self.max_MW])
            
            d_w_ring_average = ut.cal_distance(metrics_target_to_match, metrics_w_ring_to_match, self.metrics_weights)
            delta_d = d_w_ring_average - d_average
            distance_w_ring = d_w_ring_average

            if delta_d <= 0:
                P_population = P_population_w_ring
                counts_population = counts_population_w_ring
                monomer_count_population = monomer_count_population_w_ring
                MW_population = MW_population_w_ring
                metrics_population = metrics_population_w_ring
                distance_final = d_average

        # assign trajectory properties to self
        start = time.time()
        self.distance = distance
        self.distance_polymer = distance_polymer
        self.distance_w_ring = distance_w_ring
        self.distance_polymer_w_ring = distance_polymer_w_ring
        self.distance_final = distance_final

        # Characterize the entire population

        # Average the population 
        MW_population = np.expand_dims(MW_population, 1)
        monomer_count_population = np.expand_dims(monomer_count_population, 1)
        metrics_population = np.array(metrics_population)
        
        # Save population data to csv files 
        self.P_population = P_population
        if (self.trial_index==None):
            population = ch.Population(P_population, self.library_name, ResultsName=self.ResultsName)
        else:
            population = ch.Population(P_population, self.library_name, ResultsName=self.ResultsName, TrialIndex=str(self.trial_index))
        population.analyze()

        # the metrics average including branching coeff and MW
        self.metrics_current_dict = population.get_metrics_mean(additional=self.additional)
        
        end = time.time()

        # Generate plots
        if self.show_plots:
            self.monomer_count_population =  monomer_count_population
            # plot monomer counts instead of MW weighted as the later is a single value 
            metrics_population_to_plot = np.concatenate((metrics_population, monomer_count_population, MW_population), axis = 1)        
            metrics_names_to_plot = monomer_types + linkage_names + ['branching_coeff', 'monomer_count', 'MW']

            self.metrics_names_to_plot = metrics_names_to_plot
            self.metrics_population_to_plot = metrics_population_to_plot

            # Plot the metrics distribution
            ut.plot_metrics(self.metrics_target_dict, self.metrics_current_dict, metrics_population_to_plot, metrics_names_to_plot, self.ResultsPath)

            # Plot the population distance  
            ut.plot_distance_trajectory(self.distance, self.library_name, 'Population', self.ResultsPath)

            # Plot the distance trajectory for the last polymer
            p_index = len(self.distance_polymer)-1
            ut.plot_distance_trajectory(self.distance_polymer[p_index], self.library_name + '_polymer_'+str(p_index), 'Pi', self.ResultsPath)


        # final printouts for target and optimal values
        print('Target values:')
        pprint.pprint(self.metrics_target_dict, indent=4)
        print('Optimal simulated values:')
        pprint.pprint(self.metrics_current_dict, indent=4)
        
        print('\nAcceptance Rates')
        print('Monomer Acceptance: {}'.format(monomer_accepted/monomer_iterations))
        print('Polymer Acceptance: {}'.format(n_polymers/(i_step-1)))
        if (self.branching_propensity is None) or (self.branching_propensity > 0.0):
            print('Ring Acceptance: {}'.format(ring_count/ring_iterations))
        
        print('Runtime for analyzing the results : {:.2f} minutes \n'.format((end-start)/60))

        # write to output file
        file.write('Target values:\n')
        pprint.pprint(self.metrics_target_dict, indent=4, stream=file)
        file.write('Optimal simulated values:\n')
        pprint.pprint(self.metrics_current_dict, indent=4, stream=file)
        file.write('\nAcceptance Rates:\n')
        file.write('Monomer Acceptance: {}\n'.format(monomer_accepted/monomer_iterations))
        file.write('Polymer Acceptance: {}\n'.format(n_polymers/(i_step-1)))
        if (self.branching_propensity is None) or (self.branching_propensity > 0.0):
            file.write('Ring Acceptance: {}\n'.format(ring_count/ring_iterations))
        file.write('Runtime for analyzing the results : {:.2f} minutes \n'.format((end-start)/60))
        file.close()


#%% Bayesian Optimizer classes using Nextorch 
class ParameterMask():
    """
    Mask object to create full X matrix with fixed values and varying values
    """
    def __init__(
        self, 
        para_ranges : list, 
        return_1d : Optional[bool] = True):
        """Initalize the mask with input parameter ranges

        Parameters
        ----------
        para_ranges : list
            list of parameter ranges (as lists)
            or parameter values (as float)
        return_1d : Optional[bool], optional
            flag to return 1d vector, by default True
            If False, return a 2d matrix
        """        
        
        # set the return dimension
        self.return_1d = return_1d
        
        # the dimension of inputs
        self.n_dim = len(para_ranges)
        
        self.varying_axes = []
        self.para_ranges_varying = []
        self.para_fixed_values = []
        
        # Count the number of fixed inputs
        for pi, para_range_i in enumerate(para_ranges):
            if isinstance(para_range_i, list):
                self.para_ranges_varying.append(para_range_i)
                self.varying_axes.append(pi)
            else:
                self.para_fixed_values.append(para_range_i)
    
    
    def prepare_X(self, X) -> nparray:
        """Prepare the full X vector
        
        Returns
        -------
        X_full : nparray
            full X vector 
        """       
        #If 1D, make it 2D a matrix
        X_temp = copy.deepcopy(np.array(X))
        if len(X_temp.shape)<2:
            X_temp = np.expand_dims(X_temp, axis=0) #If 1D, make it 2D array
    
        n_points = X_temp.shape[0]
        self.X_temp = X_temp
        xi_list = [] # a list of the columns
        di_fixed = 0 # index for fixed value dimensions
        di_varying = 0 # index for varying x dimensions
    
        for di in range(self.n_dim):
            # Get a column from X_test
            if di in self.varying_axes:
                xi = X_temp[:, di_varying]
                di_varying += 1
            # Create a column of fix values
            else:
                fix_value_i = self.para_fixed_values[di_fixed]
                # for numerical values
                if isinstance(fix_value_i, float) or isinstance(fix_value_i, int):
                    xi = np.ones((n_points, 1)) * fix_value_i
                # for other types of values such as None or str
                else: 
                    xi = np.array([fix_value_i] * n_points)
                di_fixed += 1
    
            xi_list.append(xi)
        # Stack the columns into a matrix
        X_full = np.column_stack(xi_list)
        
        if self.return_1d:
            X_full = np.squeeze(X_full, axis = 0)
        
        return X_full
            
    
class VectorizedFunc():
    """
    Wrapper for the vectorized objective function
    """
    def __init__(self, objective_func: object):
        """Initalize with the objective function object

        Parameters
        ----------
        objective_func : object
            objective function object
        """        
        self.objective_func = objective_func
        
    def predict(self, X_real: nparray) -> nparray:
        """Vectorized objective function object
        takes a matrix and return the reponses in a matrix

        Parameters
        ----------
        X_real : nparray 
            input parameter matrix

        Returns
        -------
        Y_real : nparray
            output reponse matrix
        """      
        if len(X_real.shape) < 2:
            X_real = np.expand_dims(X_real, axis=1) #If 1D, make it 2D array
            
        Y_real = []
        for i, xi in enumerate(X_real):
            yi = self.objective_func(xi)    
            Y_real.append(yi)
                
        Y_real = np.array(Y_real)
        # Put y in a column
        Y_real = np.expand_dims(Y_real, axis=1)
            
        return Y_real 


class MaskedFunc(VectorizedFunc):
    """
    Wrapper for the objective function with fixed and varying inputs
    """
    def __init__(self, objective_func, para_ranges):
        """Initalize with the objective function object and parameter ranges

        Parameters
        ----------
        objective_func : object
            objective function object
        para_ranges : list
            list of parameter ranges (as lists)
            or parameter values (as float)
        """        
        super().__init__(objective_func)
        self.mask = ParameterMask(para_ranges)
        
    def predict(self, X_real : nparray) -> nparray:
        """Vectorized objective function object
        takes a matrix and return the reponses in a matrix

        Parameters
        ----------
        X_real : nparray 
            input parameter matrix

        Returns
        -------
        Y_real : nparray
            output reponse matrix
        """ 
        if len(X_real.shape) < 2:
            X_real = np.expand_dims(X_real, axis=1) #If 1D, make it 2D array
            
        Y_real = []
        for i, xi in enumerate(X_real):
            xi_full = self.mask.prepare_X(xi)
            #print(xi_full)
            yi = self.objective_func(xi_full)    
            Y_real.append(yi)
                
        Y_real = np.array(Y_real)
        # Put y in a column
        Y_real = np.expand_dims(Y_real, axis=1)
            
        return Y_real 
        
        
        
class BOOptimizer():
    """
    Automated BO Optimizer
    """
    def __init__(self, name: Optional[str] = 'optimizer_0'):
        """Initialize the optimizer name

        Parameters
        ----------
        name : str, Optional[str]
            optimizer name, by default 'optimizer_0'
        """
        self.name = name
        
    def optimize(
        self, 
        objective_func : object, 
        para_ranges : list, 
        n_iter : Optional[int] = 100, 
        make_plot : Optional[bool] = True, 
        log_flag : Optional[bool] = False
    ) -> Tuple[nparray, nparray, int]:
        """Training a Bayesian Optimizer

        Parameters
        ----------
        objective_func : object
            objective function object
        para_ranges : list
            list of parameter ranges (as lists)
            or parameter values (as float)
        n_iter : Optional[int], optional
            number of optimization iterations, by default 100
        make_plot : Optional[bool], optional
            flag to make discovery plot, by default True
        log_flag : Optional[bool], optional
            flag to take log of the responses, by default False

        Returns
        -------
        X_opt_full : nparray
            values of optimal parameters
        y_opt : nparray
            value of optimal reponse
        index_opt : int
            index of the sampling points where optimum is
        """
        # Vectorize the objective function
        objective_func_vectorized = MaskedFunc(objective_func, para_ranges)
        para_ranges_varying = objective_func_vectorized.mask.para_ranges_varying
        
        # Initialize a BO experiment
        Exp = bo.Experiment(self.name)
        
        # the dimension of inputs
        n_dim = len(para_ranges_varying)
        
        # Latin hypercube design with 10 initial points
        n_init = 5 * n_dim
        X_init = doe.latin_hypercube(n_dim = n_dim, n_points = n_init, seed= 1)
        
        #print(X_init)
        Y_init = bo.eval_objective_func(X_init, para_ranges_varying, objective_func_vectorized.predict)
        
        # Import the initial data
        Exp.input_data(X_init, Y_init, X_ranges=para_ranges_varying, unit_flag=True)
        Exp.set_optim_specs(objective_func=objective_func_vectorized.predict, maximize=False)
        
        # Run optimization loops        
        Exp.run_trials_auto(n_iter)
        
        # Extract the optima
        y_opt, X_opt, index_opt = Exp.get_optim()
        
        # Expand the X into the full size
        X_opt_full = objective_func_vectorized.mask.prepare_X(X_opt)
        
        if make_plot:
            # Plot the optimum discovered in each trial
            plotting.opt_per_trial_exp(Exp, log_flag=log_flag, save_fig=True)

        return X_opt_full, y_opt, index_opt



        





            
