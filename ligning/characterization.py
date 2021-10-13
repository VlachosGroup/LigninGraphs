"""
Charaterization funcions 
"""
import pprint
pp = pprint.PrettyPrinter(indent=4)
from typing import Optional, TypeVar, Union, Tuple, List
import os
import sys
import shutil

import pandas as pd
import numpy as np
from collections import Counter

#ignore divide by zero warinings
import warnings
warnings.filterwarnings("ignore", message="divide by zero encountered in divide")
warnings.filterwarnings("ignore", message="divide by zero encountered") 
warnings.filterwarnings("ignore", message="invalid value encountered")

from ligning.rules import linkage_names, monomer_types
from ligning.utils import formula_to_MW, graph_to_smile, nxgraph, molecule, nparray, graph_to_mol, smiles_to_formula
import ligning.utils as ut
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem import MolFromSmiles
from ligning.polymer import Polymer


# For each monomer
# H - 11 nodes, no OCH3 group
# G - 13 nodes, 1 OCH3 group
# S - 15 nodes, 2 OCH3 group

def get_metrics_polymer(
    P: Polymer, 
    additional: Optional[bool] = False, 
    cal_MW: Optional[bool] = False
) -> Tuple[nparray, int, float]:
    """Get the metrics and count for a polymer

    Parameters
    ----------
    P : Polymer
        polymer object
    additional : Optional[bool], optional
            include additional metrics, by default False
    cal_MW : bool, optional
            flag to calculate molecular weight, by default False

    Returns
    -------
    metrics_P : nparray
        the metrics array
    monomer_count : int
        the number of monomers
    MW : float
        the molecular weight
    """    
    MW = None
    ch_polymer = Characterize(P)
    ch_polymer.get_metrics(cal_MW=cal_MW, additional=additional)
    metrics_P = ch_polymer.metrics
    monomer_count = ch_polymer.monomer_count

    if cal_MW:
        MW = ch_polymer.MW
        
    return metrics_P, monomer_count, MW



def get_counts_polymer(
    P: Polymer, 
    additional: Optional[bool] = False, 
    cal_MW: Optional[bool] = False
) -> Tuple[nparray, int, float]:
    """Get the counts array for a polymer

    Parameters
    ----------
    P : Polymer
        polymer object
    additional : Optional[bool], optional
            include additional metrics, by default False
    cal_MW : bool, optional
            flag to calculate molecular weight, by default False

    Returns
    -------
    counts_P : nparray
        the counts array
    monomer_count : int
        the number of monomers
    MW : float
        the molecular weight
    """    
    MW = None
    ch_polymer = Characterize(P)
    ch_polymer.get_counts(cal_MW=cal_MW, additional=additional)
    counts_P = ch_polymer.counts
    monomer_count = ch_polymer.monomer_count

    if cal_MW:
        MW = ch_polymer.MW
        
    return counts_P, monomer_count, MW



class CharacterizeGraph():
    """
    polymer charaterization object
    """    
    def __init__(self, G: nxgraph):
        """Initialize with a polymer graph

        Parameters
        ----------
        G : nxgraph
            polymer graph
        """        
        self.G = G
        self.Mol = None 

        # Initialize other properties 
        self.mtype_count = None
        self.monomer_count = None
        self.linkages_count = None
        self.OCH3_count = None
        self.OH_count = None
        self.MW = None
        self.metrics = None
        self.mol = None
        self.smiles = None


    def count_types(self) -> dict:
        """count the monomer types

        Returns
        -------
        mtype_count : dict
            the monomer type counts
        """        
        mtypes = [self.G.nodes[ni]['mtype'] for ni in list(self.G.nodes)]
        
        mtype_count = {'H': mtypes.count('H')/11, # divide by the number of non-H atoms
                       'G': mtypes.count('G')/13,
                       'S': mtypes.count('S')/15}
        self.mtype_count = mtype_count

        return self.mtype_count
    

    def count_monomers(self) -> float:
        """count the total number of monomers, i.e. the polymer size

        Returns
        -------
        monomer_count : float
            the polymer size
        """       
        self.monomer_count = np.sum(list(self.mtype_count.values()))
        return self.monomer_count


    def count_linkages(self) -> dict:
        """count the type of linkages

        Returns
        -------
        linkages_count: dict
            the linkage type counts
        """        
        bonds = [self.G.edges[ei]['btype'] for ei in list(self.G.edges)]

        linkages_count = {}
        for linkage_name_i in linkage_names:
            linkages_count[linkage_name_i] = float(bonds.count(linkage_name_i)) #the bonds under the same name
        
        # adjust for double counting
        linkages_count['4-O-5'] = linkages_count['4-O-5'] /2
        linkages_count['beta-O-4'] = linkages_count['beta-O-4'] /2
        linkages_count['beta-5'] = linkages_count['beta-5'] /2
        linkages_count['beta-beta'] = linkages_count['beta-beta'] /3
        
        self.linkages_count = linkages_count
        
        return self.linkages_count
    

    def count_OCH3(self) -> float:
        """count the number of -OCH3 group

        Returns
        -------
        OCH3_count : float
            -OCH3 counts
        """        
        groups = [self.G.nodes[ni]['group'] for ni in list(self.G.nodes)]
        
        OCH3_count = float(groups.count('OCH3'))/2 # divid by 2 since C and O are both marked as 'OCH3'
        self.OCH3_count = OCH3_count
        
        return self.OCH3_count
    

    def count_OH(self) -> float:
        """count the number of available -OH3 group

        Returns
        -------
        OH_count : float
            -OH counts
        """        
        groups = [self.G.nodes[ni]['group'] for ni in list(self.G.nodes)]
        
        OH_4_indices =  [ni for ni, gi in enumerate(groups) if gi == '4OH']
        OH_9_indices =  [ni for ni, gi in enumerate(groups) if gi == '9OH']
        
        OH_4_available = [ni for ni in OH_4_indices if len(list(self.G.neighbors(ni))) == 1]
        OH_9_available = [ni for ni in OH_9_indices if len(list(self.G.neighbors(ni))) == 1]
        
        self.OH_count = float(len(OH_4_available) + len(OH_9_available))
        
        return self.OH_count
    

    def cal_MW(self) -> float:
        """calculate the molecular weight

        Returns
        -------
        float
            Molecular weight of the polymer
        """    
        self.smiles = graph_to_smile(self.G)    
        #self.Mol = MolFromSmiles(self.smiles)
        #self.MW = ExactMolWt(self.Mol)
        self.formula = smiles_to_formula(self.smiles)
        self.MW = formula_to_MW(self.formula)

        return self.MW


    def cal_all(self, cal_MW = False, print_flag = True):
        """Main count function

        Parameters
        ----------
        cal_MW : bool, optional
            flag to calculate molecular weight, by default False
        print_flag : bool, optional
            flag to print all properties, by default True
        """        
        # Count the types of monomers
        self.count_types()
        # Count the total number of monomers
        self.count_monomers()
        # Count the type of bonds
        self.count_linkages()
        # Count the number of -OCH3
        self.count_OCH3()
        # Count the number of available -OH
        self.count_OH()
        # Calculate the molecular weight
        if cal_MW:
            self.cal_MW()
        
        # Print the output as a dictionary
        if print_flag:
            pp.pprint(vars(self))
    
    
    def cal_metrics(self, cal_MW: Optional[bool] = False) -> nparray:
        """calculate the 10 by 1 metrics array 

        Returns
        -------
        metrics : nparray
            metrics array
        """        
        self.cal_all(cal_MW, print_flag = False)

        monomer_distribution_input = list(self.mtype_count.values())
        linkage_distribution_input = list(self.linkages_count.values())
    
        # Normalize the distribution to 0-1
        monomer_distribution = np.array(monomer_distribution_input)/np.sum(monomer_distribution_input)
        linkage_distribution = np.array(linkage_distribution_input)/np.sum(linkage_distribution_input)
        #print(linkage_distribution_input)
        #print(linkage_distribution)
        
        #prevent bad divsion
        if np.isnan(np.sum(linkage_distribution)):
            linkage_distribution = np.zeros_like(linkage_distribution)

        # Concatenate the metrics
        self.metrics = np.concatenate((monomer_distribution, linkage_distribution), axis= 0)


class Characterize(CharacterizeGraph):
    """
    polymer characterization object
    """
    def __init__(self, P: Polymer):
        """Initialize with a polymer object

        Parameters
        ----------
        P : Polymer
            polymer object
        """     
        super().__init__(P.G)
        self.bigG = P.bigG

        self.connections_count = None
        self.branching_coeff = None


    def count_types(self) -> dict:
        """count the monomer types

        Returns
        -------
        mtype_count : dict
            the monomer type counts
        """ 
        mtypes = [self.bigG.nodes[ni]['mtype'] for ni in list(self.bigG.nodes)]
        
        mtype_count = {'H': mtypes.count('H'),
                       'G': mtypes.count('G'),
                       'S': mtypes.count('S')}
        self.mtype_count = mtype_count

        return self.mtype_count


    def count_monomers(self) -> int:
        """count the total number of monomers, i.e. the polymer size

        Returns
        -------
        monomer_count : float
            the polymer size
        """       
        self.monomer_count = len(self.bigG)
        return self.monomer_count


    def count_linkages(self) -> dict:
        """count the type of linkages

        Returns
        -------
        linkages_count: dict
            the linkage type counts
        """     
        bonds = [self.bigG.edges[ei]['btype'] for ei in list(self.bigG.edges)]

        linkages_count = {}
        for linkage_name_i in linkage_names:
            linkages_count[linkage_name_i] = bonds.count(linkage_name_i) #the bonds under the same name
        
        self.linkages_count = linkages_count
        
        return self.linkages_count


    def count_OCH3(self) -> float:
        return super().count_OCH3()

    def count_OH(self) -> float:
        return super().count_OH()

    def cal_MW(self) -> float:
        return super().cal_MW()
    
    def count_connections(self) -> dict:
        """Count the number of connections of each monomer

        Returns
        -------
        connections_count: dict
            the number of connections
        """        

        connections = [self.bigG.degree(ni) for ni in list(self.bigG.nodes)]
        self.connections_count = dict(Counter(connections))

        return self.connections_count


    def cal_branching(self) -> float:
        """calculate the branching coefficient

        Returns
        -------
        branching_coeff : float
            the branching coefficient: 
            ratio of branched monomers to total monomers
        """        
        if self.connections_count is None:
            self.count_connections()
        if self.monomer_count is None:
            self.count_monomers()
            
        n_branched = 0
        # Count the number of branched monomers
        # a branched monomer is one that is bonded to three or more monomers
        for ki, vi in self.connections_count.items():
            if ki >= 3:
                n_branched += vi

        # number of branched monomers        
        self.n_branched = n_branched
        self.branching_coeff = n_branched/self.monomer_count

        return self.branching_coeff
    

    def cal_all(self, cal_MW = False, print_flag = True):
        """Main count function

        Parameters
        ----------
        cal_MW : bool, optional
            flag to calculate molecular weight, by default False
        print_flag : bool, optional
            flag to print all properties, by default True
        """        
        # Count the types of monomers
        self.count_types()
        # Count the total number of monomers
        self.count_monomers()
        # Count the type of bonds
        self.count_linkages()
        # Count the number of -OCH3
        self.count_OCH3()
        # Count the number of available -OH
        self.count_OH()
        # Count the connections
        self.count_connections()
        # Calculate the branching coefficient
        self.cal_branching()
        # Calculate the molecular weight
        if cal_MW:
            self.cal_MW()
        
        # Print the output as a dictionary
        if print_flag:
            pp.pprint(vars(self))
    
    
    def get_metrics(self, additional: Optional[bool] = False, cal_MW: Optional[bool] = False) -> nparray:
        """Get the metrics array for a polymer

        Parameters
        ----------
        additional : Optional[bool], optional
            include additional metrics, by default False

        Returns
        -------
        metrics : nparray
            metrics array
        """        
        self.cal_all(cal_MW, print_flag = False)

        monomer_distribution_input = list(self.mtype_count.values())
        linkage_distribution_input = list(self.linkages_count.values())
    
        # Normalize the distribution to 0-1
        monomer_distribution = np.array(monomer_distribution_input)/np.sum(monomer_distribution_input)
        linkage_distribution = np.array(linkage_distribution_input)/np.sum(linkage_distribution_input)
        #print(linkage_distribution_input)
        #print(linkage_distribution)
        
        #prevent bad divsion
        if np.isnan(np.sum(linkage_distribution)):
            linkage_distribution = np.zeros_like(linkage_distribution)

        # Concatenate the metrics
        self.metrics = np.concatenate((monomer_distribution, linkage_distribution), axis= 0)

        # to include additional features - the branching coefficient
        if additional:
            metrics_additional = np.array([self.branching_coeff])
            self.metrics = np.concatenate((self.metrics, metrics_additional), axis= 0)

        return self.metrics
    

    def get_counts(self, additional: Optional[bool] = False, cal_MW: Optional[bool] = False) -> nparray:
        """Get the count array for a polymer

        Parameters
        ----------
        additional : Optional[bool], optional
            include additional metrics, by default False

        Returns
        -------
        counts : nparray
            counts array
        """        
        self.cal_all(cal_MW, print_flag = False)

        monomer_counts = list(self.mtype_count.values())
        linkages_counts = list(self.linkages_count.values())

        # Concatenate the counts
        self.counts = np.concatenate((monomer_counts, linkages_counts), axis= 0)

        # to include additional features - the branching coefficient
        if additional:
            counts_additional = np.array([self.n_branched])
            self.counts = np.concatenate((self.counts, counts_additional), axis= 0)

        return self.counts
        
    
    
    
class Population():
    """
    Characterize a population of polymers
    """    
    def __init__(
        self, 
        population: List[Polymer],
        name: Optional[str]='lignin_x', 
        InputPath: Optional[str] =os.path.abspath(os.getcwd()), 
        ResultsName: Optional[str]='results',
        TrialIndex: Optional[str]=None):
        """initalize a population

        Parameters
        ----------
        population : List[Polymer]
            a population of polymer objects
        name : Optional[str], optional
            name of the population, by default 'lignin_x'
        InputPath : str, optional
            the input path, by default os.path.abspath(os.getcwd())
        ResultsName : str, optional
            results folder name, by default 'results'
        """        
        self.population = population
        self.name = name
        # Set the directory structures
        OutputLibrary = name + '_libray.csv'
        OutputStats = name + '_population_stats.csv'
        if TrialIndex==None:
            ResultsPath = os.path.join(InputPath, ResultsName, name)
        else:
            ResultsPath = os.path.join(InputPath, ResultsName, name, "i"+TrialIndex)


        self.OutputPathLibrary = os.path.join(ResultsPath, OutputLibrary)
        self.OutputPathStats = os.path.join(ResultsPath, OutputStats)
        
        # Set the columns names for data
        metrics_names = ['branching_coeff', 'MW', 'monomer_count', 'OH_count', 'OCH3_count',  'smiles']
        self.column_names = monomer_types + linkage_names + metrics_names
        
        # Initialize the stats and data
        self.characterization_objects = []
        self.stats = None
        self.data = None
    

    def characterize_all(self):
        """
        Characterize each individual polymer
        """        
        if len(self.characterization_objects) == 0:
            for polymer_i in self.population: 
                ch_i = Characterize(polymer_i)
                ch_i.cal_all(print_flag=False, cal_MW=True)
                self.characterization_objects.append(ch_i)
        else: 
            pass


    def analyze(self):
        """
        Analyze the population
        Output the data and stats both in csv files
        """

        # delete the previous csv files 
        ut.clean_csv_cache(self.OutputPathLibrary)
        ut.clean_csv_cache(self.OutputPathStats)
        
        self.characterize_all()

        for ch_i in self.characterization_objects: 
            
            row_i = []
            row_i += list(ch_i.mtype_count.values())
            row_i += list(ch_i.linkages_count.values())
            row_i.append(ch_i.branching_coeff)
            row_i.append(ch_i.MW)

            row_i.append(ch_i.monomer_count)
            row_i.append(ch_i.OH_count)
            row_i.append(ch_i.OCH3_count)
            row_i.append(ch_i.smiles)
            
            # write to the data output file
            ut.write_output_on_individual(row_i, self.column_names, self.OutputPathLibrary)

        # Calculate the population stats
        population_data = pd.read_csv(self.OutputPathLibrary)   
        numerical_metrics = list(population_data.mean().index)
        population_mean = np.array([population_data.mean()])
        population_std = np.array([population_data.std()])
        population_CI = np.array([ut.cal_CI(np.array(population_data[ci])) for ci in numerical_metrics])

        stats = np.concatenate((population_mean, population_std, population_CI.T), axis=0)
        pd_stats = pd.DataFrame(stats, index = ['mean', 'std', 'CI_lower', 'CI_upper'], columns=numerical_metrics)
        pd_stats.to_csv(self.OutputPathStats) #,index_label=False)
        
        # update the self
        self.stats = pd_stats
        self.data = population_data


    def get_counts(self, additional: Optional[bool] = False):
        """Get the metrics matrix for the entire population

        Parameters
        ----------
        additional : Optional[bool], optional
            include additional metrics, by default False

        Returns
        -------
        metrics : nparray
            metrics matrix
        """       
        self.characterize_all()
        counts_matrix = []
        for ch_i in self.characterization_objects:
            counts_matrix.append(list(ch_i.get_counts(additional=additional)))
        
        self.counts = np.array(counts_matrix)
        self.counts_sum = np.sum(self.counts, axis = 0)
        
        return self.counts 


    def get_metrics_mean(self, additional: Optional[bool] = False):
        """Get the metrics matrix for the entire population

        Parameters
        ----------
        additional : Optional[bool], optional
            include additional metrics, by default False

        Returns
        -------
        metrics : nparray
            metrics mean array
        """     
        self.get_counts(additional)

        MW = []
        monomer_counts = []
        # if the MW and monomer count have not been calculated 
        if self.data is None:
            for ch_i in self.characterization_objects: 
                MW.append(ch_i.MW)
                monomer_counts.append(ch_i.monomer_count)
        # Extract MW and monomer count directly from the data
        else:
            MW = self.data['MW']
            monomer_counts = self.data['monomer_count']
    
        MW = np.array(MW)
        monomer_counts = np.array(monomer_counts)

        # Compute number and weight average
        self.number_average_MW = ut.MW_array_to_number_average(MW)
        self.weight_average_MW = ut.MW_array_to_weight_average(MW)
        self.monomer_counts_average = np.mean(monomer_counts)

        # Construct an array for the means
        metric_mean = list(ut.counts_to_metrics(self.counts_sum, additional= additional))
        metric_mean += [self.number_average_MW, self.weight_average_MW, self.monomer_counts_average]

        self.metric_mean = np.array(metric_mean)

        self.metrics_mean_dict = {}
        column_names_population = monomer_types + linkage_names
        if additional:
            column_names_population +=  ['branching_coeff']
        column_names_population += ['MW', 'MW_weighted', 'monomer_count']
        
        for ci, mi in zip(column_names_population, self.metric_mean):
            self.metrics_mean_dict[ci] = mi

        return self.metrics_mean_dict






 










    
    

    


    
    
    
    



