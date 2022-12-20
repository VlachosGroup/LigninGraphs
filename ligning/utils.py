"""
Utility functions

Graph to smiles, mol conversion, visualization etc.
"""
from typing import Optional, TypeVar, Union, Tuple, List
import os 
import pandas as pd
import re
import platform
import shutil

from collections import Counter
import numpy as np
from numpy.random import RandomState
import networkx as nx
import matplotlib.pyplot as plt
from pysmiles import write_smiles, read_smiles, fill_valence
from rdkit import Chem
from scipy import stats 
from scipy.stats import norm
from scipy.stats import beta

from ligning.rules import monomer_types, monomer_select_C1_C2, linkage_names, CHO, weight_CHO


# Create a type variable for networkx graph
nxgraph = TypeVar('nxgraph')
# Create a type variable for the molecule object in rdkit
molecule = TypeVar('molecule')
# Create a type variable for arrays from numpy
nparray = TypeVar('nparray')

# matplotlib settings
import matplotlib
if platform.system() == 'Linux':
    matplotlib.use('Agg')
    print('switched')
import matplotlib.pyplot as plt   
# if platform.system() == 'Linux':
#     plt.switch_backend('agg')
    
# font = {'size'   : 20}

# matplotlib.rc('font', **font)
# matplotlib.rcParams['axes.linewidth'] = 1.5
# matplotlib.rcParams['xtick.major.size'] = 8
# matplotlib.rcParams['xtick.major.width'] = 2
# matplotlib.rcParams['ytick.major.size'] = 8
# matplotlib.rcParams['figure.dpi'] = 300.


def draw_graph(
    G: nxgraph, 
    node_labels: Optional[dict] = None,
    node_shape: Optional[str] = 'o',
    node_size: Optional[int] = 500
) -> None:
    """Plot the connected graph

    Parameters
    ----------
    G : nxgraph
        a connected graph
    node_labels : Optional[dict], optional
        dictionary of node labels to show, 
        by default None (indices are the labels)
        keys are the indices, values are the the node labels 
    node_shape : Optional[str], optional
        node marker symbol, by default 'o'
    node_size : Optional[int], optional
        node marker size, by default 300
    """    
    plt.figure(figsize = (8,8)) 
    nx.draw_networkx(G, 
                    with_labels=True, 
                    node_color = list(nx.get_node_attributes(G,'color').values()), 
                    labels=node_labels,
                    node_shape=node_shape,
                    node_size=node_size)


def draw_big_graph(G: nxgraph) -> None:
    """Plot the big graph

    Parameters
    ----------
    G : nxgraph
        a connected graph, each node is a monomer
    """    
    node_labels = {}
    for i in range(len(G.nodes)):
        node_labels[i] = G.nodes[i]['mtype']
    # big graph is marked by hexagons 
    draw_graph(G, node_labels, node_shape = 'h', node_size=1000)


def draw_atomic_graph(G: nxgraph) -> None:
    """Plot the atomic graph

    Parameters
    ----------
    G : nxgraph
        a connected graph, each node is an atom
    """   
    node_labels = {}
    for i in range(len(G.nodes)):
        node_labels[i] = G.nodes[i]['element']
    # atomic graph is marked by circles
    draw_graph(G, node_labels)


def graph_to_smile(G: nxgraph) -> str:
    """Convert the graph to smile string

    Parameters
    ----------
    G : nxgraph
        a connected graph, each node is an atom

    Returns
    -------
    smiles : str
        smiles string of the molecule
    """    
    Gcopy = G.copy() # make a copy and do not make changes on the orginal graph
    fill_valence(Gcopy, respect_hcount=True)
    smiles = write_smiles(Gcopy)
    return smiles


def graph_to_mol(
    G: nxgraph, 
    save_mol: Optional[bool]=False,
    name: Optional[str]='test',
    save_path: Optional[str]=os.getcwd()
) -> molecule:
    """Convert the graph to the mol object in RDKit
    and save it as a figure

    Parameters
    ----------
    G : nxgraph
        a connected graph, each node is an atom
    save_mol : Optional[bool], optional
        flag to save the mol as a figure, by default False
    name : Optional[str], optional
        name of the molecule, by default 'test'
    save_path : Optional[str], optional
        path to save the figure, by default os.getcwd()

    Returns
    -------
    mol : molecule
        the molecule object in RDKit
    """      
    smiles = graph_to_smile(G)
    mol = Chem.MolFromSmiles(smiles)

    # save to a local directory 
    if save_mol:
        if not os.path.exists(save_path): os.mkdir(save_path)
        filename = os.path.join(save_path, name+'.png')
        Chem.Draw.MolToFile(mol, filename, size=(500, 500))
    
    return mol
    


def smiles_to_formula(smiles: str) -> str:
    """Convert the smiles to the chemical formula in C, H and O

    Parameters
    ----------
    smiles : str
        smiles string of the molecule

    Returns
    -------
    formula : str
        chemical formula for the molecule
    """    
    # Fill the graph with H atom nodes
    G_with_H = read_smiles(smiles, explicit_hydrogen=True)
    # Count the number of elements
    count_CHO = Counter(dict(G_with_H.nodes(data='element')).values())

    formula = ''
    for ki in CHO:
        formula += ki+ str(count_CHO[ki])
    
    return formula


def formula_to_MW(formula: str) -> float:
    """Calculate molecular weight from the chemical formula

    Parameters
    ----------
    formula : str
        chemical formula for the molecule

    Returns
    -------
    MW : float
        molecular weight
    """    
    # Parse the formula into a list of atom count
    count_CHO = list(map(int, re.findall(r'\d+', formula)))
    
    MW = 0
    for i, ki in enumerate(CHO):
        MW += weight_CHO[ki] * count_CHO[i]

    return MW



def MW_array_to_number_average(MW: nparray) -> float:
    """calculate the number average molecular weight 
    from an array of molecular weight

    Parameters
    ----------
    MW : nparray
        molecular weight of a population

    Returns
    -------
    MW_mean: float
        number average molecular weight
    """    
    MW_mean = np.mean(MW)
    return MW_mean



def MW_array_to_weight_average(MW: nparray) -> float:
    """calculate the number average molecular weight 
    from an array of molecular weight

    Parameters
    ----------
    MW : nparray
        molecular weight of a population

    Returns
    -------
    MW_mean_weight : float 
        weight average molecular weight
    """    
    MW_mean_weight = np.sum(MW**2)/np.sum(MW)
    return MW_mean_weight



def join_two(G1: nxgraph, G2: nxgraph) -> nxgraph:
    """union the two graphs

    Parameters
    ----------
    G1 : nxgraph
        graph 1
    G2 : nxgraph
        graph 2

    Returns
    -------
    union : nxgraph
        the union of two graphs
    """    
    union = nx.disjoint_union(G1,G2)
    return union
    

def select_one_from_many(many: list) -> object:
    """Select a random item from a list

    Parameters
    ----------
    many : list
        list of items

    Returns
    -------
    random_item : object
        a random item
    """      
    n_m = len(many)
    index = np.random.choice(n_m, 1)[0]
    random_item = many[index]
    
    return random_item


def select_one_from_many_with_index(many: list) -> Tuple[object, int]:
    """Select a random item from a list and return the index 

    Parameters
    ----------
    many : list
        list of items

    Returns
    -------
    random_item : object
        a random item
    index : int
        the index
    """    
    n_m = len(many)
    index = np.random.choice(n_m, 1)[0]
    random_item = many[index]

    return random_item, index


def make_available(G: nxgraph, node_index: int) -> nxgraph:
    """make a single C atom available for bonding by changing the bonding property to True

    Parameters
    ----------
    G : nxgraph
        the graph
    node_index : int
        the carbon atom index

    Returns
    -------
    G : nxgraph
        the new graph
    """    
    G.nodes[node_index]['bonding'] = True 
    
    return G


def make_unavailable(G: nxgraph, node_index: int) -> nxgraph:
    """make a single C atom unavailable for bonding by changing the bonding property to False

    Parameters
    ----------
    G : nxgraph
        the graph
    node_index : int
        the carbon atom index

    Returns
    -------
    G : nxgraph
        the new graph
    """    
    G.nodes[node_index]['bonding'] = False 
    
    # change neighboring bond order to 2
    neighoring_edges = list(G.edges(node_index))
    for ei in neighoring_edges:
        G.edges[ei[0], ei[1]]['order'] = 1
    
    return G
    

def make_multi_available(G: nxgraph, monomer_type: str) -> nxgraph:
    """make a multiple C atom available for bonding in a monomer by changing the bonding property to True

    Parameters
    ----------
    G : nxgraph
        the graph
    monomer_type : str
        monomer type, must be 'H', 'G' or 'S'

    Returns
    -------
    G ： nxgraph
        the new graph
 
    Raises
    ------
    ValueError
        Input monomer type not allowed
    """
    if monomer_type not in monomer_types:
        raise ValueError("Monomer type not allowed. Must input H, G or S.")
    # Select the bonding carbon indices
    bonding_Cs = list(monomer_select_C1_C2[monomer_type].keys())
    for Ci in bonding_Cs:
        G = make_available(G, Ci - 1) #adjust for python indexing
        
    return G

def adjust_indices(G: nxgraph) -> nxgraph:
    """Adjust the indices due to node removal

    Parameters
    ----------
    G : nxgraph
        the graph

    Returns
    -------
    Gcopy : nxgraph
        the new graph
    """    
    # index mapping for before and after relabeling 
    mapping = {}
    for ni, ni_new in zip(list(G), range(len(G))): 
        mapping[ni] = ni_new
    
    Gcopy = nx.relabel_nodes(G, mapping)
    return Gcopy


def set_random_state(seed: Optional[int] = None):
    """set random state for a MC trajectory

    Parameters
    ----------
    seed : Optional[int], optional
        random seed, by default None

    Returns
    -------
    random_state : RandomState
        RandomState object in numpy
    """    

    # Set the random state
    if seed is None:
        random_state = np.random.RandomState()
    else:
        random_state = np.random.RandomState(seed)

    return random_state


def generate_random_monomer(
    monomer_distribution: nparray, 
    random_state: Optional[RandomState] = None
) -> str:
    """Generate a monomer type based on target distribution

    monomer_distribution : nparray
        monomer distribution
    random_state : Optional[RandomState], optional
        RandomState object in numpy,
        by default None

    Returns
    -------
    mtype : str
        monomer type
    """
    if random_state is None:
        random_state = np.random
    j = int(np.min(np.nonzero(random_state.rand()<np.cumsum(monomer_distribution))))
    
    mtype = monomer_types[j]
    return mtype
    

def generate_random_linkage(
    linkage_distribution: nparray, 
    random_state: Optional[RandomState] = None
) -> str:
    """Generate a random linkage type based on target distribution

    Parameters
    ----------
    linkage_distribution : nparray
        linkage distribution
    random_state : Optional[RandomState], optional
        RandomState object in numpy,
        by default None

    Returns
    -------
    linkage_name : str
        linkage name
    """        
    if random_state is None:
        random_state = np.random
    j = int(np.min(np.nonzero(random_state.rand()<np.cumsum(linkage_distribution))))
    
    linkage_name = linkage_names[j]
    return linkage_name


def generate_random_branching_state(
    branching_propensity: float, 
    random_state: Optional[RandomState] = None
) -> bool:
    """Generate a random branching state based on the propensity

    Parameters
    ----------
    branching_propensity : float
        branching propensity, the likelihood of forming a branch
    random_state : Optional[RandomState], optional
        RandomState object in numpy,
        by default None

    Returns
    -------
    branching_state : bool
        branching state, 
        True for allowing branching
    """        
    if random_state is None:
        random_state = np.random

    branching_states = [True, False]
    branching_distribution = [branching_propensity, 1-branching_propensity]

    j = int(np.min(np.nonzero(random_state.rand()<np.cumsum(branching_distribution))))
    
    branching_state = branching_states[j]
    return branching_state


def generate_random_size_from_distribution(
    mean_size: Tuple[int, float], 
    max_size: Tuple[int, float], 
    distribution_scaling: Optional[float] = 1.0,
    size_in_MW: Optional[bool] = False,
    random_state: Optional[RandomState] = None
) -> float:
    """Generate a random size from a normal distribution

    Parameters
    ----------
    mean_size : Tuple[int, float]
        the expected mean size,
        either number average molecular weight
        or expected monomer count
    max_size : Tuple[int, float]
        the max size limit,
        either number average molecular weight
        or max monomer count
    distribution_scaling : Optional[float], optional
        factor to control the width of molecular weight distribution
        between 0 and 1,
        by default 1.0
    size_in_MW : Optional[bool], optional
        size is in molecular weight,
        by default False
    random_state : Optional[RandomState], optional
        RandomState object in numpy,
        by default None

    Returns
    -------
    size : float
        random size
    """
    if random_state is None:
            random_state = np.random

    # Seed the random state
    min_size = 2
    if size_in_MW:
        min_size = 260

    lower, upper = min_size, max_size
    mu, sigma = mean_size, (max_size - min_size)/1.96/2 * distribution_scaling  # scale the range back to std
    dist = stats.truncnorm(
        (lower - mu) / sigma, (upper - mu) / sigma, loc=mu, scale=sigma)
    
    size = dist.rvs(1, random_state=random_state) 

    return size

def generate_random_size_from_beta_distribution(
    random_state: Optional[RandomState] = None
) -> float:
    if random_state is None:
            random_state = np.random
    size_in_MW = True
    # if size_in_MW:
    #     min_size = 145
    #     max_size =  26000
    a, b = 2, 7.5
    size = 108900 * beta.rvs(a, b, size=1000, random_state=random_state)

    return size


def plot_single_distribution(
    ypred_all: nparray, 
    ypred: Optional[Tuple[int, float]] = None, 
    yobj: Optional[Tuple[int, float]] = None, 
    metric_name: Optional[str] = 'x',
    save_path: Optional[str] = None):
    """Compare the normal distribution fitted by the predicted data with the object value

    Parameters
    ----------
    ypred_all : nparray
        the array for the predicted value
    ypred : Tuple[int, float]
        current value, by default None
    yobj : Tuple[int, float]
        objective (target) value, by default None
    metric_name : Optional[str], optional
        name of the metric, by default 'x'
    save_path : Optional[str], optional
        path to save the figure, by default None

    """   
    # expand the dimension
    if len(ypred_all.shape) > 1:
        ypred_all = ypred_all[:,0]

    fig, ax = plt.subplots(figsize=(4,4))
    # histogram
    ax.hist(ypred_all, density=1, bins = 10, alpha=0.5, color='steelblue', label = 'Population')
    # single current and objective values
    if ypred is not None:
        ax.axvline(ypred, color = 'C0', linestyle = 'dotted', linewidth = 2, label = 'Current')
    if yobj is not None:
        ax.axvline(yobj, color = 'C1', linestyle = '--', linewidth = 2, label = 'Target')
    # fit to a normal distribution
    mu = np.mean(ypred_all)
    sigma = np.std(ypred_all)
    x_norm = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
    ax.plot(x_norm, norm.pdf(x_norm, mu, sigma) , color='r', label = 'Gaussian fit')
    x_label = metric_name
    if (metric_name in linkage_names) or (metric_name in monomer_types):
        x_label += '(%)'
    ax.set_xlabel(x_label)
    ax.set_ylabel('Normalized Density')
    ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    ax.set_title(r'$\mu$-{:.2} $\sigma$-{:.2}'.format(mu, sigma))
    
    fig_name = 'Distribution_' + metric_name.replace(' ', '')
    if save_path is None:
        save_path = os.getcwd()
    fig.savefig(os.path.join(save_path, fig_name + '.png'), bbox_inches="tight")


def plot_metrics(
    metrics_target: dict, 
    metrics_current: dict,
    metrics_population: nparray, 
    metrics_names: list, 
    save_path: Optional[str] = None):
    """Plot the distribution for each metrics 

    Parameters
    ----------
    metrics_target : dict
        objective (target) values
    metrics_current : dict
        current values
    metrics_population: nparray
        metrics matrix for a population
    metrics_names : list
        names of the metrics for comparison
    save_path : Optional[str], optional
        path to save the figure, by default None
    """        
    for mi, name_i in enumerate(metrics_names):
        
        ypred_all= metrics_population[:, mi]

        yobj = None
        if name_i in metrics_target.keys():
            yobj=metrics_target[name_i]
            if name_i == 'branching_coeff':
                yobj = None 

        ypred = None
        if name_i in metrics_current.keys():
            ypred=metrics_current[name_i]

        plot_single_distribution(ypred_all, ypred, yobj, metric_name=name_i, save_path=save_path)


def plot_distance_trajectory(
    distances: list,
    simulation_name: Optional[str] = 'x',
    distance_name: Optional[str] = None, 
    save_path: Optional[str] = None):
    """Plot the distribution for each metrics 

    Parameters
    ----------
    Distance : list
        Distance values
    simulation_name : list
        names of the simulation
    distance_name : str
        names of the distance, population or individual
    save_path : Optional[str], optional
        path to save the figure, by default None
    """        
    fig, ax = plt.subplots(figsize=(4,4))
    ax.plot(np.arange(len(distances)), distances)
    ax.set_xlabel('Iterations')

    y_label = 'Distance' 
    if distance_name is not None:
        y_label = distance_name + '_' + 'Distance' 
    ax.set_ylabel(y_label)

    fig_name = 'Distance_' + simulation_name.replace(' ', '')
    if save_path is None:
        save_path = os.getcwd()
    fig.savefig(os.path.join(save_path, fig_name + '.png'), bbox_inches="tight")

        

def cal_CI(y: nparray, confidence: Optional[float] = 0.95) -> nparray:
    """Calculate the 95% confidence interval

    Parameters
    ----------
    y : nparray
        data array
    confidence : Optional[float], optional
        confidence level, by default 0.95

    Returns
    -------
    CI : nparray
        CI array
    """     
    CI = stats.t.interval(confidence, len(y)-1, loc=np.mean(y), scale=stats.sem(y))

    # se = stats.sem(y)
    # CI = se * stats.t.ppf((1 + confidence) / 2, len(y)-1)

    return CI


def clean_csv_cache(FilePath: str):
    """Delete the file generated previously

    Parameters
    ----------
    FilePath : str
        the file path
    """    
    if os.path.exists(FilePath): os.remove(FilePath)


def write_output_on_individual(row_data: list, column_names: list, OutputPath: str):
    """write to the output file for an individual structure
    Parameters
    ----------
    row_data : list
        data in a list
    column_names : list
        column names in a list
    OutputPath : str
        output file path

    """    
    row_dict = {}
    # Convert to a dataframe
    for ki, vi in zip(column_names, row_data):
        row_dict[ki] = [vi]
    df = pd.DataFrame(row_dict)
    
    with open(OutputPath, 'a') as f:
        df.to_csv(f, header=f.tell()==0) #, index_label=False)
    

def cal_distance(
    metrics_target : nparray, 
    metrics_P: nparray,
    metrics_weights: Optional[nparray] = None
) -> float:
    """Calculate the distance between two metrics array

    Parameters
    ----------
    metrics_target : nparray
        the metrics array of the target
    metrics_P : nparray
        the metrics array of the polymer
    metrics_weights : nparray
        the weights array of each metrics

    Returns
    -------
    d : float
        distance
    """    
    # extract nonzero metrics 
    # nonzero_indices = np.where(metrics_target > 0.0)[0]
    # normalized_diff = (metrics_P[nonzero_indices] - metrics_target[nonzero_indices])/metrics_target[nonzero_indices]
    n_metrics = len(metrics_target)
    
    if metrics_weights is None:
        metrics_weights = np.ones(n_metrics)
    # convert to numpy array
    if isinstance(metrics_weights, list):
        metrics_weights = np.array(metrics_weights)
    # normalize to sum of 1
    metrics_weights = metrics_weights/np.sum(metrics_weights)

    diff_square = (metrics_P - metrics_target)**2
    d = np.sum(diff_square * metrics_weights)
    
    # print(d)
    # print(metrics_target)
    # print(metrics_P)

    return d


def metrics_array_to_dict(
    metrics_array : nparray,
    metrics_names: List[str]
) -> dict:
    """Convert an array of metrics into dictionary

    Parameters
    ----------
    metrics_array : nparray
        metrics array
    metrics_names : List[str]
        metrics names

    Returns
    -------
    metrics_dict : dict
        metrics names as keys and values in a dictionary
    """
    metrics_dict = {}
    for ci, mi in zip(metrics_names, metrics_array):
        metrics_dict[ci] = mi

    return metrics_dict



def counts_to_metrics(
    counts: nparray, 
    additional: Optional[bool] = False
) -> nparray:
    """Convert counts array to metrics array

    Parameters
    ----------
    counts : nparray
        counts array
    additional : Optional[bool], optional
            include additional metrics, by default False

    Returns
    -------
    metrics : nparray
        metrics array in the same size
    """    
    # separate the monomer and linkage counts
    monomer_distribution_input = counts[:len(monomer_types)]
    linkage_distribution_input = counts[len(monomer_types): len(monomer_types) + len(linkage_names)]        

    # Normalize the distribution to 0-1
    monomer_distribution = np.array(monomer_distribution_input)/np.sum(monomer_distribution_input)
    linkage_distribution = np.array(linkage_distribution_input)/np.sum(linkage_distribution_input)
    
    #prevent bad divsion
    if np.isnan(np.sum(linkage_distribution)):
        linkage_distribution = np.zeros_like(linkage_distribution)

    # Concatenate the metrics
    metrics = np.concatenate((monomer_distribution, linkage_distribution), axis= 0)

    # to include additional features - the branching coefficient = n_branched/monomer_count
    if additional:
        metrics_additional = counts[len(monomer_types) + len(linkage_names):]/np.sum(monomer_distribution_input)
        metrics = np.concatenate((metrics, metrics_additional), axis= 0)

    return metrics

