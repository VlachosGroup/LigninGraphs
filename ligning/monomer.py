"""
Build lignin monomers
"""
from typing import Optional, TypeVar, Union, Tuple, List

import networkx as nx

from ligning.rules import monomer_types, default_color
from ligning.utils import nxgraph, molecule, nparray
import ligning.utils as ut


def add_OCH3(
    G:nxgraph, monomer_type: str,  
    C_index: int, 
    color: str,
    monomer_index: Optional[int] = 0
) -> nxgraph:
    """Add an OCH3 group to a monomer

    Parameters
    ----------
    G : nxgraph
        the monomer graph
    monomer_type : str
        monomer type, must be 'H' or 'G'
    C_index : int
        the bonding Carbon atom index
    color : str
        the color of the nodes in plotting
    monomer_index: Optional[str], optional
        the index of monomer in a polymer, by default 0

    Returns
    -------
    G : nxgraph
        the monomer graph
    
    Raises
    ------
    ValueError
        Input monomer type not allowed

    """
    if monomer_type not in monomer_types:
        raise ValueError("Monomer type not allowed. Must input H, G or S.")
    # if monomer_type == 'S':
    #     raise ValueError("S has two -OCH3 groups. Cannot add more.")

    O_in_OCH3_index = len(G) 
    C_in_OCH3_index = len(G) + 1
    # add O
    G.add_node(O_in_OCH3_index, element = 'O', aromatic = False, group = 'OCH3', index = O_in_OCH3_index, mtype = monomer_type, bonding = False, color = color, mi=monomer_index)
    # add C
    G.add_node(C_in_OCH3_index, element = 'C', aromatic = False, group = 'OCH3', index = C_in_OCH3_index, mtype = monomer_type, bonding = False, color = color, mi=monomer_index)

    # add C-O bond
    new_edge_1 = (C_index, O_in_OCH3_index)
    G.add_edges_from([new_edge_1], order = 1, index = new_edge_1, mtype = (monomer_type, monomer_type), btype = None)
    
    new_edge_2 = (O_in_OCH3_index, C_in_OCH3_index)
    G.add_edges_from([new_edge_2], order = 1, index = new_edge_2, mtype = (monomer_type, monomer_type), btype = None)

    return G


def monomer_graph(
    monomer_type: str, 
    color: Optional[str]  = default_color['H'],
    monomer_index: Optional[int] = 0
) -> nxgraph:
    """Generate a general monomer with 9 carbons

    Parameters
    ----------
    monomer_type : str
        monomer type, must be 'H', 'G' or 'S'
    color : Optional[str], optional
        the color of the nodes in plotting, by default 'r'
    monomer_index: Optional[str], optional
        the index of monomer in a polymer, by default 0

    Returns
    -------
    G : nxgraph
        the monomer graph

    Raises
    ------
    ValueError
        Input monomer type not allowed
    """    
    if monomer_type not in monomer_types:
        raise ValueError("Monomer type not allowed. Must input H, G or S.")

    # Initialize a networkx graph
    G = nx.Graph()
    # Add 1 to 9 carbons with 7 for alpha, 8 for beta, 9 for gamma
    n_c = 9
    for i in range(1,n_c+1): # -1 for python indexing
        aromatic_flag = False    

        if i in range(1,6+1): 
            aromatic_flag = True    
        G.add_node(i -1 , element = 'C', aromatic = aromatic_flag, group = None, index = i, mtype = monomer_type, bonding = False, color = color, mi=monomer_index)
    
    # Add O attached to 4C (index = 10) - bonding OH group
    O_1_index = n_c + 1
    G.add_node(O_1_index -1, element = 'O', aromatic = False, group = '4OH', index = O_1_index, mtype = monomer_type, bonding = True, color = color, mi=monomer_index)
    
    # Add O attached to 9C (index = 11) - nonbonding OH group
    O_2_index = 11
    G.add_node(O_2_index -1, element = 'O', aromatic = False, group = '9OH', index = O_2_index, mtype = monomer_type, bonding = False, color = color, mi=monomer_index)
    
    # make the bonding nodes available ;;;;;
    G = ut.make_multi_available(G, monomer_type)
    
    # Add edge list one by one
    # Add single bond
    edge_list_single = [(1,2), (2,3), (3,4), (4,5), (5,6), (6,1), (1,7), (7,8), (8,9), (10, 4), (9, 11)]
    for ei in edge_list_single:
        ei_index = (ei[0] -1, ei[1] -1) #adjust for python indexing
        G.add_edges_from([ei_index], order = 1, index = ei, mtype = (monomer_type, monomer_type), btype = None)
    
    # Add double bond
    edge_list_double = [(7,8)]
    for ei in edge_list_double:
        ei_index = (ei[0] -1, ei[1] -1) #adjust for python indexing
        G.add_edges_from([ei_index], order = 2, index = ei, mtype = (monomer_type, monomer_type), btype = None)
    
    return G


def monomer_H(
    color: Optional[str]  = default_color['H'],
    monomer_index: Optional[int] = 0
) -> nxgraph:
    """Generate a H monomer

    Parameters
    ----------
    color : Optional[str], optional
        the color of the nodes in plottin, by default 'r' (red)
    monomer_index: Optional[str], optional
        the index of monomer in a polymer, by default 0

    Returns
    -------
    G : nxgraph
        the monomer graph
    """    
    monomer_type = 'H'
    G = monomer_graph(monomer_type, color, monomer_index)

    return G
    

def monomer_G(
    color: Optional[str]  = default_color['G'],
    monomer_index: Optional[int] = 0
) -> nxgraph:
    """Generate a G monomer

    Parameters
    ----------
    color : Optional[str], optional
        the color of the nodes in plottin, by default 'g' (green)
    monomer_index: Optional[str], optional
        the index of monomer in a polymer, by default 0

    Returns
    -------
    G : nxgraph
        the monomer graph
    """    
    monomer_type = 'G'
    G = monomer_graph(monomer_type, color, monomer_index)
    # Add an additional OCH3 group
    C_index = 2 # the carbon bonding to OCH3
    G = add_OCH3(G, monomer_type, C_index, color, monomer_index)

    return G


def monomer_S(
    color: Optional[str]  = default_color['S'],
    monomer_index: Optional[int] = 0
) -> nxgraph:
    """Generate an S monomer

    Parameters
    ----------
    color : Optional[str], optional
        the color of the nodes in plottin, by default 'b' (blue)
    monomer_index: Optional[str], optional
        the index of monomer in a polymer, by default 0

    Returns
    -------
    G : nxgraph
        the monomer graph
    """  
    monomer_type = 'S'
    G = monomer_graph(monomer_type, color, monomer_index)
    
    # Add two additional OCH3 group
    C1_index = 2 # the carbon bonding to OCH3
    C2_index = 4 # the carbon bonding to OCH3
    G = add_OCH3(G, monomer_type, C1_index, color, monomer_index)
    G = add_OCH3(G, monomer_type, C2_index, color, monomer_index)

    return G


def select_random_monomer() -> nxgraph:
    """Generate a random monomer for initialization purpose

    Returns
    -------
    P : nxgraph
        the initial polymer graph 
    """     
    monomer_type = ut.select_one_from_many(monomer_types)
    P = Monomer(monomer_type)

    return P


class Monomer():
    """
    Lignin monomer object
    """
    def __init__( 
        self, 
        monomer_type: str, 
        monomer_index: Optional[int] = 0
    ):
        """Initialize a monomer graph from the three types: H, G, S

        Parameters
        ----------
        monomer_type : str
            monomer type, must be 'H', 'G' or 'S'
        monomer_index: Optional[str], optional
            the index of monomer in a polymer, by default 0

        """      
        self.G = None
        self.bigG = None

        self.type = monomer_type
        self.mi = monomer_index


    def create(self) -> nxgraph:
        """Initialize a monomer graph from the three types: H, G, S

        Returns
        -------
        G : nxgraph
            the monomer graph

        Raises
        ------
        ValueError
            Input monomer type not allowed
        """      
        if self.type not in monomer_types:
            raise ValueError("Monomer type not allowed. Must input H, G or S.")

        if self.type == 'H': G_m = monomer_H(monomer_index=self.mi)
        if self.type == 'G': G_m = monomer_G(monomer_index=self.mi)
        if self.type == 'S': G_m = monomer_S(monomer_index=self.mi)

        # Update the monomer graph
        self.G = G_m

        # Initialize a big graph for monomer
        bigG = nx.Graph()
        bigG.add_node(self.mi, mtype = self.type, color = default_color[self.type])
        self.bigG = bigG
        
        return G_m