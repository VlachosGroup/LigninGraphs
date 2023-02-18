"""
Build lignin polymers
"""
import networkx as nx

from typing import Optional, TypeVar, Union, Tuple, List
from copy import copy
import warnings

from ligning.rules import linkage_special_names, monomer_select_C1_C2, linkage_index_select_monomer,\
    linkage_index_to_name, linkage_name_select_C1_C2, linkage_ring
from ligning.utils import nxgraph, molecule, nparray
import ligning.utils as ut
from ligning.monomer import Monomer

class PolymerGraph():
    """
    PolymerGraph object contains all graph operations for a polymer
    """
    def __init__(self, G: nxgraph, verbose: Optional[bool] = True):
        """initialize the object with a graph

        Parameters
        ----------
        M_init : Monomer
            object of the initial monomer
        verbose : Optional[bool], optional
            the flag to control on-off for print statement, by default True
        """        
        self.G = G.copy()
        self.verbose = verbose
        # list of available C1 indices for forming linkages
        self.C1_indices_in_polymer = None
    

    def find_available_C1_in_polymer(self, branching_state: Optional[bool] = None) -> list:    
        """find available C1 node in a polymer

        Returns
        -------
        C1_indices_in_polymer : list
            the C1 node indices in the polymer
        branching_state : Optional[bool], optional
            flag to control whether the new C is branced or not
            by default None, select all Cs
            if true, select Cs in the chain
            if false, select Cs at the terminal monomers
        """            
        # check which nodes have "bonding" == True
        C1_indices_in_polymer = [n for n,v in self.G.nodes(data=True) if v['bonding']]

        if branching_state is not None:
            C1_indices_in_polymer = self.update_branching_C_in_polymer(C1_indices_in_polymer, branching_state)
        
        return C1_indices_in_polymer
    

    def update_available_C1_in_polymer(self, C1_index_in_polymer: int) -> list:
        """update available C1 node in a polymer by delete C1_index_in_polymer from the list

        Parameters
        ----------
        C1_index_in_polymer : int
            C1 indices in polymer

        Returns
        -------
        C1_indices_in_polymer : list
            the C1 node indices in the polymer
        """        
        # update the node list by taking out the one at the new linkage
        C1_indices_in_polymer = [n for n in self.C1_indices_in_polymer if not n == C1_index_in_polymer]
        
        return C1_indices_in_polymer
        

    def update_branching_C_in_polymer(self, C_indices_in_polymer: list, branching_state: Optional[bool] = False) -> list:
        """update the C indices based on branching 

        Parameters
        ----------
        C_indices_in_polymer : list
            C indices available for bonding
        branching_state : Optional[bool], optional
            find branching_state or non-branching_state monomers, by default False

        Returns
        -------
        C_available_indices_in_polymer : list
            the C indices in branching_state or non-branching_state monomers
        """        
        terminal_monomers_indices = [] 
        C_available_indices_in_polymer = [] # the results

        # find the monomer indices at the end of the chain
        for mi, node_i in enumerate(self.bigG.nodes):
            if len(self.bigG) < 3:  # for monomer and dimers, we ignore the branching 
                return C_indices_in_polymer 

            if self.bigG.degree(node_i) == 1: 
                terminal_monomers_indices.append(mi)

        # find the intersection between desired monomer and C indices 
        for ci in C_indices_in_polymer:
            monomer_index = self.G.nodes[ci]['mi'] 
            if (not branching_state) and (monomer_index in terminal_monomers_indices):
                C_available_indices_in_polymer.append(ci)
            if (branching_state) and (monomer_index not in terminal_monomers_indices):
                C_available_indices_in_polymer.append(ci)
        
        return C_available_indices_in_polymer


    def find_available_C2_in_monomer(self, C1_index_in_polymer: int, ring: Optional[bool] = False) -> list:
        """find available C2 node in monomer 
        using the linkage selection rules

        Parameters
        ----------
        C1_index_in_polymer : int
            C1 indices in polymer

        Returns
        -------
        C2_indices_in_monmer : list
            C2 index in monomer
        """        
        # select the monomer type first and then C1 node index
        C1_node = self.G.nodes[C1_index_in_polymer]
        C2_indices_in_monmer = monomer_select_C1_C2[C1_node['mtype']][C1_node['index']]

        if ring:
            C2_indices_in_monmer = [ci for ci in C2_indices_in_monmer if ci != 1]
        #print(C2_indices_in_monmer)

        return C2_indices_in_monmer


    def find_available_monomer_types(self, bond_index: tuple) -> list:
        """find available monomer types 
        using monomer selection rules

        Parameters
        ----------
        bond_index : tuple
            bonding C index pair, e.g., (4, 5)

        Returns
        -------
        monomer type list : list
            a list monomer available for forming linkages
        """        
        # Use linkage to select
        monomer_types = linkage_index_select_monomer[bond_index]
        
        return monomer_types


    def find_C2_index_in_polymer(self, C2_index_in_monomer: int) -> int:
        """Convert C2 index in monomer to its index in the polymer

        Parameters
        ----------
        C2_index_in_monomer : int
            C2 index in the monomer

        Returns
        -------
        C2_index_in_polymer : int
            C2 index in the polymer
        """           
        # -1 adjust for python indexing
        C2_index_in_polymer = C2_index_in_monomer + len(self.G) - 1 

        return C2_index_in_polymer


    def find_O_index_in_polymer(self, C_index_in_polymer: int) -> int:
        """find the oxygen neighboring to the 4th carbon in a graph

        Parameters
        ----------
        C_index_in_polymer : int
            the index of 4th carbon 

        Returns
        -------
        O_index_in_polymer : int
            the oxygen connected to the 4th carbon 
        """           
        O_index_in_polymer = [i for i in list(self.G.neighbors(C_index_in_polymer)) if (self.G.nodes[i]['element'] == 'O')][0]

        return O_index_in_polymer


    def connect_C1_C2(
        self, 
        linkage_index: Tuple[int, int], 
        C1_index_in_polymer: int, 
        C2_index_in_polymer: int
    ) -> bool:
        """Add a new linkage, connect C1 and C2 in the polymer

        Parameters
        ----------
        linkage_index : Tuple[int, int]
            (C1, C2) index tuple for the linkage
        C1_index_in_polymer : int
            C1 index in the polymer
        C2_index_in_polymer : int
            C2 index in the polymer

        Returns
        -------
        new_linkage_flag : bool
            a flag if new linkage is formed
        """    
        # a flag if new linkage is formed
        new_linkage_flag = False
        
        # Define the edge properties
        #new bond as carbon/oxygen indices
        linkage_new = linkage_index 
        print(linkage_new)
        C1_node = self.G.nodes[C1_index_in_polymer]
        C2_node = self.G.nodes[C2_index_in_polymer]
        
        # linkage monomer types
        linkage_new_types = (C1_node['mtype'], C2_node['mtype'])
        # linkage name
        linkage_new_name = linkage_index_to_name[linkage_new]
        
        # Add the edges
        linkage_new_index_list = []
        
        # additional atoms in the ring formation that need to be make unavailable
        O_index_in_polymer = None
        alpha_index_in_polymer = None
        alpha1_index_in_polymer = None
        alpha2_index_in_polymer = None
        C4_index_in_polymer = None
        
        # Check if the bond requires special bonding rules 
        if linkage_new_name in linkage_special_names: 
            
            # carbon 4 in the bond, C1-O and C2-O needs to be connected, alpha-O-4 and beta-O-4
            if 4 in linkage_new:
                
                if linkage_new[0] == 4: #C1 is the 4th carbon
                    O_index_in_polymer = self.find_O_index_in_polymer(C1_index_in_polymer)
                
                if linkage_new[1] == 4: #C2 is the 4th carbon
                    O_index_in_polymer = self.find_O_index_in_polymer(C2_index_in_polymer)
            
                linkage_new_index_C1_O = (C1_index_in_polymer, O_index_in_polymer)
                linkage_new_index_O_C2 = (O_index_in_polymer, C2_index_in_polymer)
                linkage_new_index_list = [linkage_new_index_C1_O, linkage_new_index_O_C2]
            
            # beta-5 linkage, connect beta(8th)-5 and O-alpha
            if linkage_new_name == 'beta-5':
                
                if linkage_new[0] == 8: #C1 is the beta 8th carbon, C2 is the 5th carbon
                    O_index_in_polymer = self.find_O_index_in_polymer(C2_index_in_polymer-1) #-1 to adjust to the 4th carbon
                    alpha_index_in_polymer = C1_index_in_polymer - 1
                    C4_index_in_polymer = C2_index_in_polymer - 1
                    
                if linkage_new[1] == 8: #C1 is the 5th carbon, C2 is the beta 8th carbon
                    O_index_in_polymer = self.find_O_index_in_polymer(C1_index_in_polymer-1) #-1 to adjust to the 4th carbon
                    alpha_index_in_polymer = C2_index_in_polymer - 1
                    C4_index_in_polymer = C1_index_in_polymer - 1
                
                # check bonding availablity for additional carbon in the ring
                if self.G.nodes[C4_index_in_polymer]['bonding'] and self.G.nodes[alpha_index_in_polymer]['bonding']:
                    linkage_new_index_C1_C2 = (C1_index_in_polymer, C2_index_in_polymer)
                    linkage_new_index_O_alpha = (O_index_in_polymer, alpha_index_in_polymer)
                    linkage_new_index_list = [linkage_new_index_C1_C2, linkage_new_index_O_alpha]
            
            # beta-beta linkage, 8th C to 8th, need to connect beta - beta, alpha in 1 to O and alpha in 2 to O and 
            if linkage_new_name == 'beta-beta':
                
                O1_index_in_polymer = self.find_O_index_in_polymer(C1_index_in_polymer+1) #+1 to adjust to the 9th carbon
                O2_index_in_polymer = self.find_O_index_in_polymer(C2_index_in_polymer+1) #+1 to adjust to the 9th carbon
                
                alpha1_index_in_polymer = C1_index_in_polymer - 1 #-1 to adjust to the alpha(7th) carbon
                alpha2_index_in_polymer = C2_index_in_polymer - 1 #-1 to adjust to the alpha(7th) carbon
                
                # check bonding availablity for additional carbon in the ring
                if self.G.nodes[alpha1_index_in_polymer]['bonding'] and self.G.nodes[alpha2_index_in_polymer]['bonding']:
                    linkage_new_index_C1_C2 = (C1_index_in_polymer, C2_index_in_polymer)
                    linkage_new_index_alpha1_O2 = (alpha1_index_in_polymer, O2_index_in_polymer)
                    linkage_new_index_alpha2_O1 = (alpha2_index_in_polymer, O1_index_in_polymer)
                    
                    linkage_new_index_list = [linkage_new_index_C1_C2, linkage_new_index_alpha1_O2, linkage_new_index_alpha2_O1]
            
            if linkage_new_name == 'beta-1':
                
                O_index_in_polymer = len(self.G) # index to -OH that needs to be added at alpha(7th) of the first monomer

                if linkage_new[0] == 8: #C1 is the beta 8th carbon, C2 is the 1th carbon
                    alpha1_index_in_polymer = C1_index_in_polymer-1 # the bonding alpha carbon
                    alpha2_index_in_polymer = C2_index_in_polymer-1 # the atoms connect to this needs to be deleted
                else: #C2 is the beta 8th carbon, C1 is the 1th carbon
                    alpha1_index_in_polymer = C2_index_in_polymer-1
                    alpha2_index_in_polymer = C1_index_in_polymer-1

                # Extract the bonding alpha node 
                alpha1_node = self.G.nodes[alpha1_index_in_polymer]

                # remove atoms that connects to C2
                atom_indices_to_delete = [7, 8, 9, 11]
                atom_indices_to_delete_in_polymer = [ni + alpha2_index_in_polymer for ni in atom_indices_to_delete]
                self.G.remove_nodes_from(atom_indices_to_delete_in_polymer)

                # Add the new -OH group 
                O_3_index = 16 # assign a new index 
                # Use the properties same as the bonding alpha carbon
                self.G.add_node(O_index_in_polymer, element = 'O', aromatic = False, group = '7OH', index = O_3_index, \
                    mtype = alpha1_node['mtype'], bonding = False, color = alpha1_node['color'], mi=alpha1_node['mi'])
                
                # find the old alpha beta bond
                # beta of the 1nd monomer is now the alpha of the second
                alpha2_index_in_polymer = alpha1_index_in_polymer+1 # set back to none
                linkage_old_index_alpha_beta = (alpha1_index_in_polymer, alpha2_index_in_polymer)
                # set the alpha-beta bond order back to 1
                self.G.edges[linkage_old_index_alpha_beta]['order'] = 1

                linkage_new_index_C1_C2 = (C1_index_in_polymer, C2_index_in_polymer)
                linkage_new_index_alpha1_O = (alpha1_index_in_polymer, O_index_in_polymer)
                linkage_new_index_list = [linkage_new_index_C1_C2, linkage_new_index_alpha1_O]


        # Normal bonding, connect two nodes
        else:
            linkage_new_index = (C1_index_in_polymer, C2_index_in_polymer)
            linkage_new_index_list = [linkage_new_index]
            
        # add the bond list to graph
        if len(linkage_new_index_list) > 0:
            self.G.add_edges_from(linkage_new_index_list, order= 1, index = linkage_new, mtype = linkage_new_types, btype = linkage_new_name)
           
            # change the bonding availablity
            self.G = ut.make_unavailable(self.G, C1_index_in_polymer)
            self.G = ut.make_unavailable(self.G, C2_index_in_polymer)
            
            # print out the progress
            if self.verbose: 
                print('Connect a {} to a {} unit (in polymer) via a {} bond'.format(linkage_new_types[0], linkage_new_types[1], linkage_new_name))
                
            # set the new linkage flag
            new_linkage_flag = True
                
        # make other carbon bonding with oxygen if necessary unavailable
        if not alpha_index_in_polymer == None: self.G = ut.make_unavailable(self.G, alpha_index_in_polymer)
        if not alpha1_index_in_polymer == None: self.G = ut.make_unavailable(self.G, alpha1_index_in_polymer)
        if not alpha2_index_in_polymer == None: self.G = ut.make_unavailable(self.G, alpha2_index_in_polymer)
        if not C4_index_in_polymer == None: self.G = ut.make_unavailable(self.G, C4_index_in_polymer)
        
        # For beta-1, some node indices need to be adjusted due to the node removal
        if linkage_new_name == 'beta-1' and new_linkage_flag:
            self.G = ut.adjust_indices(self.G)

        # return new_linkage_flag
        return new_linkage_flag    



class Polymer(PolymerGraph):
    """
    Polymer object for lignin growth
    """    
    def __init__(self, M_init: object, verbose: Optional[bool] = True):
        """initialize the object with a graph

        Parameters
        ----------
        M_init : Monomer or Polymer object
            object of the initial monomer/polymer
        verbose : Optional[bool], optional
            the flag to control on-off for print statement, by default True
        """        
        # make a copy of the initial object
        #M_init_copy = copy(M_init)
        # Create a graph if the monomer has none
        if M_init.G is None and isinstance(M_init, Monomer):
            M_init.create()
        
        # Initialize the graph and big graph
        super().__init__(M_init.G, verbose)
        self.bigG = M_init.bigG.copy()

        
        # Update the initial monomerindex
        # check if the input is a monomer object, initialize a new one
        if isinstance(M_init, Monomer):
            self.mi = 0 #index of the current monomer
        # else input is a polymer object
        else:
            self.mi = M_init.mi

        
    def update_units(
        self,
        linkage_index: Tuple[int, int],
        C1_index_in_polymer: int, 
        C2_index_in_polymer: int
    ):
        """Update the units dictionary in polymer object

        Parameters
        ----------
        linkage_index : Tuple[int, int]
            (C1, C2) index tuple for the linkage
        C1_index_in_polymer : int
            C1 index in the polymer
        C2_index_in_polymer : int
            C2 index in the polymer
        """    
        # linkage name
        linkage_new_name = linkage_index_to_name[linkage_index]

        # Find the monomer index of C1 and C2
        M1_index = self.G.nodes[C1_index_in_polymer]['mi']
        M2_index = self.G.nodes[C2_index_in_polymer]['mi']

        # Add the edges to monomer nodes
        print(linkage_new_name, M1_index, M2_index)
        self.bigG.add_edges_from([(M1_index, M2_index)], btype=linkage_new_name)


    def try_adding_new(
        self,
        linkage_index: Tuple[int, int],
        C1_index_in_polymer: int,
        C2_index_in_polymer: int, 
        monomer_new: Optional[Monomer] = None,
        draw: Optional[bool] = False
    ) -> bool:
        """Try adding a new monomer or linkage

        Parameters
        ----------
        linkage_index : Tuple[int, int]
            (C1, C2) index tuple for the linkage
        C1_index_in_polymer : int
            C1 index in the polymer
        C2_index_in_polymer : int
            C2 index in the polymer
        monomer_new : Optional[Monomer], optional
            monomer object for the new monomer, 
            by default None, i.e. forming a ring
        draw : Optional[bool], optional
            flag to draw the polymer graph, by default False

        Returns
        -------
        new_linkage_flag : bool
            a flag if new linkage is formed
        """        

        # Add the monomer into the graph
        # Here we create a temp PolymerGraph object to connect C1 and C2
        PG_temp = PolymerGraph(self.G.copy(), verbose=self.verbose)

        # join the monomer graph and polymer graph
        if monomer_new is not None:
            monomer_new.create()
            PG_temp.G = ut.join_two(PG_temp.G, monomer_new.G)
        print("Linkage index", linkage_index)
        # Add the linkage
        new_linkage_flag = PG_temp.connect_C1_C2(linkage_index, C1_index_in_polymer, C2_index_in_polymer)
        print("New Linkage Flag")
        # If a new linakge does form
        if new_linkage_flag: 
            # Update the polymer graph
            self.G = PG_temp.G
            # Update the polymer big graph
            if monomer_new is not None:
                if linkage_index[0] < linkage_index[1]:
                    self.bigG = ut.join_two(self.bigG, monomer_new.bigG)
                else:
                    self.bigG = ut.join_two(monomer_new.bigG, self.bigG)
            # Update the avaible C1 list after addition of a new monomer and a new linkage  
            self.C1_indices_in_polymer = self.find_available_C1_in_polymer()
            # Update the unit dictionary
            self.update_units(linkage_index, C1_index_in_polymer, C2_index_in_polymer)
            # Update the current monomer index
            self.mi += 1

        if draw: ut.draw_graph(self.G)

        return new_linkage_flag


    def add_specific_linkage(
        self, 
        linkage_type: str, 
        monomer_type: Optional[str] = None,
        branching_state: Optional[bool] = None,
        draw: Optional[bool] = False):
        """Add specific linkage to a polymer with a specific monomer if user defined 
        Otherwise add an extra linkage allowed by the linkage rules

        Parameters
        ----------
        linkage_type : str
            bond types defined in bond type selection
        monomer_type : Optional[str], optional
            monomer type, by default None
        draw : Optional[bool], optional
            flag to draw the polymer graph, by default False
        branching_state : Optional[bool], optional
            likelihood to form a branch
            by default None, select all Cs
            if true, select Cs in the chain
            if false, select Cs at the terminal monomers

        Raises
        ------
        Exception
            if input linkage type not defined
        Exception
            if input linkage does not exist
        """
        # Find available bonding nodes
        # select the linkage for the second C
        new_linkage_flag = False
        if linkage_type in linkage_name_select_C1_C2.keys():
            linkages_possible = linkage_name_select_C1_C2[linkage_type]
        else:    
            raise Exception("Input linkage type is not defined")
        
        # Select available C1 node
        if (self.C1_indices_in_polymer is None) or (branching_state is not None):
            self.C1_indices_in_polymer = self.find_available_C1_in_polymer(branching_state)
            
        # Select available C1 node and defined by the bonding
        C1_indices_in_monomer_possible = linkage_name_select_C1_C2[linkage_type].keys()
        
        # Find the intersections
        C1_indices_in_polymer = []
        for ci in self.C1_indices_in_polymer:
            if self.G.nodes[ci]['index'] in C1_indices_in_monomer_possible:
                C1_indices_in_polymer.append(ci)
            else: pass
        
        # Case 1, if there is no C1 available for bonding
        if len(C1_indices_in_polymer) == 0:
            if self.verbose: 
                warnings.warn("No more bonding atoms available")
            
        else: 
            # Case 2, select the random carbon 1 out of the suitable one
            C1_index_in_polymer = ut.select_one_from_many(C1_indices_in_polymer)
            C1_node = self.G.nodes[C1_index_in_polymer]
            C1_index_in_monomer = C1_node['index']
    
            if C1_index_in_monomer in linkages_possible.keys():
                 C2_index_in_monomer = linkages_possible[C1_index_in_monomer]
            else:    
                raise Exception("Input linkage type does not exist")
            
            linkage_index = (C1_index_in_monomer, C2_index_in_monomer)
            

            # Find the available monomer 
            monomer_types_possible = self.find_available_monomer_types(linkage_index)
            
            # Case 3, select a random monomer
            if not monomer_type == None:
                if not monomer_type in monomer_types_possible:
                    raise Exception("Input monomer type is not allowed")
                else: 
                    pass
            else: # if no user input, select randomly
                monomer_type = ut.select_one_from_many(monomer_types_possible)
            
            # Create a new monomer object at the given index
            monomer_new = Monomer(monomer_type, monomer_index=self.mi+1)
            
            # Find possible C2 index
            C2_index_in_polymer = self.find_C2_index_in_polymer(C2_index_in_monomer)


            # Try adding the monomer, linkage into the graph
            new_linkage_flag = self.try_adding_new(linkage_index, C1_index_in_polymer, C2_index_in_polymer, \
                monomer_new=monomer_new, draw=draw)

        return new_linkage_flag


    def add_specific_monomer(
        self, 
        monomer_type: str, 
        linkage_type: Optional[str] = None, 
        branching_state: Optional[bool] = None,
        draw: Optional[bool] = False):
        """Add specific linkage to a polymer with a specific monomer if user defined 
        Otherwise add an extra linkage allowed by the linkage rules

        Parameters
        ----------
        monomer_type : str
            monomer type
        linkage_type : Optional[str], optional
            linkage types defined in bond type selection, by default None
        branching_state: Optioanl[bool], optional
            likelihood to form a branch
            by default None, select all Cs
            if true, select Cs in the chain
            if false, select Cs at the terminal monomers
        draw : Optional[bool], optional
            flag to draw the polymer graph, by default False
        """        
        new_linkage_flag = False
        # create a new monomer based on its type
        monomer_new = Monomer(monomer_type, monomer_index=self.mi+1)
        Gm = monomer_new.create()
        
        # Select available C1 node in polymer
        if (self.C1_indices_in_polymer is None) or (branching_state is not None):
            self.C1_indices_in_polymer = self.find_available_C1_in_polymer(branching_state)

        C1_indices_in_polymer = self.C1_indices_in_polymer
        
        # Case 1, if there is no C1 available for bonding
        if len(C1_indices_in_polymer) == 0:
            warnings.warn("No more bonding atoms available")
            
        else: 

            # Case 2, if there is no available linkage defined the user
            # find the corresponding C1 indices in the monomers of the polymer 
            C1_indices_in_monomer_of_polymer = [self.G.nodes[ci]['index'] for ci in C1_indices_in_polymer]
            # find the C2 indices in the new monomer
            C2_indices_in_monmer_node = [n for n,v in Gm.nodes(data=True) if v['bonding']]
            C2_indices_in_monomer = [Gm.nodes[n]['index'] for n in C2_indices_in_monmer_node]
            C2_indices_in_polymer = [self.find_C2_index_in_polymer(ci) for ci in C2_indices_in_monomer]
            # combine C1 and C2 as possible linkages
            linkage_indices_in_monomer_possible = []
            linkage_indices_in_polymer_possible = []
            for c1m_i, c1p_i in zip(C1_indices_in_monomer_of_polymer, C1_indices_in_polymer):
                for c2m_i, c2p_i in zip(C2_indices_in_monomer, C2_indices_in_polymer):
                    linkage_indices_in_monomer_possible.append((c1m_i, c2m_i))
                    linkage_indices_in_polymer_possible.append((c1p_i, c2p_i))
            
            # find the corresponding linkage type
            linkage_possible_indices =  [ind for ind, li in enumerate(linkage_indices_in_monomer_possible) if li in linkage_index_to_name.keys()]
            linkage_indices_in_monomer_possible = [linkage_indices_in_monomer_possible[ind] for ind in linkage_possible_indices]
            linkage_indices_in_polymer_possible = [linkage_indices_in_polymer_possible[ind] for ind in linkage_possible_indices]
            linkage_possible = [linkage_index_to_name[li] for li in linkage_indices_in_monomer_possible]
            
            
            # select the linkage_type by user input
            linkage_indices_in_monomer_possible_input = []
            linkage_indices_in_polymer_possible_input = []

            if not linkage_type ==  None:
                for type_i, m_i, p_i in zip(linkage_possible, linkage_indices_in_monomer_possible, linkage_indices_in_polymer_possible):
                    if type_i == linkage_type:
                        linkage_indices_in_monomer_possible_input.append(m_i)
                        linkage_indices_in_polymer_possible_input.append(p_i)
            else:
                linkage_indices_in_monomer_possible_input = linkage_indices_in_monomer_possible
                linkage_indices_in_polymer_possible_input = linkage_indices_in_polymer_possible
            
            if len(linkage_indices_in_monomer_possible_input) == 0:    
                warnings.warn("Input linkage type is not supported")
                
            else:
                # Case 3, select a random linkage
                # select one linkage out of the possible ones 
                linkage_index, mi = ut.select_one_from_many_with_index(linkage_indices_in_monomer_possible_input)
                linkage_index_polymer = linkage_indices_in_polymer_possible_input[mi]
    
                C1_index_in_polymer = linkage_index_polymer[0]
                C2_index_in_polymer = linkage_index_polymer[1]
    
                
                # Try adding the monomer, linkage into the graph
                new_linkage_flag = self.try_adding_new(linkage_index, C1_index_in_polymer, C2_index_in_polymer, \
                    monomer_new=monomer_new, draw=draw)
        
        return new_linkage_flag


    def add_random_monomer(
        self, 
        draw: Optional[bool] = False,
        branching_state: Optional[bool] = None):
        """Main body - Add one random monomer and connect the new linkage
        Rejection free

        Parameters
        ----------
        branching_state: Optioanl[bool], optional
            likelihood to form a branch
            by default None, select all Cs
            if true, select Cs in the chain
            if false, select Cs at the terminal monomers
        draw : Optional[bool], optional
            flag to draw the polymer graph, by default False
        """        
        new_linkage_flag = False
        # Find available bonding nodes
        # Select available C1 node
        if (self.C1_indices_in_polymer is None) or (branching_state is not None):
            self.C1_indices_in_polymer = self.find_available_C1_in_polymer(branching_state)
            
        C1_index_in_polymer = ut.select_one_from_many(self.C1_indices_in_polymer)

        C1_node = self.G.nodes[C1_index_in_polymer]
        C1_index_in_monomer = C1_node['index']

        # select the bond for the second C
        C2_indices_in_monomer = self.find_available_C2_in_monomer(C1_index_in_polymer)
        C2_index_in_monomer = ut.select_one_from_many(C2_indices_in_monomer)
        
        linkage_index = (C1_index_in_monomer, C2_index_in_monomer)

        # Find the available monomer 
        monomer_types = self.find_available_monomer_types(linkage_index)
        monomer_type = ut.select_one_from_many(monomer_types)
        monomer_new = Monomer(monomer_type, monomer_index=self.mi+1)
        
        # Add the monomer into the graph
        C2_index_in_polymer = self.find_C2_index_in_polymer(C2_index_in_monomer)
        
        # Try adding the monomer, linkage into the graph
        new_linkage_flag = self.try_adding_new(linkage_index, C1_index_in_polymer, C2_index_in_polymer, \
            monomer_new=monomer_new, draw=draw)

        return new_linkage_flag
        
        
    def add_random_ring(
        self, 
        draw: Optional[bool] = False,
        branching_state: Optional[bool] = None):
        """Add a randome linkage inside a polymer
        a ring is formed

        Parameters
        ----------
        branching_state: Optioanl[bool], optional
            likelihood to form a branch
            by default None, select all Cs
            if true, select Cs in the chain
            if false, select Cs at the terminal monomers
        draw : Optional[bool], optional
            flag to draw the polymer graph, by default False
        """        
        new_linkage_flag = False
        # Find available bonding nodes
        # Select available C1 node
        if (self.C1_indices_in_polymer is None) or (branching_state is not None):
            self.C1_indices_in_polymer = self.find_available_C1_in_polymer(branching_state)
        
        # Case 1, if there is no C1 available for bonding
        C1_indices_in_polymer = [ci for ci in self.C1_indices_in_polymer if self.G.nodes[ci]['index'] != 1]

        if len(C1_indices_in_polymer) == 0:
            pass #raise Exception("No more bonding atoms available")
        
        else: 
            # Case 2, select the random carbon 1 out of the suitable one
            C1_index_in_polymer = ut.select_one_from_many(C1_indices_in_polymer)
            C1_node = self.G.nodes[C1_index_in_polymer]
            C1_index_in_monomer = C1_node['index']

            # select the second C from other Carbon avaible for linkage inside the polymer
            C2_indices_in_polymer = self.update_available_C1_in_polymer(C1_index_in_polymer)
            C2_indices_in_monomer_of_polymer = [self.G.nodes[ci]['index'] for ci in C2_indices_in_polymer]
            
            # base on C1, select the possible C2 index in monomer and polymer
            C2_indices_in_monomer_possible = self.find_available_C2_in_monomer(C1_index_in_polymer, ring = True)
            C2_index_in_monomer_of_polymer_possible = []
            C2_index_in_polymer_possible = []
            
            for mi, pi in zip(C2_indices_in_monomer_of_polymer, C2_indices_in_polymer):
                if mi in C2_indices_in_monomer_possible:
                    C2_index_in_monomer_of_polymer_possible.append(mi)
                    C2_index_in_polymer_possible.append(pi)
                    
            # Case 3, if there is no more available C2, pass
            if len(C2_index_in_monomer_of_polymer_possible) == 0:
                pass
            
            else:
                # Continue Case 2, add the C1-C2 linkage
                C2_index_in_monomer, pi = ut.select_one_from_many_with_index(C2_index_in_monomer_of_polymer_possible)
                linkage_index = (C1_index_in_monomer, C2_index_in_monomer)
    
                # Add the monomer into the graph
                C2_index_in_polymer = C2_index_in_polymer_possible[pi]
                

                # Try adding the linkage into the graph
                new_linkage_flag = self.try_adding_new(linkage_index, C1_index_in_polymer, C2_index_in_polymer, \
                    monomer_new=None, draw=draw)
        return new_linkage_flag


    def add_specific_ring(
        self, 
        linkage_type: str, 
        draw: Optional[bool] = False):
        """Add specific linkage to a polymer with a specific monomer if user defined 
        Otherwise add an extra linkage allowed by the linkage rules

        Parameters
        ----------
        linkage_type : str
            bond types defined in bond type selection
        monomer_type : Optional[str], optional
            monomer type, by default None
        draw : Optional[bool], optional
            flag to draw the polymer graph, by default False

        Raises
        ------
        Exception
            if input linkage type not defined
        Exception
            if input linkage does not exist
        Exception
            if input monomer is not allowed by the bonding
        """
        new_linkage_flag = False
        # Find available bonding nodes
        # select the linkage for the second C
        if linkage_type in linkage_ring:
            linkages_possible = linkage_name_select_C1_C2[linkage_type]
        else:    
            return new_linkage_flag #raise Exception("Input linkage type is not defined")
        
        # Select available C1 node
        if self.C1_indices_in_polymer == None:
            self.C1_indices_in_polymer = self.find_available_C1_in_polymer()
            
        # Select available C1 node and defined by the bonding
        C1_indices_in_monomer_possible = linkage_name_select_C1_C2[linkage_type].keys()
        
        # Find the intersections
        C1_indices_in_polymer = []
        for ci in self.C1_indices_in_polymer:
            if self.G.nodes[ci]['index'] in C1_indices_in_monomer_possible:
                C1_indices_in_polymer.append(ci)
            else: 
                #print("No C1 available for bonding")
                pass
        
        # Case 1, if there is no C1 available for bonding
        if len(C1_indices_in_polymer) == 0:
            if self.verbose: 
                warnings.warn("No more bonding atoms available")


        else: 
            # Case 2, select the random carbon 1 out of the suitable one
            C1_index_in_polymer = ut.select_one_from_many(C1_indices_in_polymer)
            C1_node = self.G.nodes[C1_index_in_polymer]
            C1_index_in_monomer = C1_node['index']
            
            
            # select the second C from other Carbon avaible for linkage inside the polymer
            C2_indices_in_polymer = self.update_available_C1_in_polymer(C1_index_in_polymer)
            C2_indices_in_monomer_of_polymer = [self.G.nodes[ci]['index'] for ci in C2_indices_in_polymer]

            C2_indices_in_monomer_possible = linkages_possible[C1_index_in_monomer] #self.find_available_C2_in_monomer(C1_index_in_polymer, ring = True)
            C2_index_in_monomer_of_polymer_possible = []
            C2_index_in_polymer_possible = []
            
            for mi, pi in zip(C2_indices_in_monomer_of_polymer, C2_indices_in_polymer):
                if mi == C2_indices_in_monomer_possible:
                    C2_index_in_monomer_of_polymer_possible.append(mi)
                    C2_index_in_polymer_possible.append(pi)
                    
            # Case 3, if there is no more available C2, pass
            if len(C2_index_in_monomer_of_polymer_possible) == 0:
                if self.verbose: 
                    warnings.warn("No C2 available for bonding")


            else:
                # Continue Case 2, add the C1-C2 linkage
                C2_index_in_monomer, pi = ut.select_one_from_many_with_index(C2_index_in_monomer_of_polymer_possible)
                linkage_index = (C1_index_in_monomer, C2_index_in_monomer)
                #print(linkage_index)
                # Add the monomer into the graph
                C2_index_in_polymer = C2_index_in_polymer_possible[pi]

                
                # Try adding the monomer, linkage into the graph
                new_linkage_flag = self.try_adding_new(linkage_index, C1_index_in_polymer, C2_index_in_polymer, \
                    monomer_new=None, draw=draw)
        
        return new_linkage_flag



               
