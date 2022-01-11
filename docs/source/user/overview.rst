===================
Overview
===================

Our framework models lignin structures across scales, from the atomistic to the monomer to the polymer to a population of chains 
(i.e., a library of structures). We map structures to undirected graphs consisting of a list of nodes and connections (edges) between nodes, 
as in traditional graph theory. 
Graphs offer efficient data storage, fast lookup operations, and intuitive visualization for each node and its neighbors as an abstract data structure. 


At the atomic scale, each C or O atom is a node, and each bond is an edge (a,b). Atom properties are conveniently stored in a node. 
These include the element (“element”: C or O), the presence of an aromatic ring (“aromatic”: True or False), the monomer type (“mtype”: H, G, or S), 
the atomic index (“index”: the numeric labels for each C/O atom in a monomer), bonding (“bonding”: True or False), etc. 
Similarly, bond properties are stored in an edge. These include the linkage type (“btype”: the common ones, α-O-4, β-O-4, 5-5, 4-O-5, β-β, and β-5) or None for intra-monomer bonds), 
the bond order (“order”: 1 for a single or 2 for a double bond), atomic indices (“index”: atomic indices of the pair), and the monomer types of the bonding atoms (“mtype” of the pair). 
The resulting “atomic graphs” represent structures at the monomer (e) and the polymer (f) scale. 


We introduce a multiscale representation of the structure by coarse-graining upon molecular structures. 
Specifically, we coarse-grain the atomic graphs into “big graphs” by aggregating all nodes in a monomer into a single node (h), 
so that each monomer is a node and each linkage is an edge. 
Big graphs are connected at the polymer scale to make up a chain. 
They allow a significantly lower computational cost and storage by counting only the number of monomers or linkages. 
Finally, at the population scale of the polymer, we compute the observables or structure metrics of the entire structure library, 
including the monomer and linkage distributions, the number average molecular weight (MW), and the branching coefficient. 
These observables can be compared to experimental values. 