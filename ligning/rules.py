"""
Define the units and bonds in liginin
"""

# Notations
#       bond - bonding between atoms
#       linkage - a speical type of bond conneting two monomers, 
#       i.e., intermonomer bonds
#       C1 - the first bonding carbon atom in a linkage
#       C2 - the second bonding carbon atom in a linkage
#       (C1, C2) - the tuple for a linkage by specifying the indices of bonding carbons

# Carbon label to index mapping
#       alpha - 7
#       beta - 8
#       gamma - 9

# Linkage mapping
# C-O linkages:
#       4-5, 5-4: 4-O-5
#       4-7, 7-4: alpha-O-4
#       4-8, 8-4: beta-O-4 
# C-C linkages:
#       5-5 :     5-5
#       5-8, 8-5: beta-5
#       8-8:      beta-beta1
#       1-8, 8-1: beta-1

# if 4 carbon is involved, the bonding includes the O
# see find_O_index_in_polymer

# the atom involved
CHO = ['C', 'H', 'O']

# atomic weights
weight_CHO = {'C': 12.011, 'H': 1.008, 'O': 15.999}

# The monomer type list 
monomer_types = ['H', 'G', 'S']

# default color list for monomers
default_color = {'H': 'lightcoral', 'G': 'lightgreen', 'S': 'lightblue'}

# linkage that involves ring formation 
linkage_ring = ['4-O-5', 'alpha-O-4', 'beta-O-4',  '5-5', 'beta-5', 'beta-beta']

# Include beta-1?
include_beta_1 = True

#%% With beta-1 rules
if include_beta_1:
      # The linkage name list
      linkage_names = ['4-O-5', 'alpha-O-4', 'beta-O-4', '5-5', 'beta-5', 'beta-beta', 'beta-1']

      # The linkage selection dictionary 
      # contains the monomer kind and corresponding
      # bonding carbons (1st and the 2nd)
      monomer_select_C1_C2 = {'H': {1: [8],
                                    4 : [5,7,8], 
                                    5 : [4,5,8], 
                                    7 : [4], 
                                    8 : [1,4,5,8]}, 
                              'G': {1: [8],
                                    4 : [5,7,8], 
                                    5 : [4,5,8], 
                                    7 : [4], 
                                    8 : [1,4,5,8]},
                              'S': {1: [8],
                                    4 : [5,7,8],              
                                    7 : [4], 
                                    8 : [1,4,5,8]}}


      # The monomer selection dictionary
      # contains the linkage tuple and the monomer kind which is allowed for the second carbon
      linkage_index_select_monomer = {(4, 5): ['H', 'G'],
                                    (4, 7): ['H', 'G', 'S'],
                                    (4, 8): ['H', 'G', 'S'],
                                    (5, 4): ['H', 'G', 'S'],
                                    (5, 5): ['H', 'G'],
                                    (5, 8): ['H', 'G', 'S'],
                                    (7, 4): ['H', 'G', 'S'],
                                    (8, 4): ['H', 'G', 'S'],
                                    (8, 5): ['H', 'G'],
                                    (8, 8): ['H', 'G', 'S'],
                                    (1, 8): ['H', 'G', 'S'],
                                    (8, 1): ['H', 'G', 'S']}


      # The linkage name selection dictionary
      # contain the linkage name and the linkage pairs
      linkage_name_select_C1_C2 = {'4-O-5': {4: 5, 5: 4},
                                    'alpha-O-4': {4: 7, 7: 4},
                                    'beta-O-4':  {4: 8, 8: 4},
                                    '5-5': {5: 5},
                                    'beta-5': {8: 5, 5: 8},
                                    'beta-beta': {8: 8},
                                    'beta-1': {1: 8, 8: 1}}


      # linkage indices with corresponding linkage names
      linkage_index_to_name = {(4, 5): '4-O-5',
                              (4, 7): 'alpha-O-4',
                              (4, 8): 'beta-O-4',
                              (5, 4): '4-O-5',
                              (5, 5): '5-5',
                              (5, 8): 'beta-5',
                              (7, 4): 'alpha-O-4',
                              (8, 4): 'beta-O-4',
                              (8, 5): 'beta-5',
                              (8, 8): 'beta-beta',
                              (1, 8): 'beta-1',
                              (8, 1): 'beta-1'}

      # the list contain special linkages which need additional rules than connecting two nodes
      # Special linkages include 4-O-5, alpha-O-4, beta-O-4 and beta-5
      # 4-O-5, alpha-O-4 and beta-O-4 puts an oxygen in between the two Cs
      # beta-5 forms a ring
      # beta-beta forms two rings
      linkage_special_names = ['4-O-5', 'alpha-O-4', 'beta-O-4', 'beta-5', 'beta-beta', 'beta-1']


      # check list for different monomer linkage combination used in testing
      monomer1_select_monomer2_linkage_name = {'H': {'H': ['alpha-O-4', 'beta-O-4', '5-5', 'beta-5', 'beta-beta', 'beta-1'],
                                                'G': ['4-O-5','alpha-O-4', 'beta-O-4', '5-5', 'beta-5', 'beta-beta', 'beta-1'],
                                                'S': ['alpha-O-4', 'beta-O-4', 'beta-beta', 'beta-1']},
                                                'G': {'H': ['4-O-5','alpha-O-4', 'beta-O-4', '5-5', 'beta-5', 'beta-beta', 'beta-1'],
                                                      'G': ['4-O-5','alpha-O-4', 'beta-O-4', '5-5', 'beta-5', 'beta-beta', 'beta-1'],
                                                      'S': ['alpha-O-4', 'beta-O-4', 'beta-beta', 'beta-1']},
                                                'S': {'H': ['alpha-O-4', 'beta-O-4',  'beta-beta', 'beta-1'],
                                                      'G': ['alpha-O-4', 'beta-O-4', 'beta-beta', 'beta-1'],
                                                      'S': ['alpha-O-4', 'beta-O-4', 'beta-beta', 'beta-1']}}

#%% No beta-1 rules
else:       
      # The linkage name list
      linkage_names = ['4-O-5', 'alpha-O-4', 'beta-O-4', '5-5', 'beta-5', 'beta-beta']

      # The linkage selection dictionary 
      # contains the monomer kind and corresponding
      # bonding carbons (1st and the 2nd)
      monomer_select_C1_C2 = {'H': {4 : [5,7,8], 
                                    5 : [4,5,8], 
                                    7 : [4], 
                                    8 : [4,5,8]}, 
                              'G': {4 : [5,7,8], 
                                    5 : [4,5,8], 
                                    7 : [4], 
                                    8 : [4,5,8]},
                              'S': {4 : [5,7,8],              
                                    7 : [4], 
                                    8 : [4,5,8]}}

      # The monomer selection dictionary
      # contains the linkage tuple and the monomer kind which is allowed for the second carbon
      linkage_index_select_monomer = {(4, 5): ['H', 'G'],
                              (4, 7): ['H', 'G', 'S'],
                              (4, 8): ['H', 'G', 'S'],
                              (5, 4): ['H', 'G', 'S'],
                              (5, 5): ['H', 'G'],
                              (5, 8): ['H', 'G', 'S'],
                              (7, 4): ['H', 'G', 'S'],
                              (8, 4): ['H', 'G', 'S'],
                              (8, 5): ['H', 'G'],
                              (8, 8): ['H', 'G', 'S']}


      # The linkage name selection dictionary
      # contain the linkage name and the linkage pairs
      linkage_name_select_C1_C2 = {'4-O-5': {4: 5, 5: 4},
                              'alpha-O-4': {4: 7, 7: 4},
                              'beta-O-4':  {4: 8, 8: 4},
                              '5-5': {5: 5},
                              'beta-5': {8: 5, 5: 8},
                              'beta-beta': {8: 8}}


      # linkage indices with corresponding linkage names
      linkage_index_to_name = {(4, 5): '4-O-5',
                        (4, 7): 'alpha-O-4',
                        (4, 8): 'beta-O-4',
                        (5, 4): '4-O-5',
                        (5, 5): '5-5',
                        (5, 8): 'beta-5',
                        (7, 4): 'alpha-O-4',
                        (8, 4): 'beta-O-4',
                        (8, 5): 'beta-5',
                        (8, 8): 'beta-beta'}


      # the list contain special linkages which need additional rules than connecting two nodes
      # Special linkages include 4-O-5, alpha-O-4, beta-O-4 and beta-5
      linkage_special_names = ['4-O-5', 'alpha-O-4', 'beta-O-4', 'beta-5', 'beta-beta']
      

      # check list for different monomer linkage combination used in testing
      monomer1_select_monomer2_linkage_name = {'H': {'H': ['alpha-O-4', 'beta-O-4', '5-5', 'beta-5', 'beta-beta'],
                                          'G': ['4-O-5','alpha-O-4', 'beta-O-4', '5-5', 'beta-5', 'beta-beta'],
                                          'S': ['alpha-O-4', 'beta-O-4', 'beta-beta']},
                                    'G': {'H': ['4-O-5','alpha-O-4', 'beta-O-4', '5-5', 'beta-5', 'beta-beta'],
                                          'G': ['4-O-5','alpha-O-4', 'beta-O-4', '5-5', 'beta-5', 'beta-beta'],
                                          'S': ['alpha-O-4', 'beta-O-4', 'beta-beta']},
                                    'S': {'H': ['alpha-O-4', 'beta-O-4',  'beta-beta'],
                                          'G': ['alpha-O-4', 'beta-O-4', 'beta-beta'],
                                          'S': ['alpha-O-4', 'beta-O-4', 'beta-beta']}}



