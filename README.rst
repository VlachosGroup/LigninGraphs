==============
LigninGraphs
==============
 

LigninGraphs is an open-source software package in Python to generate feasible lignin structures. 

Lignin is an aromatic biopolymer found in ubiquitous sources of woody biomass such as wood and bark.
Designing and optimizing lignin valorization processes requires a fundamental understanding of lignin structures.
We introduce a graph-based multiscale modeling framework for lignin structure generation and visualization. 
The framework employs accelerated rejection-free polymerization and hierarchical Metropolis Monte Carlo optimization algorithms. 
It can be used to generate feasible lignin strutcures to match experimental or literature data. 

.. image:: docs/source/logos/ligning_logo.png
    :width: 400px

Documentation
-------------

See our `documentation page`_ for examples, equations used, and docstrings. Examples notebooks can be opened directly in Google Colab.

.. image:: https://colab.research.google.com/assets/colab-badge.svg
   :alt: Open In Colab
   :target: http://colab.research.google.com/github/VlachosGroup/LigninGraphs/blob/master/

Developers
----------

-  Yifan Wang (wangyf@udel.edu)
-  Jake Kalscheur (jkalsch@udel.edu)

Dependencies
------------

-  Python >= 3.7
-  `RDKit`_ >= 2021.09.1: Used for constructing feasible chemical structures 
-  `Networkx`_ >= 1.4: Used for computational graph operations
-  `Pysmiles`_ >= 1.0.1: Used for reading and writing smiles strings from/to graphs
-  `Matplotlib`_: Used for generating plots
-  `Numpy`_: Used for vector and matrix operations
-  `Scipy`_: Used for curve fitting
-  `Pandas`_: Used to import data from Excel or CSV files
-  `openpyxl`_: Used by Pandas to import Excel files
-  `pytest`_: Used for unit tests


.. _documentation page: https://ligningraphs.readthedocs.io/en/latest/
.. _RDKit: https://www.rdkit.org/docs/Overview.html
.. _Networkx: https://networkx.org/
.. _Pysmiles: https://github.com/pckroon/pysmiles
.. _Matplotlib: https://matplotlib.org/
.. _Numpy: http://www.numpy.org/
.. _Scipy: https://www.scipy.org/
.. _Pandas: https://pandas.pydata.org/
.. _openpyxl: https://openpyxl.readthedocs.io/en/stable/
.. _pytest: https://docs.pytest.org/en/stable/



Getting Started
---------------

1. Install using pip (see documentation for full instructions)::

    pip install ligning

2. Run the unit tests.

3. Read the documentation for tutorials and examples.


License
-------

This project is licensed under the MIT License - see the LICENSE.md.
file for details.


Contributing
------------

If you have a suggestion or find a bug, please post to our `Issues` page on GitHub. 

Questions
---------

If you are having issues, please post to our `Issues` page on GitHub.

Funding
-------

This material is based upon work supported by the Department of Energy's Office 
of Energy Efficient and Renewable Energy's Advanced Manufacturing Office under 
Award Number DE-EE0007888-9.5.

Acknowledgements
------------------

-  Siyi Huang (Logo design)
  

Publications
--------------

