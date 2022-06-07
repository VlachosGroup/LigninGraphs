===================
Hyperparameters
===================


Hyperparameters are used to facilitate structure optimization. 
The values do not affect the outcome but control the runtime and accelerate convergence to the target values.
Here is a complete list of hyperparameters. The usage can be found in example 10, 11, 12, 13.
Their definition and suggested value are listed below:

Metropolis temperature
----------------------------
:code:`Tmetro`, :code:`Tmetro_out`: Greater than 0, usually smaller than 300 K. 

Control the acceptance rate for inner and outer optimization loop respectively.
The higher temperature, more likely a Monte Carlo move will be accepted.


Max number of iterations in each loop 
----------------------------------------
:code:`i_max`, :code:`i_max_out`, :code:`i_max_ring`: Usually ~1000 -10000

Limit the number of MC attempts and runtime of the simulations in each optimization loop.


Expected polymer size 
--------------------------
:code:`expected_size`: Greater than 0, usually smaller than 100. 

We consider the polymer size from a population following a normal distribution 
This value sets the estimated mean size for the polymer molecule in the inner loop. 

Maximum polymer size 
-----------------------------
:code:`max_size`: Greater than 0, usually smaller than 100. 

Set a threshold for the polymer molecule size in the inner loop. 


Size distribution scaling factor 
-------------------------------------
:code:`distribution_scaling`: Between 0 to 1.

Controls the width of this distribution and affects the weighted average molecular weight. 
The larger the k, the wider the range of polymer sizes is; the weighted average molecular weight would be 
greater than the number average molecular weight, indicating the polymer population is polydisperse. 

Branching propensity 
------------------------
:code:`branching_propensity`: Between 0 to 1.

Determine the branching probability for newly added monomers. 
The value should be similar to the experimental branching coefficient.

Hyperparameters can be passed into :code:`optimization.Simulation` class.
Code example to define all the hyperparameters: 

.. code-block::

    # set Metropolis temperatures
    Tmetro = 10
    Tmetro_out = 10

    # set the max iterations
    i_max = 1000
    i_max_out = 1000
    i_max_ring = 500

    # set the expected and max sizes
    expected_size = 7
    max_size = 100
    # size distribution scaling factor
    distribution_scaling = 0.1

    # whether the size is in Molecular weight or monomer count
    size_in_MW = False

    branching_propensity = 0
    
    
    



