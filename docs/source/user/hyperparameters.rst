===================
Hyperparameters
===================


Hyperparameters are used to facilitate structure optimization. 
The values do not affect the outcome but control the runtime and accelerate convergence to the target values.
Here is a complete list of hyperparameters. The usage can be found in example 10, 11, 12, 13.
Their definition and suggested value are listed below:

Metropolis temperature
----------------------------
`Tmetro`, `Tmetro_out`: 
Control the acceptance rate for inner and outer optimization loop respectively.
The higher temperature, more likely a Monte Carlo move will be accepted.
Greater than 0, usually smaller than 300 K. 


Max number of iterations in each loop 
----------------------------------------
`i_max`, `i_max_out`, `i_max_ring`:
Limit the number of MC attempts and runtime of the simulations in each optimization loop.
Usually ~1000 -10000


Expected polymer size 
--------------------------
`expected_size`:
We consider the polymer size from a population following a normal distribution 
This value sets the estimated mean size for the polymer molecule in the inner loop. 
Greater than 0, usually smaller than 100. 

Maximum polymer size 
-----------------------------
`max_size`:
Set a threshold for the polymer molecule size in the inner loop. 
Greater than 0, usually smaller than 100. 


Size distribution scaling factor 
-------------------------------------
`distribution_scaling`:
Controls the width of this distribution and affects the weighted average molecular weight. 
The larger the k, the wider the range of polymer sizes is; the weighted average molecular weight would be 
greater than the number average molecular weight, indicating the polymer population is polydisperse. 
Between 0 to 1.

Branching propensity 
------------------------
`branching_propensity`:`
Determine the branching probability for newly added monomers. 
The value should be similar to the experimental branching coefficient.
Between 0 to 1.

