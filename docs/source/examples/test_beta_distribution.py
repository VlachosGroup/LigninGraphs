
# set ligning path (optional if installed via pip)
import sys, os
project_path = os.path.abspath(os.path.join(os.getcwd(), '../../..'))
print(project_path)
sys.path.insert(0, project_path)

import os
from scipy.stats import beta
import matplotlib.pyplot as plt
import numpy as np

fig, ax = plt.subplots(1, 1)

a, b = 2, 7.5

r = 108900 * beta.rvs(a, b, size=1000)
# we can match the mean and variance of the distributions too
mean, var =  beta.stats(a, b, moments='mv')


ax.hist(r, density=True, histtype='stepfilled', alpha=0.2)
ax.legend(loc='best', frameon=False)
plt.show()

#
MW_lower, MW_upper = 300, 1000
n_total_population = 10000
delta_prob = beta.cdf(MW_upper/108900, a, b) - beta.cdf(MW_lower/108900, a, b)
n_simulation = delta_prob * n_total_population
