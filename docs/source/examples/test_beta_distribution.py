
# set ligning path (optional if installed via pip)
import sys, os
project_path = os.path.abspath(os.path.join(os.getcwd(), '../../..'))
print(project_path)
sys.path.insert(0, project_path)

import os
from scipy.stats import beta
import matplotlib.pyplot as plt
import numpy as np
import random
import pandas as pd

fig, ax = plt.subplots(1, 1)

df = pd.read_csv('mw_experimental.csv')
mw_exp = np.reshape(np.array(df),(df.size))
loc = np.floor(np.min(mw_exp)) #147
scale = np.ceil(np.max(mw_exp) - loc) #24907

#a, b, _, _ = beta.fit(mw_exp, floc = loc, fscale = scale)
a, b = 2, 7.5
r = beta.rvs(a, b, loc = loc, scale = scale, size=1000)

ax.hist(r, alpha=0.2, bins=40, color='b')
#ax.hist(mw_exp, density=True,  alpha=0.2, bins=40, color='r')

ax.legend(loc='best', frameon=False)
plt.show()

#
MW_lower, MW_upper = 10000, 15000
n_total_population = 10000
delta_prob = beta.cdf(MW_upper, a, b, loc, scale) - beta.cdf(MW_lower, a, b, loc, scale)
n_simulation = delta_prob * n_total_population
stop_size = np.random.uniform(300, 1000)
