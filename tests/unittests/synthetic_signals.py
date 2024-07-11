import numpy as np
import matplotlib.pyplot as plt
from agmonsynchrony import synchrony_index

ts1 = np.arange(10000)
ts2 = np.sort(ts1 + 0.1+ 0.04 * np.random.randn(ts1.size))
tau = np.linspace(-0.5, 0.5, 301)

# check the symmetry with respect to the reflection tau -> -tau

idx = 0, 1 # = reference, target
SI, pval, Nc = synchrony_index([ts1, ts2], tau=tau, window='bilateral')
SI[pval > 0.05] = np.nan
plt.plot(tau, SI[idx], c='0.6')

SI, pval, Nc = synchrony_index([ts1, ts2], tau=tau, window='unilateral')
SI[pval > 0.05] = np.nan
plt.plot(tau, SI[idx], '--')

idx = 1, 0 # = reference, target
SI, pval, Nc = synchrony_index([ts1, ts2], tau=tau, window='unilateral')
SI[pval > 0.05] = np.nan
plt.plot(-tau, SI[idx], '--')

plt.grid()
plt.show()