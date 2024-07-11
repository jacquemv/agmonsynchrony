# agmonsynchrony
Compute the jitter-based synchrony index between time series described in Agmon's paper [1].



### Syntax

```python
SI, pval, Nc = synchrony_index(timeseries, tau, Nc_max_exact=1000)
```

### Arguments

- **timeseries** (list or tuple of arrays): a list of two or more time series. Each time series is a sorted array of time instants.
- **tau** (float or array): maximum delay between events of different time series for the identification of coincidences. The width of the jitter window used for computing the expected number of coincidences is set to 2*tau. If tau is an array, all values are used consecutively
- **window** (str): if window is 'bilateral' (default), the original Agmon's method is applied; if window is 'unilateral', only coindidences occurring after (when tau > 0) or before (when tau < 0) the events of the reference time series are considered.
- **Nc_max_exact** (int): number of coincidences above which a normal approximation (Z score) is used for computing the p-values (default: 1000, which is generally appropriate)

### Outputs

- **SI** (array with 2 or 3 dimensions): matrix of synchrony indices. Values normally range between 0 (no synchrony), and 1 (perfect synchrony) but can also be negative (antisynchrony) down to -1. The synchrony index SI[i, j] is computed between time series i (the reference) and j (the target). If 'tau' is an array, SI[i, j, k] is the synchrony index for the k-th value of 'tau'.
- **pval** (array with the same shape as 'SI'): matrix of corresponding p-values for the significance of the synchrony
- **Nc** (int array with the same shape as 'SI'): matrix of the number of observed coincidences


### Example

Let's consider two time series with 4 and 2 samples respectively. We are looking for coincidences within a maximum delay of 0.1:
```python
from agmonsynchrony import synchrony_index
ts1 = [1, 2, 3, 4]
ts2 = [2.03, 3.95]
SI, pval, Nc = synchrony_index([ts1, ts2], tau=0.1)
```
The number of observed coincidences is:
```python
>> print(Nc)
[[4 2]
 [2 2]]
```
The matrix of synchrony indices is:
```python
>> print(SI)
[[1.   1.]
 [0.5  1.]]
```
The diagonal is 1 because time series are synchronized with themselves. The lower left value is 0.5 because when the time series 'ts2' is taken as a reference, only half of the samples of the time series 'ts1' are within 'tau' of a sample of 'ts2'. Therefore, the matrix is not symmetric.

### Installation

The package can be installed using the command ``pip install agmonsynchrony`` (on Windows, a compiler such as Microsoft Visual C++ is required).

If the code is downloaded from github, local installation on Linux is done by running ``make local`` and including the directory 'agmonsynchrony' in the PYTHONPATH environment variable.

Tested using Anaconda 2023.09 (python 3.11) on Linux and Windows.


### Acknowledgements

This work was supported by the Natural Sciences and Engineering Research
Council of Canada (NSERC grant RGPIN-2020-05252).


### References

1. A. Agmon. [A novel, jitter-based method for detecting and measuring spike synchrony and quantifying temporal firing precision](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3423071/). *Neural Syst. Circuits* 2012, vol. 2, article 5.

2. J.-P. Longpré, S. Salavatian, E. Beaumont, J. A. Armour, J. L. Ardell, V. Jacquemet. [Measure of synchrony in the activity of intrinsic cardiac neurons](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4019347/). *Physiol Meas.* 2014, vol. 35, no. 4, pp. 549–566.
