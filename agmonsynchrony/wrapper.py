import numpy as np
from .synchrony import (synchrony_index_matrix, 
                        synchrony_index_matrix_multiple_tau)

#------------------------------------------------------------------------------
def synchrony_index(timeseries, tau, window='bilateral', Nc_max_exact=1000):
    """Compute the jitter-based synchrony index based on Agmon's paper:
    https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3423071/

    Args:
        timeseries (list or tuple of arrays): list of time series
        tau (float or array): maximum delay between events of different time 
            series for the identification of coincidences;
            the width of the jitter window used for computing the expected 
            number of coincidences is set to 2*tau;
            if tau is an array, all values are used consecutively
        window (str): if window is 'bilateral', the original Agmon's method is 
            applied; if window is 'unilateral', only coindidences occurring 
            after (when tau > 0) or before (when tau < 0) the events of the 
            reference time series are considered
        Nc_max_exact (int): number of coincidences above which a normal 
            approximation (Z score) is used for computing the p-values 
            (default: 1000)
    
    Returns:
        SI (array with 2 or 3 dimensions): matrix of synchrony indices; 
            values normally range between 0 (no synchrony), and 1 (perfect 
            synchrony) but can also be negative (antisynchrony) down to -1;
            SI[i, j] = synchrony index between time series i (the reference)
            and j (the target); if tau is an array, SI[i, j, k] is the 
            synchrony index for the k-th value of tau
        pval (array with the same shape as SI): matrix of corresponding 
            p-values for the significance of the synchrony
        Nc (int array with the same shape as SI): matrix of the number of
            observed coincidences
    """
    to_array = lambda x: np.ascontiguousarray(x, dtype=np.float64)
    timeseries = [to_array(ts) for ts in timeseries]
    bilateral = window.startswith('bi')
    if hasattr(tau, '__iter__'): # tuple, list or array
        tau_s = to_array(tau)
        if bilateral:
            tau_s = np.abs(tau_s)
            tau_j = 2 * tau_s
        else:
            tau_j = np.abs(tau_s)
        return synchrony_index_matrix_multiple_tau(
            timeseries, tau_s, tau_j, bilateral, N_exact=Nc_max_exact
        )
    else:
        tau_s = float(tau)
        if bilateral:
            tau_s = abs(tau_s)
            tau_j = 2 * tau_s
        else:
            tau_j = abs(tau_s)
        return synchrony_index_matrix(
            timeseries, tau_s, tau_j, bilateral, N_exact=Nc_max_exact
        )