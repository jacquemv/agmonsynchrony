import numpy as np
import math

cdef extern from "windowoverlap.h":
    int time_windows_c(double* x, int n, double tau, int bilateral,
                       double* interv)
    void coincidence_distribution_c(double* p, int n, double* P)
    void overlap_c(double* x, int na, double* y, int n, double tauj, int& S, 
                   double* p)


#------------------------------------------------------------------------------
def time_windows(double[::1] t, double tau, int bilateral):
    if t.size == 0:
        return np.empty((0, 2))
    interv = np.empty(t.size*2, dtype=np.float64)
    cdef double[::1] interv_ptr = interv
    m = time_windows_c(&t[0], t.size, tau, bilateral, &interv_ptr[0])
    return interv[:2*m].reshape((-1, 2))

#------------------------------------------------------------------------------
def coincidence_distribution(double[::1] p):
    #   Thomas M, Taub A: Calculating binomial probabilities when the trial
    #   probabilities are unequal. J Statist Comput Simulation 1982,14:125-131
    if p.size == 0:
        return np.ones(1, dtype=np.float64)
    P = np.zeros(p.size+1, dtype=np.float64)
    cdef double[::1] P_ptr = P
    coincidence_distribution_c(&p[0], p.size, &P_ptr[0])
    return P

#------------------------------------------------------------------------------
def synchrony_index(double[::1] x, double[::1] y, double tau_s, double tau_j, 
                    int bilateral, int N_exact, jitter=None):
    if tau_s == 0:
        return 0., 1., 0
    if jitter is None:
        jitter = time_windows(x, tau_s, bilateral).ravel()
    p = np.zeros(y.size, dtype=np.float64)
    cdef:
        double[::1] jitter_ptr = jitter
        double[::1] p_ptr = p
        int Nc=0, n1
        double Nc_jitter, std_Nc_jitter, Z, beta
    overlap_c(&jitter_ptr[0], jitter.size // 2, &y[0], y.size, tau_j, Nc, 
              &p_ptr[0])
    Nc_jitter = p.sum()
    n1 = y.size
    if n1 == 0:
        SI, pval = 0., 1.
    else:
        if bilateral:
            beta = min(2, tau_j/(tau_j-abs(tau_s)))
        else:
            beta = min(2, tau_j/(tau_j-abs(tau_s)/2))
        SI = beta * (Nc - Nc_jitter) / n1
        pnz = p[p > 0]
        if pnz.size > N_exact:
            std_Nc_jitter = np.sqrt(np.sum(pnz*(1-pnz)))
            Z = (Nc-Nc_jitter) / std_Nc_jitter
            pval = 0.5 * math.erfc(math.fabs(Z) / math.sqrt(2))
        else:
            P = coincidence_distribution(pnz)
            if Nc > Nc_jitter:
                # P[k] = Prob{Nc=k}; pval = Prob{k>Nc)
                pval = np.sum(P[Nc+1:])
            else:
                pval = np.sum(P[:Nc-1])
    return SI, pval, Nc

#------------------------------------------------------------------------------
def synchrony_index_matrix(T, double tau_s, double tau_j, int bilateral,
                           int N_exact):
    cdef:
        int n = len(T)
        int i, j
    SI = np.zeros((n, n), dtype=np.float64)
    pval = np.ones((n, n), dtype=np.float64)
    Nc = np.zeros((n, n), dtype=np.int32)
    for i in range(n):
        jitter = time_windows(T[i], tau_s, bilateral).ravel()
        for j in range(n):
            if T[i].size == 0 or T[j].size == 0:
                continue
            if i == j:
                SI[i, j] = 1
                pval[i, j] = 0
                Nc[i, j] = T[i].size
                continue
            out = synchrony_index(T[i], T[j], tau_s, tau_j,
                                  bilateral, N_exact, jitter)
            SI[i, j], pval[i, j], Nc[i, j] = out
    return SI, pval, Nc

#------------------------------------------------------------------------------
def synchrony_index_matrix_multiple_tau(T, list_tau_s, list_tau_j, 
                                        int bilateral, int N_exact):
    cdef:
        int n = len(T)
        int m = len(list_tau_s)
        int i, j
        double tau_s
    SI = np.zeros((n, n, m), dtype=np.float64)
    pval = np.ones((n, n, m), dtype=np.float64)
    Nc = np.zeros((n, n, m), dtype=np.int32)
    for k in range(m):
        tau_s = list_tau_s[k]
        tau_j = list_tau_j[k]
        for i in range(n):
            jitter = time_windows(T[i], tau_s, bilateral).ravel()
            for j in range(n):
                if T[i].size == 0 or T[j].size == 0:
                    continue
                if i == j:
                    SI[i, j, k] = 1
                    pval[i, j] = 0
                    Nc[i, j, k] = T[i].size
                    continue
                out = synchrony_index(T[i], T[j], tau_s, tau_j, 
                                      bilateral, N_exact, jitter)
                SI[i, j, k], pval[i, j, k], Nc[i, j, k] = out
    return SI, pval, Nc
