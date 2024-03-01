import numpy as np
from agmonsynchrony import synchrony_index

def test_coincidence_positive_tau():
    t1 = [1, 2, 3, 4, 5]
    t2 = [2.05, 2.99, 4.012, 5.5]
    SI, p, Nc = synchrony_index([t1, t2], 0.1, bilateral=False)
    assert np.all(Nc == [[5, 2], [1, 4]])

def test_coincidence_negative_tau():
    t1 = [1, 2, 3, 4, 5]
    t2 = [2.05, 2.99, 4.012, 5.5]
    SI, p, Nc = synchrony_index([t1, t2], -0.1, bilateral=False)
    assert np.all(Nc == [[5, 1], [2, 4]])

def _interval_count(t1, t2, tau):
    if tau == 0:
        return 0
    S = np.zeros_like(t2, np.int8)
    if tau > 0:
        for t in t1:
            S |= (t2 >= t) & (t2 < t+tau)
    else:
        for t in t1:
            S |= (t2 >= t+tau) & (t2 < t)
    return S.sum()

def test_coincidence_count():
    np.random.seed(28937456)
    t1 = np.random.random(100).cumsum()
    t2 = np.random.random(100).cumsum()
    tau = 0.3
    SI, p, Nc_pos = synchrony_index([t1, t2], tau, bilateral=False)
    SI, p, Nc_neg = synchrony_index([t1, t2], -tau, bilateral=False)
    assert Nc_pos[0, 1] ==  _interval_count(t1, t2, tau)
    assert Nc_neg[0, 1] ==  _interval_count(t1, t2, -tau)
    assert Nc_pos[1, 0] ==  _interval_count(t2, t1, tau)
    assert Nc_neg[1, 0] ==  _interval_count(t2, t1, -tau)

def test_bilateral_equivalent():
    np.random.seed(56856785)
    t1 = np.random.random(100).cumsum()
    t2 = np.random.random(100).cumsum()
    for tau in 0.3, 0.5, -0.3, -0.5:
        SI1, p1, Nc1 = synchrony_index([t1, t2], tau, bilateral=False)
        SI2, p2, Nc2 = synchrony_index([t1, t2-tau/2], abs(tau)/2)
        assert np.isclose(SI1[0, 1], SI2[0, 1])
        assert np.isclose(p1[0, 1], p2[0, 1])
        assert Nc1[0, 1] == Nc2[0, 1]
        SI2, p2, Nc2 = synchrony_index([t1-tau/2, t2], abs(tau)/2)
        assert np.isclose(SI1[1, 0], SI2[1, 0])
        assert np.isclose(p1[1, 0], p2[1, 0])
        assert Nc1[1, 0] == Nc2[1, 0]
