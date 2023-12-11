import numpy as np
from agmonsynchrony import synchrony_index

def load_timeseries(n):
    return np.fromfile(f'../data/timeseries{n}.txt', sep=' ')

T = [load_timeseries(n) for n in range(1, 5)]


def test_empty_list():
    SI, p, Nc = synchrony_index([], 0.03)
    assert SI.size == 0 and p.size == 0 and Nc.size == 0

def test_one_time_series():
    SI, p, Nc = synchrony_index([T[0]], 0.03)
    assert SI == 1
    assert p == 0
    assert Nc == T[0].size

def test_empty_time_series():
    SI, p, Nc = synchrony_index([T[0], np.empty(0)], 0.03)
    assert np.all(SI == [[1, 0], [0, 0]])
    assert np.all(p == [[0, 1], [1, 1]])
    assert np.all(Nc == [[len(T[0]), 0], [0, 0]])

def test_integer_time_series():
    SI, p, Nc = synchrony_index([T[0], np.array([1, 2, 3, 4])], 0.03)

def test_list_time_series():
    SI, p, Nc = synchrony_index([[2, 4.1], [1, 2, 3, 4]], 0.03)

def test_noncontiguous_array():
    T = np.linspace(0, 50, 200)
    SI, p, Nc = synchrony_index([T[::2], [1, 2, 3, 4]], 0.03)

def test_known_SI():
    tau = np.array([0.01, 0.02, 0.03])
    SI, p, Nc = synchrony_index(T, tau[2])
    shape = len(T), len(T)
    assert SI.shape == shape and p.shape == shape and Nc.shape == shape
    assert np.isclose(SI[0, 1], 0.23)
    assert np.isclose(SI[2, 3], 0.47)

def test_multiple_tau():
    tau_list = [0.01, 0.02, 0.03]
    SI, p, Nc = synchrony_index(T, tau_list)
    shape = len(T), len(T), len(tau_list)
    assert SI.shape == shape and p.shape == shape and Nc.shape == shape
