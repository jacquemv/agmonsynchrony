from agmonsynchrony import synchrony_index

ts1 = [1, 2, 3, 4]
ts2 = [2.03, 3.95]
SI, pval, Nc = synchrony_index([ts1, ts2], tau=0.1)
print(SI)
print(pval)
print(Nc)