#!/usr/bin/env python3

import statsmodels.api as sm
import numpy as np

a = np.ndarray(shape=(2,2,3), dtype = float)
for i in range(3):
    for j in range(2):
        for k in range(2):
            a[j,k,i] = (100.0*i) + (10.0*j) + (1.0*k)
print(a)

b = sm.stats.StratifiedTable(a)
#b = sm.stats.contingency_tables.StratifiedTable(a)
print(b)

c = b.summary()
print(c)
print(dir(c))
print(dir(b.test_null_odds()))
print(b.test_null_odds().pvalue)
