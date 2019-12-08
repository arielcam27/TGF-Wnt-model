from __future__ import division

import math

import numpy as np


# Non-monotonic Ishigami Function (3 parameters)
# First-order indices:
# x1: 0.3139
# x2: 0.4424
# x3: 0.0

def evaluate(values):
    Y = np.zeros([values.shape[0]])
    
    aC = values[:,0]
    bC = values[:,1]
    bCT = values[:,2]
    
    aBW = values[:,3]
    bB = values[:,4]
    
    aT = values[:,5]
    bT = values[:,6]
    
    aW = values[:,7]
    bW = values[:,8]
    
    xWeq = bB/aBW
    xTeq = np.sqrt(bW*aT*xWeq/(aW*bT))
    Y = bT*xTeq/aT
#    xBeq = aC/(xCeq*(bC+bCT*xTeq))
    
#    print(Y)
    
#    Y[0] = xCeq
#    Y[1] = xBeq
#    Y[2] = xTeq
#    Y[3] = xWeq

    return Y