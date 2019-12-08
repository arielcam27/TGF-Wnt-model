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
    bC = 1.0e-3
    bCT = values[:,1]
    
    aBW = values[:,2]
    bB = values[:,3]
    
    aT = values[:,4]
    bT = 1.0
    
    aW = values[:,5]
    bW = 1.0
    
    KC = values[:,6]
    KB = values[:,7]
    aCM = values[:,8]
    
    KM = 1.0
    
    Y = bB/aBW
    
    KM = 1.0
    
    xWeq = bB/aBW
    xTeq = np.sqrt(bW*aT*xWeq/(aW*bT*(1.0-KB)))
    Y = bT*xTeq/aT
#    xBeq = aC/(xCeq*(bC+bCT*xTeq))
    
#    print(Y)
    
#    Y[0] = xCeq
#    Y[1] = xBeq
#    Y[2] = xTeq
#    Y[3] = xWeq

    return Y