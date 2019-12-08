from __future__ import division

import math

import numpy as np


# Non-monotonic Ishigami Function (3 parameters)
# First-order indices:
# x1: 0.3139
# x2: 0.4424
# x3: 0.0

# experiment 1: transient
#'aC': 1.0e0
#'bCT': 1.0e-2,
#'aBW': 2.0e-2,
#'aT': 5.0e1,
#'aW': 1.0e0,
#'bW': 5.0e-1, 
#'KC': 0.0,
#'KB': 0.9,
#'aCM': 1.5,
#'aM' : 1.0e-3,
#'KM' : 1.0,
#'aMT': 2.0e-2

#x1 aC 0.01 10.0
#x2 bCT 1.0 2.0
#x3 aBW 0.001 0.1
#x4 aT 0.001 0.1
#x5 aW 0.001 0.1
#x6 KC 0.0 0.9
#x7 KB 0.0 0.9
#x8 aCM 0.01 0.1

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
    
    xWeq = bB/aBW
    xTeq = np.sqrt(bW*aT*xWeq/(aW*bT*(1.0-KB)))
    xCeq = bT*xTeq/aT
    Y = aC*(1.0-KC)/(xCeq*(bC+bCT*xTeq-aCM*KM))
    
#    print(Y)
    
#    Y[0] = xCeq
#    Y[1] = xBeq
#    Y[2] = xTeq
#    Y[3] = xWeq

    return Y