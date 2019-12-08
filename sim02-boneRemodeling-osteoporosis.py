#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Bone remodeling model

State variables:
    - xC: osteoclast
    - xB: osteoblast
    - xT: TGF
    - xW: Wnt
    
Scenarios:
    - Homeostasis
    - Osteoporosis
    - Oscillations
"""

import PyDSTool
import numpy as np
import pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import sympy as sp

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# name of system
DSargs = PyDSTool.args(name='ode')

#--- PARAMETERS ---#

# with reference
bCpar = 0.2 # komarova
bTpar = 499. # farhat
bBpar = 0.02 # komarova

# estimates
bCTpar = 2.0e-2
bWpar = 0.1 

# SCENARIO 2: osteoporosis
string_case = "osteoporosis"
string_title = "osteoporosis"

aTpar = 0.5*bTpar

aBWpar = bBpar
aWpar = bWpar*aTpar/bTpar

#--- PYDSTOOL ---#
DSargs.pars = {
'aC': bCpar + bCTpar,
'bC': bCpar,
'bCT': bCTpar,
'aBW': aBWpar,
'bB': aBWpar,
'aT': aTpar,
'bT': bTpar,
'aW': aWpar,
'bW': bWpar
}

# model equations                   
DSargs.varspecs = {
'xC': 'aC*xB^(-1.0) - bC*xC - bCT*xC*xT',
'xB': 'aBW*xB*xW - bB*xB',
'xT': 'aT*xC - bT*xT',
'xW': 'aW*xT*xC - bW*xW'}

# baseline initial conditions
DSargs.ics      = {'xC': 0.5, 
                   'xB': 0.2, 
                   'xT': 0.0,
                   'xW': 0.0}

aC = DSargs.pars['aC']
bC = DSargs.pars['bC']
bCT = DSargs.pars['bCT']

aBW = DSargs.pars['aBW']
bB = DSargs.pars['bB']

aT = DSargs.pars['aT']
bT = DSargs.pars['bT']

aW = DSargs.pars['aW']
bW = DSargs.pars['bW']

xWeq = bB/aBW
xTeq = np.sqrt(bW*aT*xWeq/(aW*bT))
xCeq = bT*xTeq/aT
xBeq = aC/(xCeq*(bC+bCT*xTeq))

print('Steady-state:')
eqs = [xCeq, xBeq, xTeq, xWeq]
print(eqs)

# time
DSargs.tdomain = [0,150]

# solve
ode  = PyDSTool.Generator.Vode_ODEsystem(DSargs)
traj = ode.compute('odeSol')
pd   = traj.sample(dt=0.1)

# plot
plt.figure(figsize=(3,1.75))
plt.plot(pd['t'], pd['xC'],color='#1f77b4',linewidth=1)
plt.plot(pd['t'], pd['xB'],color='#ff7f0e',linewidth=1)

plt.legend(['OCs','OBs'],
           ncol=5,fontsize=7,loc='upper center')
plt.xlabel('Time (days)',fontsize=10)
plt.ylabel('Density',fontsize=10)
plt.title('\\bf ' + string_title)

plt.tight_layout()

plt.savefig(string_case + "_cells.pdf")

# plot
plt.figure(figsize=(3,1.75))
plt.plot(pd['t'], pd['xT'],color='#2ca02c',linewidth=1)
plt.plot(pd['t'], pd['xW'],color='#d62728',linewidth=1)

plt.legend(['TGFb','Wnt'],
           ncol=5,fontsize=7,loc='upper center')
plt.xlabel('Time (days)',fontsize=10)
plt.ylabel('Density',fontsize=10)
plt.title('\\bf ' + string_title)

plt.tight_layout()

plt.savefig(string_case + "_mols.pdf")