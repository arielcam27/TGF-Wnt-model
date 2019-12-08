#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Bone metastasis model

State variables:
    - xC: osteoclast
    - xB: osteoblast
    - xT: TGF
    - xW: Wnt
    - xM: cancer cell
    
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

DSargs = PyDSTool.args(name='ode')

#--- HOMEOSTASIS
bCpar = 0.2 
bTpar = 499.
bBpar = 0.02
bCTpar = 2.0e-2
bWpar = 0.1 

aTpar = bTpar
aBWpar = bBpar
aWpar = bWpar*aTpar/bTpar

string_case = "osteoclastic1"
string_title = "osteoclastic lesion1"

KCpar = 0.0
KBpar = 0.9
aCMpar = 1.0e-1
aMTpar = 3.0e-2

KCpar = 0.9
KBpar = 0.0
aCMpar = 1.0e-1
aMTpar = 3.0e-2

# PYDSTOOL
DSargs.pars = {
'aC': bCpar + bCTpar,
'bC': bCpar,
'bCT': bCTpar,
'aBW': aBWpar,
'bB': aBWpar,
'aT': aTpar,
'bT': bTpar,
'aW': aWpar,
'bW': bWpar,
'KC': KCpar,
'KB': KBpar,
'aCM': aCMpar,
'aM' : 1.0e-3,
'KM' : 1.0,
'aMT': aMTpar
}

# model equations                   
DSargs.varspecs = {
'xC': '(1.0 - KC*xM/KM)*aC*xB^(-1.0) - bC*xC - bCT*xC*xT + aCM*xM',
'xB': 'aBW*xB*xW - bB*xB',
'xT': 'aT*xC - bT*xT',
'xW': '(1.0 - KB*xM/KM)*aW*xT*xC - bW*xW',
'xM': 'xM*(aM + aMT*xT)*(1.0 - xM/KM)'}

DSargs.ics      = {'xC': 1.0, 
                   'xB': 1.0, 
                   'xT': 1.0,
                   'xW': 1.0,
                   'xM': 0.1}

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

KC = DSargs.pars['KC']
KB = DSargs.pars['KB']
KM = DSargs.pars['KM']
aCM = DSargs.pars['aCM']

xWeq2 = bB/aBW
xTeq2 = np.sqrt(bW*aT*xWeq2/(aW*(1.0-KB)*bT))
xCeq2 = bT*xTeq2/aT
xBeq2 = aC*(1-KC)/(xCeq2*(bC+bCT*xTeq2)-aCM)

eqs  = [xCeq,  xBeq,  xTeq,  xWeq,  0.0]
eqs2 = [xCeq2, xBeq2, xTeq2, xWeq2, 1.0]
print(eqs)
print(eqs2)

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
plt.plot(pd['t'], pd['xM'],color='#9467bd',linewidth=1)

plt.legend(['OCs','OBs','CCs'],
           ncol=5,fontsize=7,loc='upper center')
plt.xlabel('Time (days)',fontsize=10)
plt.ylabel('Density',fontsize=10)
plt.ylim([0,2.5])
plt.title('\\bf ' + string_title)

plt.tight_layout()

plt.savefig(string_case + "_cells_cancer.pdf")

# plot
plt.figure(figsize=(3,1.75))
plt.plot(pd['t'], pd['xT'],color='#2ca02c',linewidth=1)
plt.plot(pd['t'], pd['xW'],color='#d62728',linewidth=1)
#plt.plot(pd['t'], pd['xM'],color='#9467bd',linewidth=1)

plt.legend(['TGFb','Wnt'],
           ncol=5,fontsize=7,loc='upper center')
plt.xlabel('Time (days)',fontsize=10)
plt.ylabel('Density',fontsize=10)
#plt.ylim([0,15])
plt.title('\\bf ' + string_title)

plt.tight_layout()

plt.savefig(string_case + "_mols_cancer.pdf")
