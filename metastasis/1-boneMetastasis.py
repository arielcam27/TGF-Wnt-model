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

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

DSargs = PyDSTool.args(name='ode')



#--- HOMEOSTASIS
[p0,p1,p2,p3,p4,p5,p6,p7,p8,p9] = [
        2.72e-01, 1.9e-01, 8.18e-02, 
        9.38e-03, 9.37e-03, 
        5.88e-01, 1.51e+01, 
        1.0, 
        2.57e+01, 2.57e+01]

# osteoclast
#aCMpar = 1.0e0
#aBMpar = 0.0e-3
#bBMpar  = 0.5
#aMTpar = 5.0e-2
#string_case = "osteoclastic"
#string_title = "osteoclastic lesion"

# osteoblast
#aCMpar = 0.0e-1
#aBMpar = 2.0e-2
#bBMpar  = 0.0
#aMTpar = 5.0e-2
#string_case = "osteoblastic"
#string_title = "osteoblastic lesion"

## mixed?
aCMpar = 4e0
aBMpar = 1.5e-1
bBMpar  = 1e1
aMTpar = 5.0e-2
string_case = "mixed"
string_title = "mixed lesion"

# parameters
DSargs.pars = {
'aC': p0, 'bC': p1, 'bCT': p2, 
    
'aBW': p3, 'bB': p4,
    
'aT': p5, 'bT': p6,
    
'aW': p7, 'bW': 1,
    
'kB': p8, 'kC': p9,

'aCM': aCMpar,
'aBM': aBMpar,
'bBM': bBMpar,
'aM' : -5.0e-2,
'KM' : 1.0,
'aMT': aMTpar
}

# model equations                   
DSargs.varspecs = {
'xC': '(aC + aCM*xM)*xT/xB - bC*xC - bCT*xT*xC',
'xB': '(aBW + aBM*xM)*xB*xW/xT- bB*xB' ,
'xT': 'aT*kC*xC*BM - bT*xT',
'xW': 'aW*xC - (bW + bBM*xM)*xW',
'BM': 'kB*xB - kC*xC*BM',
'xM': 'xM*(aM + aMT*xT)*(1.0 - xM/KM)'
}

DSargs.ics      = {'xC': 1.0, 
                   'xB': 1.0, 
                   'xT': 1.0,
                   'xW': 1.0,
                   'BM': 1.0,
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

bBM = DSargs.pars['bBM']
KM = DSargs.pars['KM']
aCM = DSargs.pars['aCM']

# time
DSargs.tdomain = [0,400]

# solve
ode  = PyDSTool.Generator.Vode_ODEsystem(DSargs)
traj = ode.compute('odeSol')
pd   = traj.sample(dt=0.1)

# plot
plt.close('all')
plt.figure(figsize=(3,1.75))
plt.plot(pd['t'], pd['xC'],color='#1f77b4',linewidth=1)
plt.plot(pd['t'], pd['xB'],color='#ff7f0e',linewidth=1)
plt.plot(pd['t'], pd['xM'],color='#ba34eb',linewidth=1)
plt.plot([0,400], [1,1],color='#000000',linestyle='dashed',linewidth=1)
         
plt.legend(['OCs','OBs','CCs'],
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

# plot
plt.figure(figsize=(3,1.75))
plt.plot(pd['t'], pd['BM'],color='#000000',linewidth=1)

plt.legend(['BM'],
           ncol=5,fontsize=7,loc='best')
plt.xlabel('Time (days)',fontsize=10)
plt.ylabel('Density',fontsize=10)
plt.title('\\bf ' + string_title)

plt.tight_layout()

plt.savefig(string_case + "_bm.pdf")