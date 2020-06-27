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

#--- PARAMETERS ---#

# SCENARIO 1: HOMEOSTASIS
string_case = "homeostasis"
string_title = "homeostasis"

# we must give a name
DSargs = PyDSTool.args(name='ode')

#[p0,p1,p2,p3,p4,p5,p6,p7,p8,p9] = [2.72228935e-01, 1.90347758e-01, 8.18795014e-02, 9.38001758e-03,
# 9.37995984e-03, 5.88839726e-01, 1.51564419e+01, 9.99987690e-01, 2.57395030e+01, 2.57393446e+01]

[p0,p1,p2,p3,p4,p5,p6,p7,p8,p9] = [
        2.72e-01, 1.9e-01, 8.18e-02, 
        9.38e-03, 9.37e-03, 
        5.88e-01, 1.51e+01, 
        1.0, 
        2.57e+01, 2.57e+01]


# parameters
DSargs.pars = {
'aC': p0, 'bC': p1, 'bCT': p2, 
    
'aBW': p3, 'bB': p4,
    
'aT': p5, 'bT': p6,
    
'aW': p7, 'bW': 1,
    
'kB': p8, 'kC': p9
}

# initial conditions
DSargs.ics      = {'xC': 0.9, 
                   'xB': 0.9, 
                   'xT': 0.9,
                   'xW': 0.9,
                   'BM': 0.9}

# model equations                   
DSargs.varspecs = {
'xC': 'aC*xT/xB - bC*xC - bCT*xT*xC',
'xB': 'aBW*xB*xW/xT - bB*xB',
'xT': 'aT*kC*xC*BM - bT*xT',
'xW': 'aW*xC - bW*xW',
'BM': 'kB*xB - kC*xC*BM'}

# numerical simulation
DSargs.tdomain = [0,250]

ode  = PyDSTool.Generator.Vode_ODEsystem(DSargs)
traj = ode.compute('odeSol')

pd   = traj.sample(dt=0.1)

#plt.legend(['OCs','OBs','TGF-b','Wnt'])
#plt.xlabel('Time (days)',fontsize=20)
#plt.ylabel('Density',fontsize=20)

# time
DSargs.tdomain = [0,250]

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