#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Bone remodeling model

Bifurcation diagram for bT
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

# homeostasis
#[p0,p1,p2,p3,p4,p5,p6,p7,p8,p9] = [
#        2.72e-01, 1.9e-01, 8.18e-02, 
#        9.38e-03, 9.37e-03, 
#        5.88e-01, 1.51e+01, 
#        1.0, 
#        2.57e+01, 2.57e+01]

# osteoporosis
[p0,p1,p2,p3,p4,p5,p6,p7,p8,p9] = [
        2.72e-01, 1.9e-01, 8.18e-02, 
        9.38e-03, 9.37e-03, 
        1, 1.51e+01, 
        1.0, 
        2.57e+01, 2.57e+01]

# parameters
DSargs.pars = {
'aC': p0, 'bC': p1, 'bCT': p2, 
    
'aBW': p3, 'bB': p4,
    
'aT': p5, 'bT': p6, 'TG1': 0.0,
    
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
'xT': 'aT*kC*xC*BM - (bT+TG1)*xT',
'xW': 'aW*xC - bW*xW',
'BM': 'kB*xB - kC*xC*BM'}

# time
DSargs.tdomain = [0,10]

# solve
ode  = PyDSTool.Generator.Vode_ODEsystem(DSargs)
traj = ode.compute('odeSol')
pd   = traj.sample(dt=0.1)

# plot
#plt.figure(figsize=(3,1.75))
#plt.plot(pd['t'], pd['xC'],color='#1f77b4',linewidth=1)
#plt.plot(pd['t'], pd['xB'],color='#ff7f0e',linewidth=1)
#
#plt.legend(['OCs','OBs'],
#           ncol=5,fontsize=7,loc='upper center')
#plt.xlabel('Time (days)',fontsize=10)
#plt.ylabel('Density',fontsize=10)
##plt.ylim([0,15])
#plt.tight_layout()

#%%

# pars values
aC  = DSargs.pars['aC']
bC  = DSargs.pars['bC']
bCT = DSargs.pars['bCT']

aBW = DSargs.pars['aBW']
bB  = DSargs.pars['bB']

aT = DSargs.pars['aT']
bT  = DSargs.pars['bT']

aW = DSargs.pars['aW']
bW   = DSargs.pars['bW']

kB = DSargs.pars['kB']
kC = DSargs.pars['kC']

DSargs.ics      = {
'xC': (-bB*bC*bT*bW + np.sqrt(bB*bT*bW*(4*aC*aT*aW*aBW*bCT*kB + bB*bC**2*bT*bW)))/(2*aW*aBW*bT*bCT),
'xB': (-bB*bC*bT*bW + np.sqrt(bB*bT*bW*(4*aC*aT*aW*aBW*bCT*kB + bB*bC**2*bT*bW)))/(2*aT*bB*bCT*bW*kB),
'xT': (-bB*bC*bT*bW + np.sqrt(bB*bT*bW*(4*aC*aT*aW*aBW*bCT*kB + bB*bC**2*bT*bW)))/(2*bB*bT*bCT*bW),
'xW': (-bB*bC*bT*bW + np.sqrt(bB*bT*bW*(4*aC*aT*aW*aBW*bCT*kB + bB*bC**2*bT*bW)))/(2*aBW*bT*bCT*bW),
'BM': aW*aBW*bT/(aT*bB*bW*kC)}

ode  = PyDSTool.Generator.Vode_ODEsystem(DSargs)

PyCont = PyDSTool.ContClass(ode)

#---BIFURCATION PARAMETER---#
bifPar = 'TG1'

step = 0.5
PCargs = PyDSTool.args(name='EQ1', type='EP-C')

PCargs.freepars     = [bifPar]
PCargs.StepSize     = step
PCargs.MaxNumPoints = 100
PCargs.MinStepSize  = 0.01*step
PCargs.MaxStepSize  = 1*step
PCargs.FuncTol      = 1e-12
PCargs.VarTol       = 1e-12
PCargs.LocBifPoints = 'all'
PCargs.SaveEigen    = True

PyCont.newCurve(PCargs)
PyCont['EQ1'].forward()
#PyCont['EQ1'].backward()

#PyCont.computeEigen()


#%%
PyCont['EQ1'].display((bifPar,'xC'),axes=(3,1,1),stability=True)
#PyCont['LC1'].display((bifPar,'xC_min'),axes=(3,1,1),stability=True)
#PyCont['LC1'].display((bifPar,'xC_max'),axes=(3,1,1),stability=True)
plt.title('\\bf TGF$\\beta$ inhibiton')
plt.xlabel('')
plt.ylabel('$x_C^*$',fontsize=18)
plt.xlim([0,50])
plt.ylim([0,2])

PyCont['EQ1'].display((bifPar,'xB'),axes=(3,1,2),stability=True)
#PyCont['LC1'].display((bifPar,'xB_min'),axes=(3,1,2),stability=True)
#PyCont['LC1'].display((bifPar,'xB_max'),axes=(3,1,2),stability=True)
plt.title('')
plt.xlabel('')
plt.ylabel('$x_B^*$',fontsize=18)
plt.xlim([0,50])
plt.ylim([0.5,1.5])

PyCont['EQ1'].display((bifPar,'BM'),axes=(3,1,3),stability=True)
#PyCont['LC1'].display((bifPar,'BM_min'),axes=(3,1,3),stability=True)
#PyCont['LC1'].display((bifPar,'BM_max'),axes=(3,1,3),stability=True)
plt.title('')
plt.ylabel('$z^*$',fontsize=18)
plt.xlim([0,50])
plt.ylim([0,3])

plt.xlabel('$u_{tg1}$',fontsize=18)

PyCont.plot.toggleLabels('off')
PyCont.plot.togglePoints('off')
PyCont.plot.togglePoints(visible='on', bylabel='H1')
PyCont.plot.toggleLabels('on','H1')
PyCont.plot.setLabels('$H1$', bylabel='H1')

fig = plt.gcf()
fig.set_size_inches(3, 3.75)
plt.tight_layout()
plt.savefig("rem-osteoporosis-bif-bT.pdf")