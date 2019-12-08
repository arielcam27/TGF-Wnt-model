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

bCpar = 0.8
bTpar = 499.
bBpar = 0.48
bCTpar = 2.0e-2
bWpar = 0.1 

aBWpar = bBpar
aTpar = bTpar
aWpar = bWpar*aTpar/bTpar

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
DSargs.ics      = {'xC': 0.8, 
                   'xB': 0.8, 
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
#plt.ylim([0,15])
plt.tight_layout()

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

eqs = [xCeq, xBeq, xTeq, xWeq]
print(eqs)

DSargs.ics      = {
'xC': xCeq,
'xB': xBeq,
'xT': xTeq,
'xW': xWeq}

ode  = PyDSTool.Generator.Vode_ODEsystem(DSargs)

PyCont = PyDSTool.ContClass(ode)

# bifurcation parameter
bifPar = 'bT'

PCargs = PyDSTool.args(name='EQ1', type='EP-C')
PCargs.freepars = [bifPar]
PCargs.StepSize = 1e-1
PCargs.MaxNumPoints = 3000
PCargs.MaxStepSize = 1e-1
PCargs.LocBifPoints = 'all'
PCargs.SaveEigen = True

print("Calculating EQ1 curve ...")
PyCont.newCurve(PCargs)
PyCont['EQ1'].forward()
PyCont['EQ1'].backward()

#%%
PCargs.name = 'LC1'
PCargs.type = 'LC-C'
PCargs.initpoint = 'EQ1:H1'
PCargs.StepSize = 0.1
PCargs.MinStepSize = 0.1
PCargs.MaxStepSize = 0.1
PCargs.force = True
PCargs.NumSPOut = 10
PCargs.verbosity = 1
PCargs.SolutionMeasures = 'all'
PCargs.LocBifPoints = 'all'
PCargs.FuncTol = 1e-3
PCargs.VarTol = 1e-3
PCargs.TestTol = 1e-3
PCargs.MaxNumPoints = 5000
PCargs.SaveEigen = True

print("Calculating EQ1:LC1 curve ...")
PyCont.newCurve(PCargs)
PyCont['LC1'].backward()
PyCont['LC1'].forward()

#%%
# cancer-invasion bifurcation curve
plt.figure(figsize=(4,2));

PyCont['EQ1'].display((bifPar,'xC'),axes=(2,1,1),stability=True)
PyCont['LC1'].display((bifPar,'xC_min'),axes=(2,1,1),stability=True)
PyCont['LC1'].display((bifPar,'xC_max'),axes=(2,1,1),stability=True)
plt.title('\\bf TGF$\\beta$ inhibiton')
plt.xlabel('')
plt.ylabel('$x_C^*$',fontsize=12)
plt.xlim([200,800])
#plt.ylim([0,10])

PyCont['EQ1'].display((bifPar,'xB'),axes=(2,1,2),stability=True)
PyCont['LC1'].display((bifPar,'xB_min'),axes=(2,1,2),stability=True)
PyCont['LC1'].display((bifPar,'xB_max'),axes=(2,1,2),stability=True)
plt.title('')
plt.ylabel('$x_B^*$',fontsize=12)
plt.xlim([200,800])
#plt.ylim([0,10])

plt.xlabel('$\\beta_{T}$')

PyCont.plot.toggleLabels('off')
PyCont.plot.togglePoints('off')
PyCont.plot.togglePoints(visible='on', bylabel='H1')
PyCont.plot.toggleLabels('on','H1')
PyCont.plot.setLabels('$H1$', bylabel='H1')

fig = plt.gcf()
fig.set_size_inches(3, 2.5)
plt.tight_layout()
plt.savefig("v2_rem_bif_bT.pdf")
