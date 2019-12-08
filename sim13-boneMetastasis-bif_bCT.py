#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Bone metastasis model

Bifurcation diagram for bCT
"""

import PyDSTool
import numpy as np
import pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import sympy as sp

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

DSargs = PyDSTool.args(name='ode')

bCpar = 0.8 
bTpar = 499. 
bBpar = 0.48
bCTpar = 2.0e-2
bWpar = 0.1 
aBWpar = bBpar
aTpar = bTpar
aWpar = bWpar*aTpar/bTpar
string_case = "dummie"
string_title = "dummie"
KCpar = 0.0
KBpar = 0.0
aCMpar = 1.0e-1
aMTpar = 3.0e-2

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

# baseline initial conditions
# homeo-osteo
DSargs.ics      = {'xC': 0.5, 
                   'xB': 0.2, 
                   'xT': 0.0,
                   'xW': 0.0,
                   'xM': 1.0e-1}

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
DSargs.tdomain = [0,1]

DSargs.ics      = {
'xC': xCeq,
'xB': xBeq,
'xT': xTeq,
'xW': xWeq,
'xM': 0.0}

ode  = PyDSTool.Generator.Vode_ODEsystem(DSargs)

PyCont = PyDSTool.ContClass(ode)
# bifurcation parameter
bifPar = 'bCT'

PCargs = PyDSTool.args(name='EQ1', type='EP-C')
PCargs.freepars = [bifPar]
PCargs.StepSize = 1e-2
PCargs.MaxNumPoints = 100
#PCargs.MaxNumPoints = 5000
PCargs.MaxStepSize = 1e-2
PCargs.LocBifPoints = 'all'
PCargs.SaveEigen = True

print("Calculating EQ1 curve ...")
PyCont.newCurve(PCargs)
PyCont['EQ1'].forward()
PyCont['EQ1'].backward()

PCargs.name = 'LC1'
PCargs.type = 'LC-C'
PCargs.initpoint = 'EQ1:H1'
PCargs.StepSize = 1.e-1
PCargs.MinStepSize = 1.e-1
PCargs.MaxStepSize = 1.e-1
PCargs.force = True
PCargs.NumSPOut = 10
PCargs.verbosity = 1
PCargs.SolutionMeasures = 'all'
PCargs.LocBifPoints = 'all'
PCargs.FuncTol = 1e-1
PCargs.VarTol = 1e-1
PCargs.TestTol = 1e-1
PCargs.MaxNumPoints = 500
PCargs.SaveEigen = True

print("Calculating EQ1:LC1 curve ...")
PyCont.newCurve(PCargs)
PyCont['LC1'].backward()
PyCont['LC1'].forward()

KC = DSargs.pars['KC']
KB = DSargs.pars['KB']
KM = DSargs.pars['KM']
aCM = DSargs.pars['aCM']

xWeq2 = bB/aBW
xTeq2 = np.sqrt(bW*aT*xWeq2/(aW*(1.0-KB)*bT))
xCeq2 = bT*xTeq2/aT
xBeq2 = aC*(1.0-KC)/(xCeq2*(bC+bCT*xTeq2-aCM*KM))

PCargs.initpoint      = {
'xC': xCeq2,
'xB': xBeq2,
'xT': xTeq2,
'xW': xWeq2,
'xM': 1.0}

PCargs.name = 'EQ2'
PCargs.type = 'EP-C'
PCargs.freepars = [bifPar]
PCargs.StepSize = 1e-2
PCargs.MaxNumPoints = 100
PCargs.MaxStepSize = 1e-2
PCargs.LocBifPoints = 'all'
PCargs.SaveEigen = True

print("Calculating EQ2 curve ...")
PyCont.newCurve(PCargs)
PyCont['EQ2'].forward()
PyCont['EQ2'].backward()

PCargs.name = 'LC2'
PCargs.type = 'LC-C'
PCargs.initpoint = 'EQ2:H1'
PCargs.StepSize = 0.2
PCargs.MinStepSize = 0.2
PCargs.MaxStepSize = 0.2
PCargs.force = True
PCargs.NumSPOut = 10
PCargs.verbosity = 1
PCargs.SolutionMeasures = 'all'
PCargs.LocBifPoints = 'all'
PCargs.FuncTol = 1e-2
PCargs.VarTol = 1e-2
PCargs.TestTol = 1e-2
PCargs.MaxNumPoints = 6000
PCargs.SaveEigen = True

print("Calculating EQ2:LC2 curve ...")
PyCont.newCurve(PCargs)
PyCont['LC2'].backward()
PyCont['LC2'].forward()


PyCont['EQ1'].display((bifPar,'xC'),axes=(2,1,1),stability=True)
PyCont['LC1'].display((bifPar,'xC_min'),axes=(2,1,1),stability=True)
PyCont['LC1'].display((bifPar,'xC_max'),axes=(2,1,1),stability=True)
PyCont['EQ2'].display((bifPar,'xC'),axes=(2,1,1),stability=True)
PyCont['LC2'].display((bifPar,'xC_min'),axes=(2,1,1),stability=True)
PyCont['LC2'].display((bifPar,'xC_max'),axes=(2,1,1),stability=True)
plt.title('\\bf TGF$\\beta$ receptor inhib.')
plt.ylabel('$x_C^*$',fontsize=18)
plt.xlim([-0.1,0.1])
plt.ylim([0.,3.])

PyCont['EQ1'].display((bifPar,'xB'),axes=(2,1,2),stability=True)
PyCont['LC1'].display((bifPar,'xB_min'),axes=(2,1,2),stability=True)
PyCont['LC1'].display((bifPar,'xB_max'),axes=(2,1,2),stability=True)
PyCont['EQ2'].display((bifPar,'xB'),axes=(2,1,2),stability=True)
PyCont['LC2'].display((bifPar,'xB_min'),axes=(2,1,2),stability=True)
PyCont['LC2'].display((bifPar,'xB_max'),axes=(2,1,2),stability=True)
plt.title('')
plt.ylabel('$x_B^*$',fontsize=18)
plt.xlim([-0.1,0.1])
plt.ylim([0.,3.])

plt.xlabel('$\\beta_{CT}$',fontsize=18)

PyCont.plot.toggleLabels('off')
PyCont.plot.togglePoints('off')
PyCont.plot.togglePoints(visible='on', bylabel='H1')

fig = plt.gcf()
fig.set_size_inches(3, 2.5)
plt.tight_layout()
