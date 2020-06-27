#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Bone metastasis model

Bifurcation diagram for bT
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
[p0,p1,p2,p3,p4,p5,p6,p7,p8,p9] = [2.72e-01, 1.9e-01, 8.18e-02, 9.38e-03,
 9.37e-03, 5.88e-01, 1.51e+01, 1.0, 2.57e+01, 2.57e+01]

# osteoclast
#aCMpar = 1.0e-1
#aCMpar = 1.0e0
#aBMpar = 0.0e-3
#bBMpar  = 0.5
#aMTpar = 5.0e-2
#string_case = "osteoclastic"
#string_title = "osteoclastic lesion"

## osteoblast
bBMpar  = 0.0
aCMpar = 0.0e-1
#aCMpar = 9.0e-1
aBMpar = 2.0e-2
aMTpar = 5.0e-2
string_case = "osteoblastic"
string_title = "osteoblastic lesion"

## mixed?
#bBMpar  = 1e1
#aCMpar = 4e0
#aBMpar = 1.5e-1
#aMTpar = 5.0e-2
#string_case = "mixed"
#string_title = "mixed lesion"

# parameters
DSargs.pars = {
'aC': p0, 'bC': p1, 'bCT': p2, 
    
'aBW': p3, 'bB': p4,
    
'aT': p5, 'bT': p6,
    
'aW': p7, 'bW': 1,
    
'kB': p8, 'kC': p9,

'bBM': bBMpar,
'aCM': aCMpar,
'aBM': aBMpar,
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

kC = DSargs.pars['kC']
kB = DSargs.pars['kB']

bBM = DSargs.pars['bBM']
KM = DSargs.pars['KM']
aCM = DSargs.pars['aCM']
aBM = DSargs.pars['aBM']
aM = DSargs.pars['aM']
aMT = DSargs.pars['aMT']

xM = 0

xCeq = -(bBM*xM + bW)*(bB*bC*bT*(aBM*bBM*xM**2 + aBM*bW*xM + bBM*xM + bW) - np.sqrt(bB*bT*(bBM*xM + bW)*(4*aBM*aBW*aC*aT*aW*bCT*kB*xM + 4*aBM*aBW*aCM*aT*aW*bCT*kB*xM**2 + 4*aBW*aC*aT*aW*bCT*kB + 4*aBW*aCM*aT*aW*bCT*kB*xM + bB*bBM*bC**2*bT*xM + bB*bC**2*bT*bW))*(aBM*xM + 1))/(2*aBW*aW*bCT*bT*(aBM*xM + 1)*(aBM*bBM*xM**2 + aBM*bW*xM + bBM*xM + bW))

xBeq = (-bB*bC*bT*(aBM*bBM*xM**2 + aBM*bW*xM + bBM*xM + bW) + np.sqrt(bB*bT*(bBM*xM + bW)*(4*aBM*aBW*aC*aT*aW*bCT*kB*xM + 4*aBM*aBW*aCM*aT*aW*bCT*kB*xM**2 + 4*aBW*aC*aT*aW*bCT*kB + 4*aBW*aCM*aT*aW*bCT*kB*xM + bB*bBM*bC**2*bT*xM + bB*bC**2*bT*bW))*(aBM*xM + 1))/(2*aT*bB*bCT*kB*(aBM*bBM*xM**2 + aBM*bW*xM + bBM*xM + bW))

xTeq = (-bB*bC*bT*(aBM*bBM*xM**2 + aBM*bW*xM + bBM*xM + bW) + np.sqrt(bB*bT*(bBM*xM + bW)*(4*aBM*aBW*aC*aT*aW*bCT*kB*xM + 4*aBM*aBW*aCM*aT*aW*bCT*kB*xM**2 + 4*aBW*aC*aT*aW*bCT*kB + 4*aBW*aCM*aT*aW*bCT*kB*xM + bB*bBM*bC**2*bT*xM + bB*bC**2*bT*bW))*(aBM*xM + 1))/(2*bB*bCT*bT*(aBM*bBM*xM**2 + aBM*bW*xM + bBM*xM + bW))

xWeq = (-bB*bC*bT*(aBM*bBM*xM**2 + aBM*bW*xM + bBM*xM + bW) + np.sqrt(bB*bT*(bBM*xM + bW)*(4*aBM*aBW*aC*aT*aW*bCT*kB*xM + 4*aBM*aBW*aCM*aT*aW*bCT*kB*xM**2 + 4*aBW*aC*aT*aW*bCT*kB + 4*aBW*aCM*aT*aW*bCT*kB*xM + bB*bBM*bC**2*bT*xM + bB*bC**2*bT*bW))*(aBM*xM + 1))/(2*aBW*bCT*bT*(aBM*xM + 1)*(aBM*bBM*xM**2 + aBM*bW*xM + bBM*xM + bW))

BMeq = aBW*aW*bT*(aBM*xM + 1)/(aT*bB*kC*(bBM*xM + bW))

DSargs.ics      = {
'xC': xCeq,
'xB': xBeq,
'xT': xTeq,
'xW': xWeq,
'BM': BMeq,
'xM': 0.0}

ode  = PyDSTool.Generator.Vode_ODEsystem(DSargs)

PyCont = PyDSTool.ContClass(ode)

bifPar = 'bT'

PCargs = PyDSTool.args(name='EQ1', type='EP-C')
PCargs.freepars = [bifPar]
PCargs.FuncTol      = 1e-8
PCargs.VarTol       = 1e-8
PCargs.TestTol      = 1e-8
PCargs.MaxNumPoints = 100
PCargs.MinStepSize  = 0.1
PCargs.StepSize     = 0.2
PCargs.MaxStepSize  = 0.5

PCargs.LocBifPoints = 'all'
PCargs.SaveEigen = True

print("Calculating EQ1 curve ...")
PyCont.newCurve(PCargs)
PyCont['EQ1'].forward()
PyCont['EQ1'].backward()
#PyCont.computeEigen()

#%%
xM = 1

xCeq = -(bBM*xM + bW)*(bB*bC*bT*(aBM*bBM*xM**2 + aBM*bW*xM + bBM*xM + bW) - np.sqrt(bB*bT*(bBM*xM + bW)*(4*aBM*aBW*aC*aT*aW*bCT*kB*xM + 4*aBM*aBW*aCM*aT*aW*bCT*kB*xM**2 + 4*aBW*aC*aT*aW*bCT*kB + 4*aBW*aCM*aT*aW*bCT*kB*xM + bB*bBM*bC**2*bT*xM + bB*bC**2*bT*bW))*(aBM*xM + 1))/(2*aBW*aW*bCT*bT*(aBM*xM + 1)*(aBM*bBM*xM**2 + aBM*bW*xM + bBM*xM + bW))

xBeq = (-bB*bC*bT*(aBM*bBM*xM**2 + aBM*bW*xM + bBM*xM + bW) + np.sqrt(bB*bT*(bBM*xM + bW)*(4*aBM*aBW*aC*aT*aW*bCT*kB*xM + 4*aBM*aBW*aCM*aT*aW*bCT*kB*xM**2 + 4*aBW*aC*aT*aW*bCT*kB + 4*aBW*aCM*aT*aW*bCT*kB*xM + bB*bBM*bC**2*bT*xM + bB*bC**2*bT*bW))*(aBM*xM + 1))/(2*aT*bB*bCT*kB*(aBM*bBM*xM**2 + aBM*bW*xM + bBM*xM + bW))

xTeq = (-bB*bC*bT*(aBM*bBM*xM**2 + aBM*bW*xM + bBM*xM + bW) + np.sqrt(bB*bT*(bBM*xM + bW)*(4*aBM*aBW*aC*aT*aW*bCT*kB*xM + 4*aBM*aBW*aCM*aT*aW*bCT*kB*xM**2 + 4*aBW*aC*aT*aW*bCT*kB + 4*aBW*aCM*aT*aW*bCT*kB*xM + bB*bBM*bC**2*bT*xM + bB*bC**2*bT*bW))*(aBM*xM + 1))/(2*bB*bCT*bT*(aBM*bBM*xM**2 + aBM*bW*xM + bBM*xM + bW))

xWeq = (-bB*bC*bT*(aBM*bBM*xM**2 + aBM*bW*xM + bBM*xM + bW) + np.sqrt(bB*bT*(bBM*xM + bW)*(4*aBM*aBW*aC*aT*aW*bCT*kB*xM + 4*aBM*aBW*aCM*aT*aW*bCT*kB*xM**2 + 4*aBW*aC*aT*aW*bCT*kB + 4*aBW*aCM*aT*aW*bCT*kB*xM + bB*bBM*bC**2*bT*xM + bB*bC**2*bT*bW))*(aBM*xM + 1))/(2*aBW*bCT*bT*(aBM*xM + 1)*(aBM*bBM*xM**2 + aBM*bW*xM + bBM*xM + bW))

BMeq = aBW*aW*bT*(aBM*xM + 1)/(aT*bB*kC*(bBM*xM + bW))

PCargs.initpoint      = {
'xC': xCeq,
'xB': xBeq,
'xT': xTeq,
'xW': xWeq,
'BM': BMeq,
'xM': 1.0}

PCargs.name = 'EQ2'
PCargs.type = 'EP-C'
PCargs.LocBifPoints = 'all'
#PCargs.SaveEigen = True

PCargs.FuncTol      = 1e-8
PCargs.VarTol       = 1e-8
PCargs.TestTol      = 1e-8
PCargs.MaxNumPoints = 500
PCargs.MinStepSize  = 0.1
PCargs.StepSize     = 0.1
PCargs.MaxStepSize  = 0.1

print("Calculating EQ2 curve ...")
PyCont.newCurve(PCargs)
PyCont['EQ2'].forward()
PyCont['EQ2'].backward()
#PyCont.computeEigen()

#%%

#PCargs.name = 'EQ3' # second branch
#PCargs.initpoint = 'EQ1:BP2'
##PCargs.initdirec = PyCont['EQ1'].getSpecialPoint('BP1').labels['BP']['data'].branch
#PCargs.force = True
## Compute second branch
#PyCont.newCurve(PCargs)
#PyCont['EQ3'].forward()
#PyCont['EQ3'].backward()
#%%

little = 0.00001

xT = -aM/aMT

xCeq = (aBM*aC*aT*aW*kB - aBW*aCM*aT*aW*kB + bB*bBM*bC*bT*xT + bB*bBM*bCT*bT*xT**2 + np.sqrt(aBM**2*aC**2*aT**2*aW**2*kB**2 - 2*aBM*aBW*aC*aCM*aT**2*aW**2*kB**2 - 2*aBM*aC*aT*aW*bB*bBM*bC*bT*kB*xT - 2*aBM*aC*aT*aW*bB*bBM*bCT*bT*kB*xT**2 + 4*aBM*aCM*aT*aW*bB*bC*bT*bW*kB*xT + 4*aBM*aCM*aT*aW*bB*bCT*bT*bW*kB*xT**2 + aBW**2*aCM**2*aT**2*aW**2*kB**2 - 2*aBW*aCM*aT*aW*bB*bBM*bC*bT*kB*xT - 2*aBW*aCM*aT*aW*bB*bBM*bCT*bT*kB*xT**2 + bB**2*bBM**2*bC**2*bT**2*xT**2 + 2*bB**2*bBM**2*bC*bCT*bT**2*xT**3 + bB**2*bBM**2*bCT**2*bT**2*xT**4))/(2*(aBM+little)*aW*bT*(bC + bCT*xT)) 

xBeq =  bT*xT/(aT*kB) 

xTeq =  xT

xWeq = (aBM*aC*aT*aW*kB - aBW*aCM*aT*aW*kB - bB*bBM*bC*bT*xT - bB*bBM*bCT*bT*xT**2 + np.sqrt(aBM**2*aC**2*aT**2*aW**2*kB**2 - 2*aBM*aBW*aC*aCM*aT**2*aW**2*kB**2 - 2*aBM*aC*aT*aW*bB*bBM*bC*bT*kB*xT - 2*aBM*aC*aT*aW*bB*bBM*bCT*bT*kB*xT**2 + 4*aBM*aCM*aT*aW*bB*bC*bT*bW*kB*xT + 4*aBM*aCM*aT*aW*bB*bCT*bT*bW*kB*xT**2 + aBW**2*aCM**2*aT**2*aW**2*kB**2 - 2*aBW*aCM*aT*aW*bB*bBM*bC*bT*kB*xT - 2*aBW*aCM*aT*aW*bB*bBM*bCT*bT*kB*xT**2 + bB**2*bBM**2*bC**2*bT**2*xT**2 + 2*bB**2*bBM**2*bC*bCT*bT**2*xT**3 + bB**2*bBM**2*bCT**2*bT**2*xT**4))/(2*bT*(aBM*bC*bW + aBM*bCT*bW*xT - aBW*bBM*bC - aBW*bBM*bCT*xT)) 

BMeq =  bT*(aBM*aC*aT*aW*kB - aBW*aCM*aT*aW*kB + bB*bBM*bC*bT*xT + bB*bBM*bCT*bT*xT**2 - np.sqrt(aBM**2*aC**2*aT**2*aW**2*kB**2 - 2*aBM*aBW*aC*aCM*aT**2*aW**2*kB**2 - 2*aBM*aC*aT*aW*bB*bBM*bC*bT*kB*xT - 2*aBM*aC*aT*aW*bB*bBM*bCT*bT*kB*xT**2 + 4*aBM*aCM*aT*aW*bB*bC*bT*bW*kB*xT + 4*aBM*aCM*aT*aW*bB*bCT*bT*bW*kB*xT**2 + aBW**2*aCM**2*aT**2*aW**2*kB**2 - 2*aBW*aCM*aT*aW*bB*bBM*bC*bT*kB*xT - 2*aBW*aCM*aT*aW*bB*bBM*bCT*bT*kB*xT**2 + bB**2*bBM**2*bC**2*bT**2*xT**2 + 2*bB**2*bBM**2*bC*bCT*bT**2*xT**3 + bB**2*bBM**2*bCT**2*bT**2*xT**4))/(2*aT**2*bB*kB*kC*(aC*(bBM+little) - (aCM+little)*bW)) 

xMeq = (-aBM*aC*aT*aW*kB - aBW*aCM*aT*aW*kB + bB*bBM*bC*bT*xT + bB*bBM*bCT*bT*xT**2 + np.sqrt(aBM**2*aC**2*aT**2*aW**2*kB**2 - 2*aBM*aBW*aC*aCM*aT**2*aW**2*kB**2 - 2*aBM*aC*aT*aW*bB*bBM*bC*bT*kB*xT - 2*aBM*aC*aT*aW*bB*bBM*bCT*bT*kB*xT**2 + 4*aBM*aCM*aT*aW*bB*bC*bT*bW*kB*xT + 4*aBM*aCM*aT*aW*bB*bCT*bT*bW*kB*xT**2 + aBW**2*aCM**2*aT**2*aW**2*kB**2 - 2*aBW*aCM*aT*aW*bB*bBM*bC*bT*kB*xT - 2*aBW*aCM*aT*aW*bB*bBM*bCT*bT*kB*xT**2 + bB**2*bBM**2*bC**2*bT**2*xT**2 + 2*bB**2*bBM**2*bC*bCT*bT**2*xT**3 + bB**2*bBM**2*bCT**2*bT**2*xT**4))/(2*(aBM+little)*(aCM+little)*aT*aW*kB)

PCargs.initpoint      = {
'xC': xCeq+little,
'xB': xBeq+little,
'xT': xTeq,
'xW': xWeq+little,
'BM': BMeq+little,
'xM': xMeq+little}

PCargs.name = 'EQ3'
PCargs.type = 'EP-C'
PCargs.LocBifPoints = 'all'
PCargs.SaveEigen = True
PCargs.force = True

print("Calculating EQ3 curve ...")
PyCont.newCurve(PCargs)
PyCont['EQ3'].forward()
PyCont['EQ3'].backward()
#%%
#PyCont['EQ1'].display((bifPar,'xM'),axes=(1,1,1),stability=True)
#PyCont['EQ2'].display((bifPar,'xM'),axes=(1,1,1),stability=True)

#%%
PyCont['EQ1'].display((bifPar,'xC'),axes=(4,1,1),stability=True)
PyCont['EQ2'].display((bifPar,'xC'),axes=(4,1,1),stability=True)
PyCont['EQ3'].display((bifPar,'xC'),axes=(4,1,1),stability=True)
plt.title('\\bf TGF$\\beta$ inhibition')
plt.ylabel('$x_C^*$',fontsize=18)
plt.xlabel('')
plt.xlim([0,20])
plt.ylim([0,4])

PyCont['EQ1'].display((bifPar,'xB'),axes=(4,1,2),stability=True)
PyCont['EQ2'].display((bifPar,'xB'),axes=(4,1,2),stability=True)
PyCont['EQ3'].display((bifPar,'xB'),axes=(4,1,2),stability=True)
plt.title('')
plt.ylabel('$x_B^*$',fontsize=18)
plt.xlabel('')
plt.xlim([0,20])
plt.ylim([0,1.5])

PyCont['EQ1'].display((bifPar,'BM'),axes=(4,1,3),stability=True)
PyCont['EQ2'].display((bifPar,'BM'),axes=(4,1,3),stability=True)
PyCont['EQ3'].display((bifPar,'BM'),axes=(4,1,3),stability=True)
plt.title('')
plt.xlabel('')
plt.ylabel('$z^*$',fontsize=18)
plt.xlim([0,20])
plt.ylim([0,1.5])

PyCont['EQ1'].display((bifPar,'xM'),axes=(4,1,4),stability=True)
PyCont['EQ2'].display((bifPar,'xM'),axes=(4,1,4),stability=True)
PyCont['EQ3'].display((bifPar,'xM'),axes=(4,1,4),stability=True)
plt.title('')
plt.ylabel('$x_M^*$',fontsize=18)
plt.xlim([0,20])
plt.ylim([-0.1,1.1])

plt.xlabel('$\\beta_T$',fontsize=18)

PyCont.plot.toggleLabels('off')
PyCont.plot.togglePoints('off')
PyCont.plot.togglePoints(visible='on', bylabel='H1')

fig = plt.gcf()
fig.set_size_inches(3, 5)
plt.tight_layout()
plt.savefig("met-bif-bT-osteoblastic.pdf")