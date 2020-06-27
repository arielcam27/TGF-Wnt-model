import numpy as np
import scipy as sc
import pylab as plt
import csv
from scipy import integrate
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#---------#
# initial time
t0 = 0.0
# final time
tf = 250.
# sub-intervals
N = 10
# time windows
time = np.linspace(t0, tf, int((tf-t0)*N))


#--- HOMEOSTASIS
[p0,p1,p2,p3,p4,p5,p6,p7,p8,p9] = [
        2.72e-01, 1.9e-01, 8.18e-02, 
        9.38e-03, 9.37e-03, 
        5.88e-01, 1.51e+01, 
        1.0, 
        2.57e+01, 2.57e+01]

# osteoclast
aCMpar = 1.0e0
aBMpar = 0.0e-3
bBMpar  = 0.5
aMTpar = 5.0e-2
aMpar = -5.0e-2
KMpar = 1
# day ~200
xC0 = 2.2
xB0 = 1.4
xT0 = 1.4
xW0 = 1.8
BM0 = 0.6
xM0 = 0.4
string_case = "osteoclastic"
string_title = "osteoclastic lesion"

## osteoblast
#bBMpar  = 0.0
#aCMpar = 0.0e-1
##aCMpar = 9.0e-1
#aBMpar = 2.0e-2
#aMTpar = 5.0e-2
#string_case = "osteoblastic"
#string_title = "osteoblastic lesion"

## mixed?
#bBMpar  = 1e1
#aCMpar = 4e0
#aBMpar = 1.5e-1
#aMTpar = 5.0e-2
#string_case = "mixed"
#string_title = "mixed lesion"




aC = p0 
bC = p1
bCT = p2     
aBW = p3
bB = p4
aT = p5
bT = p6    
aW = p7
bW = 1 
kB = p8
kC = p9

def myControl(intON, intOFF, bisPAR, TG1PAR, winPAR, woutPAR, chmPAR):
    # MODEL
    def rhsOFF(x,t):
        # MODEL PARAMETERS
        aC = p0 
        bC = p1
        bCT = p2     
        aBW = p3
        bB = p4
        aT = p5
        bT = p6    
        aW = p7
        bW = 1 
        kB = p8
        kC = p9
        aCM = aCMpar
        aBM = aBMpar
        bBM = bBMpar
        aM = aMpar
        aMT = aMTpar
        KM = KMpar
        
        # CONTROL PARAMETERS
        bis = 0.0
        TG1 = 0.0
        win = 0.0
        wout = 0.0
        chm = 0.0
        
        # STATE VARIABLES        
        xC = x[0]
        xB = x[1]
        xT = x[2]
        xW = x[3]
        BM = x[4]
        xM = x[5]
        
        xCdot = (aC + aCM*xM)*xT/xB - bC*xC - bCT*xT*xC - bis*xC
        xBdot = (aBW + aBM*xM)*xB*xW/xT - bB*xB
        xTdot = aT*kC*xC*BM - bT*xT - TG1*xT
        xWdot = aW*xC - (bW + bBM*xM + wout)*xW + win
        BMdot = kB*xB - kC*xC*BM
        xMdot = xM*(aM + aMT*xT)*(1.0 - xM/KM) - chm*xM
        
        fx = np.array([xCdot, xBdot, xTdot, xWdot, BMdot, xMdot])
        return fx
    
    # MODEL
    def rhsON(x,t):
        # MODEL PARAMETERS
        aC = p0 
        bC = p1
        bCT = p2     
        aBW = p3
        bB = p4
        aT = p5
        bT = p6    
        aW = p7
        bW = 1 
        kB = p8
        kC = p9
        aCM = aCMpar
        aBM = aBMpar
        bBM = bBMpar
        aM = aMpar
        aMT = aMTpar
        KM = KMpar
        
        # CONTROL PARAMETERS
        bis = bisPAR
        TG1 = TG1PAR
        win = winPAR
        wout = woutPAR
        chm = chmPAR
        
        # STATE VARIABLES        
        xC = x[0]
        xB = x[1]
        xT = x[2]
        xW = x[3]
        BM = x[4]
        xM = x[5]
        
        xCdot = (aC + aCM*xM)*xT/xB - bC*xC - bCT*xT*xC - bis*xC
        xBdot = (aBW + aBM*xM)*xB*xW/xT - bB*xB
        xTdot = aT*kC*xC*BM - bT*xT - TG1*xT
        xWdot = aW*xC - (bW + bBM*xM + wout)*xW + win
        BMdot = kB*xB - kC*xC*BM
        xMdot = xM*(aM + aMT*xT)*(1.0 - xM/KM) - chm*xM
        
        fx = np.array([xCdot, xBdot, xTdot, xWdot, BMdot, xMdot])
        return fx
    
    x0 = np.array([xC0, xB0, xT0, xW0, BM0, xM0])
    
    t0Aux = 0
    tfAux = intON
    i = 0
    
    while tfAux < tf:
        timeAux = np.linspace(t0Aux, tfAux, int((tfAux-t0Aux)*N))
        
        if i%2==0:
            mySol = integrate.odeint(rhsON,x0,timeAux)
            t0Aux = tfAux
            tfAux += intOFF
            bisAux = bisPAR*np.ones(len(timeAux))
            TG1Aux = TG1PAR*np.ones(len(timeAux))
            winAux = winPAR*np.ones(len(timeAux))
            woutAux = woutPAR*np.ones(len(timeAux))
            chmAux = chmPAR*np.ones(len(timeAux))
            controlAux = np.array([[bisAux],[TG1Aux],[winAux],[woutAux],[chmAux]])
        if i%2==1:
            mySol = integrate.odeint(rhsOFF,x0,timeAux)
            t0Aux = tfAux
            tfAux += intON
            bisAux = bisPAR*np.zeros(len(timeAux))
            TG1Aux = TG1PAR*np.zeros(len(timeAux))
            winAux = winPAR*np.zeros(len(timeAux))
            woutAux = woutPAR*np.zeros(len(timeAux))
            chmAux = chmPAR*np.zeros(len(timeAux))
            controlAux = np.array([[bisAux],[TG1Aux],[winAux],[woutAux],[chmAux]])
        x0 = np.array([mySol[-1, 0],
                       mySol[-1, 1], 
                       mySol[-1, 2], 
                       mySol[-1, 3], 
                       mySol[-1, 4],
                       mySol[-1, 5]])
        if i==0:
            mySolFinal = mySol
            myTimeFinal = timeAux
            bisFinal = bisAux
            TG1Final = TG1Aux
            winFinal = winAux
            woutFinal = woutAux
            chmFinal = chmAux
            myControlFinal = np.array([[bisFinal],[TG1Final],[winFinal],[woutFinal],[chmFinal]])
        else:
            mySolFinal = np.concatenate((mySolFinal, mySol))
            myTimeFinal = np.concatenate((myTimeFinal, timeAux))
            myControlFinal = np.concatenate((myControlFinal, controlAux),axis=2)
        i += 1
        
    timeAux = np.linspace(t0Aux, tf, int((tf-t0Aux)*N))
    
    if i%2==0:
        mySol = integrate.odeint(rhsON,x0,timeAux)
        t0Aux = tfAux
        tfAux += intOFF
        bisAux = bisPAR*np.ones(len(timeAux))
        TG1Aux = TG1PAR*np.ones(len(timeAux))
        winAux = winPAR*np.ones(len(timeAux))
        woutAux = woutPAR*np.ones(len(timeAux))
        chmAux = chmPAR*np.ones(len(timeAux))
        controlAux = np.array([[bisAux],[TG1Aux],[winAux],[woutAux],[chmAux]])
    if i%2==1:
        mySol = integrate.odeint(rhsOFF,x0,timeAux)
        t0Aux = tfAux
        tfAux += intON
        bisAux = bisPAR*np.zeros(len(timeAux))
        TG1Aux = TG1PAR*np.zeros(len(timeAux))
        winAux = winPAR*np.zeros(len(timeAux))
        woutAux = woutPAR*np.zeros(len(timeAux))
        chmAux = chmPAR*np.zeros(len(timeAux))
        controlAux = np.array([[bisAux],[TG1Aux],[winAux],[woutAux],[chmAux]])
        
    x0 = np.array([mySol[-1, 0],
                   mySol[-1, 1], 
                   mySol[-1, 2], 
                   mySol[-1, 3], 
                   mySol[-1, 4],
                   mySol[-1, 5]])
    
    mySolFinal = np.concatenate((mySolFinal, mySol))
    myTimeFinal = np.concatenate((myTimeFinal, timeAux))
    myControlFinal = np.concatenate((myControlFinal, controlAux),axis=2)
    
    return mySolFinal, myControlFinal, myTimeFinal

def extract_ctrl(folder):
    x  = []
    bis = []
    tgf = []
    wnt = []
    chm = []
    
    folderName = './'+folder+'/'
    
    with open(folderName + 'discretization_times.export','r') as csvfile:
        plots = csv.reader(csvfile)
        for row in plots:
            x.append(float(row[0]))
    
    with open(folderName + 'bis.export','r') as csvfile:
        plots = csv.reader(csvfile)
        for row in plots:
            bis.append(float(row[0]))
    
    with open(folderName + 'antitgf.export','r') as csvfile:
        plots = csv.reader(csvfile)
        for row in plots:
            tgf.append(float(row[0]))
            
    with open(folderName + 'wnt.export','r') as csvfile:
        plots = csv.reader(csvfile)
        for row in plots:
            wnt.append(float(row[0]))
            
    with open(folderName + 'chemo.export','r') as csvfile:
        plots = csv.reader(csvfile)
        for row in plots:
            chm.append(float(row[0]))
            
            
    xC = []
    
    with open(folderName + 'xC.export','r') as csvfile:
        plots = csv.reader(csvfile)
        for row in plots:
            xC.append(float(row[0]))
    
    xB = []
    
    with open(folderName + 'xB.export','r') as csvfile:
        plots = csv.reader(csvfile)
        for row in plots:
            xB.append(float(row[0]))
    
    xT = []
    
    with open(folderName + 'xT.export','r') as csvfile:
        plots = csv.reader(csvfile)
        for row in plots:
            xT.append(float(row[0]))
    
    xW = []
    
    with open(folderName + 'xW.export','r') as csvfile:
        plots = csv.reader(csvfile)
        for row in plots:
            xW.append(float(row[0]))
            
    BM = []
    
    with open(folderName + 'BM.export','r') as csvfile:
        plots = csv.reader(csvfile)
        for row in plots:
            BM.append(float(row[0]))
            
    xM = []
    
    with open(folderName + 'xM.export','r') as csvfile:
        plots = csv.reader(csvfile)
        for row in plots:
            xM.append(float(row[0]))


            
    return x, bis, tgf, wnt, chm, xC, xB, xT, xW, BM, xM



    
#-------------------#
#
# control(TimePeriod, alpha, theta, d2, m)
#
#-------------------#


def plotResults(state,control,time):
    
    # no control
    NOstate, NOcontrol, NOtime = myControl(1, 1, 0, 0, 0, 0, 0)
    
    
    plt.figure(figsize=(8,3.5))
    
    plt.subplot(2,4,1)
    plt.plot(NOtime,NOstate[:,0],':',color='#1f77b4',linewidth=1)
    plt.plot(time,state[:,0],'#1f77b4',linewidth=1)
    plt.ylabel('OCs',fontsize=12)
    plt.xticks([0, 250])
    plt.xlabel('Time (days)',fontsize=12)
    #plt.ylim([2,7])
    #plt.yticks([2,4.5,7])
    plt.subplots_adjust(wspace=0, hspace=0)
    
    plt.subplot(2,4,2)
    plt.plot(NOtime,NOstate[:,1],':',color='#ff7f0e',linewidth=1)
    plt.plot(time,state[:,1],'#ff7f0e',linewidth=1)
    plt.ylabel('OBs',fontsize=12)
    plt.xticks([0, 250])
    plt.xlabel('Time (days)',fontsize=12)
    #plt.ylim([0,9])
    #plt.yticks([0,4.5,9])
    plt.subplots_adjust(wspace=0, hspace=0)
    
    plt.subplot(2,4,3)
    plt.plot(NOtime,NOstate[:,5],':',color='#ba34eb',linewidth=1)
    plt.plot(time,state[:,5],'#ba34eb',linewidth=1)
    plt.ylabel('CCs',fontsize=12)
    plt.xticks([0, 250])
    plt.xlabel('Time (days)',fontsize=12)
    #plt.ylim([0,9])
    #plt.yticks([0,4.5,9])
    plt.subplots_adjust(wspace=0, hspace=0)
    
    plt.subplot(2,4,4)
    plt.plot(NOtime,NOstate[:,4],':',color='b',alpha=0.5,linewidth=1)
    plt.plot(time,state[:,4],'#000000',linewidth=1)
    plt.xticks([0, 250])
    plt.xlabel('Time (days)',fontsize=12)
    #plt.xlim([2,7])
    #plt.xticks([2,4.5,7])
    plt.ylim([0,2])
    #plt.yticks([0,1,4])
    plt.ylabel('Bone Mass',fontsize=12)
    plt.subplots_adjust(wspace=0, hspace=0)
    
    plt.subplot(2,4,5)
    plt.plot(time,control[0,0,:],'k',linewidth=1)
    #plt.ylim([0,0.5])
    #plt.yticks([0.0, 0.25, 0.5])
    #plt.ylim([0,0.01])
    #plt.yticks([0.0, 0.005, 0.01])
    plt.fill_between(time,control[0,0,:],color='k',alpha=0.5)
    plt.xlabel('Time (days)',fontsize=12)
    plt.ylabel('Bisphosphonate',fontsize=12)
    plt.xticks([0, 250])
    #plt.ylim([-0.05, 0.5])
    plt.yticks([np.min(control[0,0,:]), np.max(control[0,0,:])], ['off','on'])
    plt.xlabel('Time (days)',fontsize=12)
    plt.subplots_adjust(wspace=0, hspace=0)
    
    plt.subplot(2,4,6)
    plt.plot(time,control[1,0,:],'c',linewidth=1)
    #plt.ylim([0,100])
    #plt.yticks([0, 50, 100])
    plt.fill_between(time,control[1,0,:],color='c',alpha=0.5)
    plt.xlabel('Time (days)',fontsize=12)
    plt.ylabel('TGF inhibition',fontsize=12)
    plt.xticks([0, 250])
    plt.yticks([np.min(control[1,0,:]), np.max(control[1,0,:])], ['off','on'])
    #plt.ylim([-100, 1000])
    plt.xlabel('Time (days)',fontsize=12)
    plt.subplots_adjust(wspace=0, hspace=0)
    
#    plt.subplot(2,4,7)
#    plt.plot(time,control[2,0,:],'m',linewidth=1)
#    #plt.ylim([0,2])
#    #plt.yticks([0, 1, 2])
#    plt.fill_between(time,control[2,0,:],color='m',alpha=0.5)
#    plt.xlabel('Time (days)',fontsize=12)
#    plt.ylabel('Wnt input',fontsize=12)
#    plt.yticks([np.min(control[2,0,:]), np.max(control[2,0,:])], ['off','on'])
#    plt.xticks([0, 250])
#    #plt.ylim([-0.1,1])
#    plt.xlabel('Time (days)',fontsize=12)
    
    plt.subplot(2,4,7)
    plt.plot(time,control[3,0,:],color='#b3d345',linewidth=1)
    #plt.ylim([0,2])
    #plt.yticks([0, 1, 2])
    plt.fill_between(time,control[3,0,:],color='#b3d345',alpha=0.5)
    plt.xlabel('Time (days)',fontsize=12)
    plt.ylabel('Wnt inhibition',fontsize=12)
    plt.yticks([np.min(control[3,0,:]), np.max(control[3,0,:])], ['off','on'])
    plt.xticks([0, 250])
    #plt.ylim([-0.1,1])
    plt.xlabel('Time (days)',fontsize=12)
    
    plt.subplot(2,4,8)
    plt.plot(time,control[4,0,:],'r',linewidth=1)
    #plt.ylim([0,2])
    #plt.yticks([0, 1, 2])
    plt.fill_between(time,control[4,0,:],color='r',alpha=0.5)
    plt.xlabel('Time (days)',fontsize=12)
    plt.ylabel('Chemotherapy',fontsize=12)
    plt.yticks([np.min(control[4,0,:]), np.max(control[4,0,:])], ['off','on'])
    plt.xticks([0, 250])
    #plt.ylim([-0.1,1])
    plt.xlabel('Time (days)',fontsize=12)
    
    plt.subplots_adjust(wspace=0, hspace=0)
    plt.tight_layout()

plt.close('all')


# (intON, intOFF, bis, tg1, win, wout, chm)
state, control, time = myControl(249, 1, 0.2, 1, 0, 0, 0.05)
plotResults(state, control, time)
plt.savefig("control-met-1-inter-1.pdf")

state, control, time = myControl(2, 12, 0.2, 1, 0, 0, 0.05)
plotResults(state, control, time)
plt.savefig("control-met-1-inter-2.pdf")

x, bis, tgf, wnt, chm, xC, xB, xT, xW, BM, xM = extract_ctrl('met-1-sol')
optimalState = np.transpose(np.array([[xC],[xB],[xT],[xW],[BM],[xM]]))
optimalControl = np.array([[bis],[tgf],[wnt],[wnt],[chm]])
optimalTime = np.array(x)
plotResults(optimalState[:,0,:], optimalControl, optimalTime)
plt.savefig("control-met-1-optimal-1.pdf")

# tg1
fig = plt.figure(figsize=(4,2))
ax = plt.subplot(111)
my_soln, my_control, time = myControl(2, 2, 0.0, 0.0, 0.0, 0.0, 0.0)
plt.plot(time, my_soln[:,4],color='black', linestyle='dashed', label='no control')#
my_soln, my_control, time = myControl(2, 12, 0.0, 1.0, 0.0, 0.0, 0.0)
plt.plot(time, my_soln[:,4], label=r'$u_{tg1}=1$')
my_soln, my_control, time = myControl(2, 12, 0.0, 10, 0.0, 0.0, 0.0)
plt.plot(time, my_soln[:,4], label=r'$u_{tg1}=10$')
my_soln, my_control, time = myControl(2, 12, 0.0, 100, 0.0, 0.0, 0.0)
plt.plot(time, my_soln[:,4], label=r'$u_{tg1}=100$')
plt.ylim([0.5,2.5])
plt.xlabel('Time (days)',fontsize=12)
plt.ylabel('Bone Mass',fontsize=12)
plt.legend()
plt.tight_layout()
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.65, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig("control-met-1-tg1.pdf")

# chm
fig = plt.figure(figsize=(4,2))
ax = plt.subplot(111)
my_soln, my_control, time = myControl(2, 2, 0.0, 0.0, 0.0, 0.0, 0.0)
plt.plot(time, my_soln[:,4],color='black', linestyle='dashed', label='no control')#
my_soln, my_control, time = myControl(2, 12, 0.0, 0.0, 0.0, 0.0, 0.1)
plt.plot(time, my_soln[:,4], label=r'$u_{chm}=0.1$')
my_soln, my_control, time = myControl(2, 12, 0.0, 0, 0.0, 0.0, 0.2)
plt.plot(time, my_soln[:,4], label=r'$u_{chm}=0.2$')
my_soln, my_control, time = myControl(2, 12, 0.0, 0, 0.0, 0.0, 0.5)
plt.plot(time, my_soln[:,4], label=r'$u_{chm}=0.5$')
plt.ylim([0.5,2.5])
plt.xlabel('Time (days)',fontsize=12)
plt.ylabel('Bone Mass',fontsize=12)
plt.legend()
plt.tight_layout()
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.65, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig("control-met-1-chm.pdf")

# bis
fig = plt.figure(figsize=(4,2))
ax = plt.subplot(111)
my_soln, my_control, time = myControl(2, 2, 0.0, 0.0, 0.0, 0.0, 0.0)
plt.plot(time, my_soln[:,4],color='black', linestyle='dashed', label='no control')#
my_soln, my_control, time = myControl(2, 12, 0.1, 0.0, 0.0, 0.0, 0.0)
plt.plot(time, my_soln[:,4], label=r'$u_{bis}=0.1$')
my_soln, my_control, time = myControl(2, 12, 0.2, 0, 0.0, 0.0, 0.0)
plt.plot(time, my_soln[:,4], label=r'$u_{bis}=0.2$')
my_soln, my_control, time = myControl(2, 12, 0.5, 0, 0.0, 0.0, 0.0)
plt.plot(time, my_soln[:,4], label=r'$u_{bis}=0.5$')
plt.ylim([0.5,2.5])
plt.xlabel('Time (days)',fontsize=12)
plt.ylabel('Bone Mass',fontsize=12)
plt.legend()
plt.tight_layout()
box = ax.get_position()
ax.set_position([box.x0, box.y0, box.width * 0.65, box.height])
ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
plt.savefig("control-met-1-bis.pdf")

