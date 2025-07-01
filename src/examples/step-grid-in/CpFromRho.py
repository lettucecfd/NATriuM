import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib import rc
rc('font', **{'size'   : 9, 'family': 'serif', 'serif': ['Charis SIL']})
rc('text', usetex=False)

basepath = "/mnt/c/Users/phili/Desktop/"
jobids = ["9682930", "dw872713"]
vtupaths = ["9682930/output", "NACA/dw872713-natrium/step-grid-in-old/Re10000-Ma1.5-reflevel0-time1695800330"]

imgpath = basepath + "images/"
if not os.path.exists(imgpath):
    os.mkdir(imgpath)

cs = 0.563436
rho0LU = 1
gamma = 1.4
Ma = 1.5#*np.sqrt(1.4)
dx = 0.000508780
dt = 9.21616e-05
L = 1
a = 343
Tref = 1
RLU = 1

ref = np.loadtxt('/home/philipp/NATriuM/NATriuM/src/examples/step-grid-in/refs/reference.txt',delimiter=';',skiprows=1)
frap = np.loadtxt('/home/philipp/NATriuM/NATriuM/src/examples/step-grid-in/refs/frapolli.txt',delimiter=';',skiprows=1)

    
for jobid, vtupath in zip(jobids, vtupaths):
    results = np.loadtxt(imgpath + "rhoUxUy_" + jobid + ".txt")

    rhoLU = results[:,2]
    pLU = rhoLU*cs*cs
    p0LU = rho0LU*cs*cs
    uxLU = results[:,3]
    uyLU = results[:,4]
    TLU = results[:,5]
    
    UmagLU = np.sqrt(uxLU*uxLU+uyLU*uyLU)
    MaLocal = UmagLU*1.5*1.5*1.5*gamma*Tref/3

    fig, ax = plt.subplots()
    ax.plot(results[:,0], MaLocal)
    fig.savefig(imgpath + "MaLocal_" + jobid + ".png")

    # V1
    CpLU = (pLU-p0LU)/(0.5*gamma*p0LU*Ma*Ma)

    # V2
    pPU = pLU / ((dx*dx)/(L*L)) * (cs*cs)/(a*a) * rhoLU/rho0LU
    p0PU = p0LU / ((dx*dx)/(L*L)) * (cs*cs)/(a*a) * rho0LU/rho0LU
    CpPU = (pPU-p0PU)/(0.5*gamma*p0PU*Ma*Ma)

    # V3
    pPU = pLU * (dx*dx)/(dt*dt) / rho0LU
    p0PU = p0LU * (dx*dx)/(dt*dt) / rho0LU
    CpPU = (pPU-p0PU)/(0.5*gamma*p0PU*Ma*Ma)

    # V4
    CrhoLU = (rhoLU-rho0LU)/(0.5*gamma*rho0LU*Ma*Ma)

    # V5
    CpLocalLU = (rhoLU-rho0LU)/(0.5*gamma*rho0LU*MaLocal*MaLocal)

    # V6
    pLU = rhoLU*RLU*TLU
    p0LU = rho0LU*RLU*Tref
    CpTLU = (pLU-p0LU)/(0.5*gamma*p0LU*Ma*Ma)

    Cp = CpPU
    results = np.hstack((results,np.expand_dims(Cp,1)))
    np.savetxt(imgpath + "Cp_" + jobid + ".txt", results)

    fig, ax = plt.subplots()
    ax.plot(ref[:,0],ref[:,1],'o',label='Latt et al.',color='royalblue')
    ax.plot(frap[:,0],frap[:,1],'*',label='Frapolli et al.',color='red')
    ax.set_ylim(1.7,-0.3)
    ax.set_xlim(-1,1.5)
    ax.set_xlabel('x/C')
    ax.set_ylabel(r'$C_P$')
    # ax.plot(results[:,0], CpLU,     label='SLLBM LU')
    # ax.plot(results[:,0], CpPU,     label='SLLBM PU')
    # ax.plot(results[:,0], CrhoLU,   label='SLLBM CpRhoLU')
    # ax.plot(results[:,0], CpLocalLU,label='SLLBM Cp from MaLocal')
    # ax.plot(results[:,0], CpTLU,    label='SLLBM CpTLU')
    ax.plot(results[:,0], CpTLU,    color='black',    label='SLLBM')
    ax.legend(frameon=False, loc='lower right')
    fig.savefig(imgpath + "Cp_" + jobid + ".png")

    print(f"Finished jobid {jobid}")

