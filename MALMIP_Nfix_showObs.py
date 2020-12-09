# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 16:56:48 2020
@author: Benjamin
"""

import os

#Geophysical libraries
import pybert as pb
from pybert import tdip
import pygimli as pg

# My own library
from utils_rhizo import fct_utils as FU
from MALM4RootsRhizo_class import MALM4RootsRhizo_class as MR

# General libraries
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors

# Visualisation libraries
import pyvista as pv

# MALM Inversion library
from icsd3d_class import iCSD3d as i3d
from plotters.mpl_plot import showObs2d
from importers.read import *

plt.close('all')

#%% Define survey parameters (see excel file for correspondance dates/files/parameters)

waterRes = 21.23 # water resistivity (Ohm.m)
rec = True # reciprocal analysis
A = 72-1 # Id electrode A
B = 65-1 # Id electrode B
injection_duration = 1 # time of injection
date = '0112' # 0112 '1310' '1611' # (ddmm)
# 'MALMIP1201.bin' 'MALMIP1013' 'MALMIP1116.bin' # filename in raw data folder
inputfile = 'MALMIP1201.bin' #  
plotCC = False # show Cole-Cole fitted parameters (to compare with literature)
split_Nfix = [True, 71-1]
Nfix = 71-1 #71-1  #None 71-1 # Id electrode N , put None if N is varying

rmvInvalid = True # True if you want to write filtered files/ False for raw data
rmv_outliers = True
rmv_id= 61-1 #None # 61-1 o None

all_gates= True
if not all_gates: 
    gateIP = 3 
else: gateIP = None # id extracted time window in mV/V (?)

icsd = False
#%% define PATHS 
# or set you working directory local path here 
main = os.getcwd()
os.chdir(main)

geomPath= './geom/'
meshPath= './mesh/'
icsdPath= './icsd/'

icsdPath += date 
if rmvInvalid is False: 
    icsdPath += '_raw' 
if all_gates: 
    icsdPath += '_AllG' 
if Nfix: 
    icsdPath += '_Nfix' 
if gateIP: 
    icsdPath += '_M' + str(gateIP)
icsdPath += '/'
    
#%% LOAD geometry and mesh
RemLineNb, Injection, coordE, pointsE= load_geom(geomPath) # geometry file containing electrodes position includinf remotes 
mesh3d, sensors = FU.mesh_import(meshPath + 'BaseRhizo_Vrte.msh')

#%% 

def filterTDIP(dataTDIP,id2rmv):
    valid = np.ones(len(dataTDIP.data('m')))
    valid[id2rmv] = 0
    dataTDIP.data.set('valid',valid)
    dataTDIP.MA = dataTDIP.MA[:, dataTDIP.data['valid'].array()==1]
    dataTDIP.data.removeInvalid()

    return dataTDIP, valid

#%% Import data
IPcurves = tdip.TDIPdata('./raw_data/' + inputfile) # e.g. ABEM or Syscal TXT export
valid = np.ones(len(IPcurves.data('m')))

#%% split and filter
IPcurves_f = tdip.TDIPdata('./raw_data/' + inputfile)

if split_Nfix[0]:

    IPcurves.data('m').array()
    idn = np.where(IPcurves.data('m')==split_Nfix[1])[0]
    idm = np.where(IPcurves.data('n')==split_Nfix[1])[0]
    idfix = list(idn) + list(idm)
    
    IPcurves.data('a')[idfix].array()
    IPcurves.data('b')[idfix].array()
    id_elec_2rmv = idfix # remove Nfix electrodes
    
    if Nfix is not None:
        a = set(list(range(0, len(IPcurves.data('m')))))
        id_elec_2rmv = a.difference(set(idfix))
        id_elec_2rmv = list(id_elec_2rmv)
    
    IPcurves_f, valid_split = filterTDIP(IPcurves_f,id_elec_2rmv)

#%% Import data and visualisation

if rmv_outliers: 
    id_outliers = np.where(abs(IPcurves_f.data['M1'])>100)[0]
    id_elec_2rmv= list(id_outliers)
    if rmv_id: 
        id_elec_2rmv.append(rmv_id)
    IPcurves_f, valid_outliers = filterTDIP(IPcurves_f,id_elec_2rmv)

#%%
j=0
# Ensure variable is defined
try:
    valid_split
except NameError:
    valid_split = np.ones(len(IPcurves.data('m')))
    
for i, v in enumerate(valid_split):
    if v==1: # if valid in split
        # check if valid in outliers
        if rmv_outliers:
            if valid_outliers[j] == 0:
                valid[i] = 0
                print(i)
            else:
                valid[i] = 1
            j = j + 1
    else: 
        valid[i] = 0    
        
np.count_nonzero(valid==0)



#%%

# IPcurves.filter(m=11, n=12)
# IPcurves.filter(m=46, n=53)
# IPcurves.filter(m=53, n=44)
# IPcurves.filter(m=48, n=56)
# IPcurves.filter(m=51, n=58)
# IPcurves.filter(m=55, n=31)
# IPcurves.filter(m=60, n=61)
# IPcurves.filter(m=69, n=60)
# IPcurves.filter(m=54, n=45)
# IPcurves.filter(m=59, n=51)
# IPcurves.filter(m=45, n=52)
# IPcurves.filter(m=68, n=59)

# idmn_2_rmv = np.array(  [[11,12],[46,53],[53,44],[48,56],
#                         [51,58],[55,31],[60,61],[69,60],
#                         [54,45],[59,51],[45,52],[68,59]]
#                         )

# for i,j in enumerate(idmn_2_rmv):
#     if (IPcurves_f.data['m']==j[0]).any():
#         if (IPcurves_f.data['n']==j[1]).any(): 
#             print('add')
#             bool2rmv[i]= 1
#         else: 
#             bool2rmv[i]= 0

# IPcurves_f.data['valid'] = bool2rmv
# IPcurves_f.MA = IPcurves_f.MA[:, IPcurves_f.data['valid'].array()==1]
# IPcurves_f.data.removeInvalid()

#id_elec_2rmv = [42,44,61] # remove some noisy electrodes
 #  = [] # remove some noisy electrodes
#for i in id_elec_2rmv:
#   bool2rmv = IPcurves_f.data['m']==i
#   if bool2rmv.any() == False: 
#       print('search in n')
#       bool2rmv = IPcurves_f.data['n']==i
#   IPcurves_f.data.markInvalid(bool2rmv)
#  IPcurves.data.markInvalid(bool2rmv)
#IPcurves_f.MA = IPcurves_f.MA[:, IPcurves_f.data['valid'].array()==1]
#IPcurves_f.data.removeInvalid()


# bool_outliers = np.abs(IPcurves.MA[0])>10
# IPcurves_f.data.markInvalid(bool_outliers)
# IPcurves_f.MA = IPcurves_f.MA[:, IPcurves_f.data['valid'].array()==1]
# IPcurves_f.data.removeInvalid()

#id_elec_2rmv = [42,44,61] # remove some noisy electrodes
#for i in id_elec_2rmv:
#    bool2rmv = IPcurves_f.data['m']==i
#    if bool2rmv.any() == False: 
#        print('search in n')
#        bool2rmv = IPcurves_f.data['n']==i
#    IPcurves_f.data.markInvalid(bool2rmv)
#    IPcurves.data.markInvalid(bool2rmv)
#IPcurves_f.MA = IPcurves_f.MA[:, IPcurves_f.data['valid'].array()==1]
#IPcurves_f.data.removeInvalid()

IPcurves_f.showDecay(nr=np.arange(0,len(IPcurves_f.data['a'])), showFit=False, 
                   yscale='linear',xscale='linear')


if rmvInvalid is False:
    IPcurves_f = tdip.TDIPdata('./raw_data/' + inputfile)
    
#%% plot Cole Cole model parameters 
if plotCC:

    #m0,tau,fit = IPcurves_f.fitDecays(show=True)
    #IPcurves_f.showDecay(nr=np.arange(0,len(IPcurves_f.data['a'])), showFit=False, 
    #                   yscale='linear',xscale='linear')
    #m0_02,tau_02,fit_02 = IPcurves_f.fitDecays(show=True, tmin=0.2)
    #IPcurves_f.showDecay(nr=np.arange(0,len(IPcurves_f.data['a'])), showFit=False, 
    #                   yscale='linear',xscale='linear')

    fig, ax = plt.subplots(nrows=1, ncols=4)
    ax[0].violinplot(m0)
    ax[1].violinplot(tau)
    ax[2].violinplot(fit)
    ax[3].violinplot(abs(IPcurves.data('r')))
    # ax.set_title('Default violin plot')
    ax[0].set_ylabel('Observed values')
    ax[0].set_title('m0 (mV/V)')
    ax[1].set_title('Tau (s)')
    ax[2].set_title('fit rms (log)')
    ax[3].set_title('rhoa')
    
    
    fig, ax = plt.subplots(nrows=1, ncols=4)
    ax[0].violinplot(m0_02)
    ax[1].violinplot(tau_02)
    ax[2].violinplot(fit_02)
    ax[3].violinplot(abs(IPcurves.data('r')))
    # ax.set_title('Default violin plot')
    ax[0].set_ylabel('Observed values')
    ax[0].set_title('m0 (mV/V)')
    ax[1].set_title('Tau (s)')
    ax[2].set_title('fit rms (log)')
    ax[3].set_title('rhoa')


#%% select observation data - in preparaion for ICSD analysis --> 3 files required (same folder)
# 1st file needed = virtual electrode position
VRTEpos= FU.VRTEpos(mesh=mesh3d,dim=3) # return VRTE file format icsd MarkerVRTE=991
#plt.scatter(VRTEpos[:,0],VRTEpos[:,1],color='b')

# 2nd file = observations data
Obs_raw = pb.importer.importSyscalPro('./raw_data/' + inputfile) 
Obs, dataABMN=  FU.PrepareMALMData('./raw_data/' + inputfile, Rec=False, DevErr=1,
                           MinV=1, MaxRc=1, Kfact=1, MinMaxAppRes=1, 
                           SwE=False, 
                           valid=valid)

# 2nd file = observations data = TDIP
Obs_IP, dataABMN=  FU.PrepareMALMData('./raw_data/' + inputfile, Rec=False, DevErr=1,
                           MinV=1, MaxRc=1, Kfact=1, MinMaxAppRes=1, 
                           SwE=False, gIP=gateIP, 
                           valid=valid) # remove reciprocal from original dataset

if all_gates:
    for i, g in enumerate(range(1,20,1)):
        Obs_IP, dataABMN=  FU.PrepareMALMData('./raw_data/' + inputfile, Rec=False, DevErr=1,
                                   MinV=1, MaxRc=1, Kfact=1, MinMaxAppRes=1, 
                                   SwE=False, gIP=g, 
                                   valid=valid) # remove reciprocal from original dataset

#Obs.setSensorPositions(coordE[:,1:3])
#Obs.save('obs.data')

# filter out wrong electrodes 

coordE_f = []
for i, mi in enumerate(Obs['m']):
    if mi==Nfix:
       mi=Obs['n'][i]
    id_coordE_f = np.where(mi+1==coordE[:,0])[0]
    #if len(id_coordE_f) > 1:
    coordE_f.append(coordE[id_coordE_f[0],:])
coordE_f = np.array(coordE_f)

# filter out wrong observation 
if (rmvInvalid or split_Nfix[0]) is False:
    print('raw analysis')
    Obs, dataABMN=  FU.PrepareMALMData('./raw_data/' + inputfile, Rec=False, DevErr=1,
                               MinV=1, MaxRc=1, Kfact=1, MinMaxAppRes=1, 
                               SwE=False,
                               valid=None)
    Obs_IP, dataABMN=  FU.PrepareMALMData('./raw_data/' + inputfile, Rec=False, DevErr=1,
                               MinV=1, MaxRc=1, Kfact=1, MinMaxAppRes=1, 
                               SwE=False, gIP=gateIP, 
                               valid=None) # remove reciprocal from original dataset
#len(dataABMN)

    coordE_f = []
    for i, mi in enumerate(Obs['m']):
        if mi==Nfix:
           mi=Obs['n'][i]
        id_coordE_f = np.where(mi+1==coordE[:,0])[0]
        coordE_f.append(coordE[id_coordE_f[0],:])
    coordE_f = np.array(coordE_f)
    


#%% check consistency

IPcurves_f
Obs['M'+str(gIP)]
abs(Obs('r').array())[60]
Obs('M1').array()

#%% foward modelling of Green's functions
# mk_full_malm --> 3rd file = simulated green functions (simulated data to compare with observation data)
SeqFullMALM= MR.mk_full_malm(dataABMN-1, VRTe = range(len(sensors), len(sensors)+len(VRTEpos)),
                             mesh=mesh3d, R3=False) # output .shm with sensors

rhomap = [[1, waterRes],
          [2, waterRes]]

rho = pg.solver.parseArgToArray(rhomap, mesh3d.cellCount(), mesh3d)
MR.SimulateGreenFcts(mesh_VRTs=mesh3d,rhomapVRTs=rho,schemeVrts=SeqFullMALM, 
                     Name='VRTeSim')

#%% copy file to icsd folder
import glob
import shutil

for fpath in glob.glob('OMALM*.txt'):  # this returns a list of the CURRENT contents. Therefore, there is no need to sanitize for the directories that we will later create
    if not os.path.isdir(icsdPath):
        os.mkdir(icsdPath)
    #shutil.move(fpath, icsdPath)
    shutil.move(os.path.join(os.getcwd(), fpath), os.path.join(icsdPath, fpath))

shutil.copy('VRTeSim.txt', icsdPath)
shutil.copy('VRTeCoord.txt', icsdPath)


#%% Quiver plot (gradient*conductivity)
if Nfix is not None:

    FU.streamlines(coordE_f, abs(Obs('r').array()), waterRes,
                   sensors=sensors, A=A, B=B, Nfix=Nfix)
    
    fig, ax = plt.subplots(nrows=1, ncols=4)
    for i, g in enumerate(range(1,20,5)):
        FU.streamlines(coordE_f, Obs('M'+str(g)).array(), waterRes,
                       sensors=sensors, A=A, B=B, Nfix=Nfix, ax=ax[i],
                       vmin=-100, vmax=10)
        ax[i].set_title('Gate t:' + str(IPcurves.t[g]) + 's')
        
    fig, ax = plt.subplots(nrows=1, ncols=4)
    for i, g in enumerate(range(1,20,5)):
        FU.streamlines(coordE_f, Obs('M'+str(g)).array(), waterRes,
                       sensors=sensors, A=A, B=B, Nfix=Nfix, ax=ax[i]
                       )
        ax[i].set_title('Gate t:' + str(IPcurves.t[g]) + 's')

        
#%%
if icsd:
    from icsd3d_class import iCSD3d as i3d
    from plotters import mpl_plot
    
    #mpl_plot.showObs2d(path2files+'/')
    sep = '_'
    NameSave = 'O'+ os.path.basename(inputfile).split(sep, 1)[0] + '.txt'
    
    path2files= icsdPath
    icsd=i3d(dirName=path2files +'/')   
    icsd.createSurvey(fname_obs=NameSave,fname_sim='VRTeSim.txt')
    icsd.regMesh='strc'
    icsd.type='2d'
    #icsd.mesh='invdir/fwd/forward_model.vtk'
    icsd.obs_err='const' # sqrt choose between constant weight and w = 1/sqrt(abs(obs))
    icsd.wr=1 #weight regularization
    icsd.alphaSxy=False
    icsd.x0_prior=False
    icsd.x0_ini_guess=False # initial guess
    # icsd.plotElecs=False
    icsd.method_m0='F1'
    m0 = icsd.estimateM0(method_m0='F1',show=True)
    
    # icsd.clim=[0,0.1]
    # icsd.run_single()
    icsd.alphax0=1
    
    sol= icsd.invert(wr=1)
    
    fig, ax = plt.subplots()
    icsd.showResults(ax=ax)
    plt.show()
    
    
    if all_gates:
        
        NameSave = 'O'+ os.path.basename(inputfile).split(sep, 1)[0] + '.txt'
        fig, ax = plt.subplots()
        icsd=i3d(dirName=path2files +'/')   
        icsd.createSurvey(fname_obs=NameSave,fname_sim='VRTeSim.txt')
        icsd.method_m0='F1'
        m0 = icsd.estimateM0(method_m0='F1',show=True,ax=ax)
        
        fig, ax = plt.subplots(nrows=1, ncols=3)
        for i, g in enumerate(range(1,20,8)):
            NameSave = 'O'+ os.path.basename(inputfile).split(sep, 1)[0] + 'M'+ str(g) + '.txt'
            print(NameSave)
            icsd=i3d(dirName=path2files +'/')   
            icsd.createSurvey(fname_obs=NameSave,fname_sim='VRTeSim.txt')
            m0 = icsd.estimateM0(method_m0='F1',show=True, ax=ax[i])
            ax[i].set_title('Gate t:' + str(IPcurves.t[g]) + 's')
            plt.tight_layout()
    
    
        fig, ax = plt.subplots(nrows=1, ncols=4)
        for i, g in enumerate(range(1,20,5)):
            NameSave = 'O'+ os.path.basename(inputfile).split(sep, 1)[0] + 'M'+ str(g) + '.txt'
            icsd=i3d(dirName=path2files +'/')   
            icsd.createSurvey(fname_obs=NameSave,fname_sim='VRTeSim.txt')
            sol= icsd.invert(wr=1)
            icsd.showResults(ax=ax[i])
            plt.show()
            ax[i].set_title('Gate t:' + str(IPcurves.t[g]) + 's')
    
    #    fig, ax = plt.subplots(nrows=1, ncols=4)
    #    icsd.createTimeLapseSurvey(fname_obs=NameSave+'M'+str(g),fname_sim='VRTeSim.txt')
    #    for i, g in enumerate(range(1,20,5)):
    #        icsd=i3d(dirName=path2files +'/')   
    #        icsd.icsd_init(icsd)
    #        sol= icsd.invert(wr=1)
    #        icsd.showResults(ax=ax[i])
    #        plt.show()
    #        ax[i].set_title('Gate t:' + str(IPcurves.t[g]) + 's')

