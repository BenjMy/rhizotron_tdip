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
from run_ICSD import icsd_TDIP_plant
from run_invert_pygimli_ERT import invert_Resipy_ERT, invert_pygimli_ERT
# General libraries
import numpy as np
import matplotlib.pyplot as plt

# Visualisation libraries
import pyvista as pv

# MALM Inversion library
from icsd3d_class import iCSD3d as i3d
from plotters.mpl_plot import showObs2d
from importers.read import *

plt.close('all')

#%% Define survey parameters (see excel file for correspondance dates/files/parameters)
invERT = False
waterRes = (1/767)*1e4 #21.23 # water resistivity (Ohm.m) converted from microS/cm 767, 477 
rec = True # reciprocal analysis
A = 72-1 # Id electrode A
B = 65-1 # Id electrode B
injection_duration = 2 # time of injection
date = '0113' #  '1712' 0112 '1310' '1611' # (ddmm)
# 'MALMIP1217' 'MALMIP1201.bin' 'MALMIP1013' 'MALMIP1116.bin' # filename in raw data folder
inputfile = 'MALMIP_0113.bin' #  
plotCC = False # show Cole-Cole fitted parameters (to compare with literature)
split_Nfix = [True, 71-1]
Nfix = 71-1 #71-1  #None 71-1 # Id electrode N , put None if N is varying

rmvInvalid = False # True if you want to write filtered files/ False for raw data
rmv_outliers = False
rmv_id= None #None # 61-1 o None

all_gates= True
if not all_gates: 
    gateIP = 3 
else: gateIP = None # id extracted time window in mV/V (?)

icsd = True

main = os.getcwd()
os.chdir(main)
geomPath, meshPath, icsdPath = FU.definePath(main,date)

    
#%% LOAD geometry and mesh
RemLineNb, Injection, coordE, pointsE= load_geom(geomPath) # geometry file containing electrodes position includinf remotes 
mesh3d, sensors = FU.mesh_import(meshPath + 'BaseRhizo_Vrte.msh')

#%% INVERT ERT data

if invERT:
    # inputfileERT = 'ERT_0113.csv'
    # k = invert_Resipy_ERT(inputfileERT)

    inputfileERT = 'ERT_0113.bin'
    model = invert_pygimli_ERT(inputfileERT,sensors,mesh3d)
    
#%% filterTDIP

def filterTDIP(dataTDIP,id2rmv):
    valid = np.ones(len(dataTDIP.data('m')))
    valid[id2rmv] = 0
    dataTDIP.data.set('valid',valid)
    dataTDIP.MA = dataTDIP.MA[:, dataTDIP.data['valid'].array()==1]
    dataTDIP.data.removeInvalid()

    return dataTDIP, valid

#%% Import data TDIP
IPcurves = tdip.TDIPdata('./raw_data/' + inputfile) # e.g. ABEM or Syscal TXT export
valid = np.ones(len(IPcurves.data('m')))

#%% split and filter
IPcurves_f = tdip.TDIPdata('./raw_data/' + inputfile)

if split_Nfix[0]:

    IPcurves.data('m').array()
    idn = np.where(IPcurves.data('m')==split_Nfix[1])[0]
    idm = np.where(IPcurves.data('n')==split_Nfix[1])[0]
    idfix = list(idn)  #+ list(idm)
    
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

IPcurves_f.showDecay(nr=np.arange(0,len(IPcurves_f.data['a'])), showFit=False, 
                   yscale='linear',xscale='linear')


if rmvInvalid is False:
    IPcurves_f = tdip.TDIPdata('./raw_data/' + inputfile)
    
#%% plot Cole Cole model parameters 
if plotCC:
   plot_CC_violin(m0,tau,fit,r)

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

    
#%%
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

    FU.streamlines(coordE_f, Obs('r').array(), waterRes,
                   sensors=sensors, A=A, B=B, Nfix=Nfix)
    
    fig, ax = plt.subplots(nrows=1, ncols=4)
    for i, g in enumerate(range(1,20,5)):
        FU.streamlines(coordE_f, Obs('M'+str(g)).array(), waterRes,
                       sensors=sensors, A=A, B=B, Nfix=Nfix, ax=ax[i],
                       vmin=-10, vmax=10)
        ax[i].set_title('Gate t:' + str(IPcurves.t[g-1]) + 's')
        
    fig, ax = plt.subplots(nrows=1, ncols=4)
    for i, g in enumerate(range(1,20,5)):
        FU.streamlines(coordE_f, Obs('M'+str(g)).array(), waterRes,
                       sensors=sensors, A=A, B=B, Nfix=Nfix, ax=ax[i]
                       )
        ax[i].set_title('Gate t:' + str(IPcurves.t[g-1]) + 's')

        
#%%
if icsd:
    icsd_TDIP_plant(icsdPath,inputfile,all_gates)


