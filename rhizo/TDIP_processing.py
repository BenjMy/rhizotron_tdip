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
#from run_ICSD import icsd_TDIP_plant
from run_ERT_invert_plant import invert_Resipy_ERT, invert_pygimli_ERT
# General libraries
import numpy as np
import matplotlib.pyplot as plt

# Visualisation libraries
import pyvista as pv

# MALM Inversion library
#from icsd3d_class import iCSD3d as i3d
#from plotters.mpl_plot import showObs2d
from importers.read import *

plt.close('all')

#%% Define survey parameters (see excel file for correspondance dates/files/parameters)
invERT = True
waterRes = (1/855)*1e4 #21.23 # water resistivity (Ohm.m) converted from microS/cm 767, 477 
rec = True # reciprocal analysis
A = 72-1 # Id electrode A
B = 65-1 # Id electrode B
injection_duration = 2 # time of injection
# MALMIP_0122 MALMIP_0113 'MALMIP1217' 'MALMIP1201.bin' 'MALMIP1013' 'MALMIP1116.bin' # filename in raw data folder
date = '0218' #  '1712' 0112 '1310' '1611' # (ddmm)
inputfileMALM = 'MALMIP_' + date + '.bin' #  
inputfileERT = 'ERT_' + date + '.bin' #
plotCC = False # show Cole-Cole fitted parameters (to compare with literature)
split_Nfix = [True, 71-1]
Nfix = 71-1 #71-1  #None 71-1 # Id electrode N , put None if N is varying

rmvInvalid = False # True if you want to write filtered files/ False for raw data
rmv_outliers = True
rmv_id= None #None # 61-1 o None
rmv_dup = True

all_gates= True
if not all_gates: 
    gateIP = 3 
else: gateIP = None # id extracted time window in mV/V (?)

# icsd = True

main = os.getcwd()
os.chdir(main)
geomPath, meshPath, icsdPath, figpath = FU.definePath(main,date)

    
#%% LOAD geometry and mesh
RemLineNb, Injection, coordE, pointsE= load_geom(geomPath) # geometry file containing electrodes position includinf remotes 
mesh3d_fwd, sensors = FU.mesh_import(meshPath + 'BaseRhizo_Vrte.msh')
mesh3d_inv, sensors = FU.mesh_import(meshPath + 'BaseRhizo_Vrte_inv.msh')

#%% INVERT ERT data

#if invERT:
#    inputfileERTcsv = 'ERT_' + date + '.csv'
#    k = invert_Resipy_ERT(date, inputfileERTcsv)
#    k.showResults()
#    k.saveInvPlots(figpath + 'figs')
       
    
    
if invERT:
    #inputfileERT = 'ERT_0122.bin'
    model = invert_pygimli_ERT(inputfileERT,sensors,mesh3d_inv,date)
    pg.show(mesh3d_inv,data=model,notebook=True)
    mesh3d_inv.addData('model',model)
    mesh3d_inv.exportVTK(figpath +'model' + date + '.vtk')
    
    plotter, _ = pg.show(mesh3d_inv, data=model,
                         alpha=0.9, hold=True, notebook=True)
    plotter.view_xy()
    #plotter.clim([20, 60])
    plotter.show()
    plotter.screenshot(figpath + 'model' + date + '.png')
#plt.imshow(plotter.image)
#plt.show()

#%% Import data TDIP
IPcurves = tdip.TDIPdata('./raw_data/' + inputfileMALM) # e.g. ABEM or Syscal TXT export
valid = np.ones(len(IPcurves.data('m')))

fig, ax = plt.subplots()
ax = IPcurves.showDecay(nr=np.arange(0,len(IPcurves.data['a'])), showFit=False, 
                   yscale='linear',xscale='linear', ax=ax)
plt.savefig(figpath + 'rawdecay' + date + '.png')

#%% split and filter to Nfix
# ----------------
IPcurves_f = tdip.TDIPdata('./raw_data/' + inputfileMALM)

if split_Nfix[0]:

    IPcurves.data('m').array()
    idn = np.where(IPcurves.data('m')==split_Nfix[1])[0]
    idm = np.where(IPcurves.data('n')==split_Nfix[1])[0]
    print(idn)
    print(idm)
    idfix = list(idn) #+ list(idm) #list(idn)  #+ list(idm)
    
    IPcurves.data('a')[idfix].array()
    IPcurves.data('b')[idfix].array()
    id_elec_2rmv = idfix # remove Nfix electrodes
    
    if Nfix is not None:
        a = set(list(range(0, len(IPcurves.data('m')))))
        id_elec_2rmv = a.difference(set(idfix))
        id_elec_2rmv = list(id_elec_2rmv)
    
    IPcurves_f, valid_split = FU.filterTDIP(IPcurves,id_elec_2rmv)

# remove duplicates 
# ----------------
if rmv_dup:
    a = IPcurves.data('a').array()
    abmn = np.zeros((len(a),4))
    abmn[:,0] = IPcurves.data('a').array()
    abmn[:,1] = IPcurves.data('b').array()
    abmn[:,2] = IPcurves.data('m').array()
    abmn[:,3] = IPcurves.data('n').array()
    unique, index_unique = np.unique(abmn, axis=0, return_index=True)
    test = set(list(range(0, len(IPcurves.data('a').array()))))
    id_elec_2rmv = test.difference(set(index_unique))
    id_elec_2rmv = list(id_elec_2rmv)
    IPcurves_f, valid_dup = FU.filterTDIP(IPcurves,id_elec_2rmv)

# remove outliers 
# ----------------
if rmv_outliers: 
    id_outliers = np.where(abs(IPcurves.data['M1'])<10)[0]
    test = set(list(range(0, len(IPcurves.data('a').array()))))
    id_elec_2rmv = test.difference(set(id_outliers))
    id_elec_2rmv = list(id_elec_2rmv)
    IPcurves_f, valid_outliers = FU.filterTDIP(IPcurves,id_elec_2rmv)

#if rmv_id: 
IPcurves.data['M1'].array()
        
# try:
#     valid_split
# except NameError:
#     valid_split = np.ones(len(IPcurves.data('m')))
    
j=0
for i, v in enumerate(valid):
    if split_Nfix[0]:
        if valid_split[j] == 0:
            valid[i] = 0
    if rmv_outliers:
        if valid_outliers[j] == 0:
            valid[i] = 0
    if rmv_dup:
        if valid_dup[j] == 0:
            valid[i] = 0
    j = j + 1

np.count_nonzero(valid==0)

IPcurves_f.data.set('valid',valid)
IPcurves_f.MA = IPcurves_f.MA[:, IPcurves_f.data['valid'].array()==1]
IPcurves_f.data.removeInvalid()
IPcurves_f.data.save('TDIP_filtered.data')

#%% Show decay
# ----------------
fig, ax = plt.subplots()
ax = IPcurves_f.showDecay(nr=np.arange(0,len(IPcurves_f.data['a'])), showFit=False, 
                   yscale='linear',xscale='linear', ax=ax)
plt.savefig(figpath + 'filtered_decay' + date + '.png')

if rmvInvalid is False:
    IPcurves_f = tdip.TDIPdata('./raw_data/' + inputfileMALM)
    
#%% plot Cole Cole model parameters 
# ----------------
if plotCC:
   plot_CC_violin(m0,tau,fit,r)

#%% select observation data - in preparaion for ICSD analysis --> 3 files required (same folder)
# 1st file needed = virtual electrode position
VRTEpos= FU.VRTEpos(mesh=mesh3d_fwd,dim=3) # return VRTE file format icsd MarkerVRTE=991
#plt.scatter(VRTEpos[:,0],VRTEpos[:,1],color='b')

# 2nd file = observations data
Obs_raw = pb.importer.importSyscalPro('./raw_data/' + inputfileMALM) 
Obs, dataABMN=  FU.PrepareMALMData('./raw_data/' + inputfileMALM, Rec=False, DevErr=1,
                           MinV=1, MaxRc=1, Kfact=1, MinMaxAppRes=1, 
                           SwE=False, 
                           valid=valid,savefile=True)

# 2nd file = observations data = TDIP
Obs_IP, dataABMN=  FU.PrepareMALMData('./raw_data/' + inputfileMALM, Rec=False, DevErr=1,
                           MinV=1, MaxRc=1, Kfact=1, MinMaxAppRes=1, 
                           SwE=False, gIP=gateIP, 
                           valid=valid) # remove reciprocal from original dataset

if all_gates:
    for i, g in enumerate(range(1,20,1)):
        Obs_IP, dataABMN=  FU.PrepareMALMData('./raw_data/' + inputfileMALM, Rec=False, DevErr=1,
                                   MinV=1, MaxRc=1, Kfact=1, MinMaxAppRes=1, 
                                   SwE=False, gIP=g, 
                                   valid=valid) # remove reciprocal from original dataset

coordE_f = []
for i, mi in enumerate(Obs['m']):
    if mi==Nfix:
       mi=Obs['n'][i]
    id_coordE_f = np.where(mi+1==coordE[:,0])[0]
    #if len(id_coordE_f) > 1:
    coordE_f.append(coordE[id_coordE_f[0],:])
coordE_f = np.array(coordE_f)

# filter out wrong observation 
# if (rmvInvalid or split_Nfix[0]) is False:
#     print('raw analysis')
#     Obs, dataABMN=  FU.PrepareMALMData('./raw_data/' + inputfileMALM, Rec=False, DevErr=1,
#                                MinV=1, MaxRc=1, Kfact=1, MinMaxAppRes=1, 
#                                SwE=False,
#                                valid=None)
#     Obs_IP, dataABMN=  FU.PrepareMALMData('./raw_data/' + inputfileMALM, Rec=False, DevErr=1,
#                                MinV=1, MaxRc=1, Kfact=1, MinMaxAppRes=1, 
#                                SwE=False, gIP=gateIP, 
#                                valid=None) # remove reciprocal from original dataset
#     coordE_f = []
#     for i, mi in enumerate(Obs['m']):
#         if mi==Nfix:
#            mi=Obs['n'][i]
#         id_coordE_f = np.where(mi+1==coordE[:,0])[0]
#         coordE_f.append(coordE[id_coordE_f[0],:])
#     coordE_f = np.array(coordE_f)

#%% foward modelling of homogeneous medium
# interpolate to MALM fwd mesh
if invERT:
    rho_fwd = pg.interpolate(mesh3d_fwd, mesh3d_inv, model, method='spline')
    rhomap = rho_fwd.array()
    id0 = np.where(rho_fwd.array() <=0)
    rhomap[id0] = 0.1
else:
    model = waterRes
    rhomap = [[1, waterRes],
              [2, waterRes]]


#%% MR.SimulateSynthMALMSol(mesh=mesh3d_fwd,SeqSol=dataABMN-1,rhomap=model)

#%% foward modelling of Green's functions
# mk_full_malm --> 3rd file = simulated green functions (simulated data to compare with observation data)
SeqFullMALM= MR.mk_full_malm(dataABMN-1, 
                             VRTe = range(len(sensors),
                                          len(sensors)+len(VRTEpos)),
                             mesh=mesh3d_fwd, 
                             R3=False) # output .shm with sensors

#pg.show(mesh3d_fwd,data=rhomap,notebook=True)
#pg.show(mesh3d_inv,data=rhomap,notebook=True)

MR.SimulateGreenFcts(mesh_VRTs=mesh3d_fwd,rhomapVRTs=rhomap,schemeVrts=SeqFullMALM, 
                     Name='VRTeSim')


#%% Quiver plot (gradient*conductivity)
if Nfix is not None:

    FU.streamlines(coordE_f, Obs('r').array(), model,
                   sensors=sensors, A=A, B=B, Nfix=Nfix,
                    vmin=-300, vmax=200, mesh_inv=mesh3d_inv)
    plt.savefig(figpath +'streamlines_PV.png')

    fig, ax = plt.subplots(nrows=1, ncols=4,figsize=(20,5))
    for i, g in enumerate(range(1,20,5)):
        FU.streamlines(coordE_f, Obs('M'+str(g)).array(), model,
                       sensors=sensors, A=A, B=B, Nfix=Nfix, ax=ax[i],
                       vmin=-10, vmax=10, mesh_inv=mesh3d_inv)
        ax[i].set_title('Gate t:' + str(IPcurves.t[g-1]) + 's')
        plt.tight_layout()
        plt.savefig(figpath +'streamlines_transients.png')

    fig, ax = plt.subplots(nrows=1, ncols=4,figsize=(20,5))
    for i, g in enumerate(range(1,20,5)):
        FU.streamlines(coordE_f, Obs('M'+str(g)).array(), model,
                       sensors=sensors, A=A, B=B, Nfix=Nfix, ax=ax[i]
                       , mesh_inv=mesh3d_inv)
        ax[i].set_title('Gate t:' + str(IPcurves.t[g-1]) + 's')
        plt.tight_layout()
        plt.savefig(figpath + 'streamlines_transients2.png')


#%%  

#if icsd:
    #%% copy file to icsd folder
    #import glob
    #import shutil
    
    #for fpath in glob.glob('OMALM*.txt'):  # this returns a list of the CURRENT contents. Therefore, there is no need to sanitize for the directories that we will later create
    #    if not os.path.isdir(icsdPath):
    #        os.mkdir(icsdPath)
        #shutil.move(fpath, icsdPath)
    #    shutil.move(os.path.join(os.getcwd(), fpath), os.path.join(icsdPath, fpath))
    
    #shutil.copy('VRTeSim.txt', icsdPath)
    #shutil.copy('VRTeCoord.txt', icsdPath)

#    icsd_TDIP_plant(icsdPath,inputfile,all_gates,IPcurves)


