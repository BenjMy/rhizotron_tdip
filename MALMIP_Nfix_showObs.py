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

#%% Define survey parameters

waterRes = 21.23 # water resistivity (Ohm.m)
rec = True # reciprocal analysis
A = 72-1 # Id electrode A
B = 65-1 # Id electrode B
Nfix = 71-1 # Id electrode N
gateIP = 3 # id extracted time window in mV/V (?)


#%% define PATHS (see excel file for correspondance dates/files)
# or set you working directory local path here 
main = os.getcwd()
os.chdir(main)

geomPath= './geom/'
meshPath= './mesh/'
icsdPath= './icsd/'

date = '1611'
inputfile = 'MALMIP1116.bin'
if gateIP:
    icsdPath += date +'_M' + str(gateIP) 
else:
    icsdPath += date 

#%% LOAD geometry and mesh
RemLineNb, Injection, coordE, pointsE= load_geom(geomPath) # geometry file containing electrodes position includinf remotes 
mesh3d, sensors = FU.mesh_import(meshPath + 'BaseRhizo_Vrte.msh')

#%% Import data and visualisation
IPcurves = tdip.TDIPdata('./raw_data/' + inputfile) # e.g. ABEM or Syscal TXT export
#IPcurves.data.setSensorPositions(coordE[:,1:3])
#IPcurves.filter(m=[71],n=[62])
IPcurves_f = tdip.TDIPdata('./raw_data/' + inputfile)

id_elec_2rmv = [42,44,61] # remove some noisy electrodes

for i in id_elec_2rmv:
    bool2rmv = IPcurves_f.data['m']==i
    if bool2rmv.any() == False: 
        print('search in n')
        bool2rmv = IPcurves_f.data['n']==i
    IPcurves_f.data.markInvalid(bool2rmv)
    IPcurves.data.markInvalid(bool2rmv)
IPcurves_f.MA = IPcurves_f.MA[:, IPcurves_f.data['valid'].array()==1]
IPcurves_f.data.removeInvalid()

IPcurves_f.showDecay(nr=np.arange(0,len(IPcurves_f.data['a'])), showFit=False, 
                   yscale='linear',xscale='linear')

m0,tau,fit = IPcurves_f.fitDecays(show=True)
IPcurves_f.showDecay(nr=np.arange(0,len(IPcurves_f.data['a'])), showFit=True, 
                   yscale='linear',xscale='linear')
m0_02,tau_02,fit_02 = IPcurves_f.fitDecays(show=True, tmin=0.2)


#%% plot Cole Cole model parameters 
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

#%% plot secondary voltage distribuion along a given profile
abmn=[]
for nn in range(len(IPcurves_f.data['a'])):
    abmn.append([int(IPcurves_f.data(t)[nn]+1) for t in ['a', 'b', 'm', 'n']])
abmn = np.vstack(abmn)

# define the profile
pelecs = np.array([5,13,21,29,37,45,53,61,69])
pelecs = np.array([0,1,2,3,4,5,6,7]) +32

idp = np.zeros(len(pelecs))
Vg =  np.zeros(20)
Vi = np.zeros([20,len(pelecs)])
for i,p in enumerate(pelecs):
    idp[i]= np.where(p==abmn)[0][0]
    for g in range(20):
        gatestr = 'M' + str(g+1)
        Vg[g] = IPcurves_f.data[gatestr].array()[int(idp[i])]
    Vi[:,i] = Vg

fig, ax = plt.subplots()
for g in range(3,20,5):
    plt.plot(pelecs+1,Vi[g,:],'o-', alpha=0.5, label='M'+str(g))
plt.legend()    
plt.xlabel('# Electrode')
plt.ylabel('Secondary voltage (mV)')
plt.xticks(pelecs+1)

#%%
fig, ax = plt.subplots()
ax.scatter(sensors[:,0],sensors[:,1],color='k',marker='.',label='pot. elec')
ax.scatter(sensors[B,0],sensors[B,1],color='b',marker='v',label='B. elec')
ax.scatter(sensors[A,0],sensors[A,1],color='r',marker='v',label='A. elec')
ax.scatter(sensors[Nfix,0],sensors[Nfix,1],color='g',marker='v',label='Nfix. elec')
for i in range(len(sensors[:,0])):
    ax.annotate(str(i+1), (sensors[i,0],sensors[i,1]))

ax.legend(loc="upper right")    
ax.set_xlabel('x (m)')
ax.set_ylabel('y (m)')

#%% select observation data
VRTEpos= FU.VRTEpos(mesh=mesh3d,dim=3) # return VRTE file format icsd MarkerVRTE=991
plt.scatter(VRTEpos[:,0],VRTEpos[:,1],color='b')

Obs_raw = pb.importer.importSyscalPro('./raw_data/' + inputfile) 
Obs, dataABMN=  FU.PrepareMALMData('./raw_data/' + inputfile, Rec=True, DevErr=1,
                           MinV=1, MaxRc=1, Kfact=1, MinMaxAppRes=1, 
                           SwE=False) # remove reciprocal from original dataset

Obs_f = Obs
Obs_f['valid'] = IPcurves.data['valid']
Obs_f.removeInvalid()
len(Obs_f['M1'])
dataABMN_f = dataABMN[IPcurves.data['valid'].array()==1,:]

Obs.setSensorPositions(coordE[:,1:3])
Obs.save('obs.data')

Obs_f['m'] 
Obs_f['n'] 

coordE_f = []
for i, mi in enumerate(Obs_f['m']):
    if mi==Nfix:
       mi=Obs_f['n'][i]
    id_coordE_f = np.where(mi+1==coordE[:,0])[0]
    coordE_f.append(coordE[id_coordE_f[0],:])
coordE_f = np.array(coordE_f)
        
#%%

fig, ax = plt.subplots(nrows=1, ncols=2)
sc=ax[0].scatter(coordE_f[:,1], coordE_f[:,2], c=abs(Obs_f('r').array()), 
              cmap ='coolwarm',s=5e2, vmin=0.1, norm=matplotlib.colors.Normalize())
cbar = plt.colorbar(sc,ax=ax[0])
cbar.set_label('rhoa (Ohm.m)')   
ax[0].scatter(sensors[:,0],sensors[:,1],color='k',marker='.',label='pot. elec')
ax[0].scatter(sensors[B,0],sensors[B,1],color='b',marker='v',label='B. elec')
ax[0].scatter(sensors[A,0],sensors[A,1],color='r',marker='v',label='A. elec')
ax[0].scatter(sensors[Nfix,0],sensors[Nfix,1],color='g',marker='v',label='Nfix. elec')
for i in range(len(sensors[:,0])):
    ax[0].annotate(str(i+1), (sensors[i,0],sensors[i,1]))   
ax[0].set_ylabel('y [m]',fontsize=15)
ax[0].set_xlabel('x [m]',fontsize=15)
sc=ax[1].scatter(coordE_f[:,1], coordE_f[:,2], c=abs(Obs_f('M'+str(gateIP)).array()), 
              cmap ='coolwarm',s=5e2, norm=matplotlib.colors.Normalize())
cbar = plt.colorbar(sc,ax=ax[1]) 
    


#%%

fig, ax = plt.subplots(nrows=1, ncols=4)
for i, g in enumerate(range(3,20,5)):
    ax[i].set_ylabel('y [m]',fontsize=15)
    ax[i].set_xlabel('x [m]',fontsize=15)
    sc=ax[i].scatter(coordE_f[:,1], coordE_f[:,2], c=abs(Obs_f('M'+str(g)).array()), 
                  cmap ='coolwarm',s=5e2)
    cbar = plt.colorbar(sc,ax=ax[i]) 
    ax[i].set_title('Gate t:' + str(IPcurves.t[g]) + 's')
    for ie in range(len(sensors[:,0])):
        ax[i].annotate(str(ie+1), (sensors[ie,0],sensors[ie,1]))

    