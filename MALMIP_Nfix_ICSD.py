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
from pygimli.meshtools import readGmsh

# My own library
from MALM4RootsRhizo_class import MALM4RootsRhizo_class as MR
#from utils import fct_utils as FU

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
#%% 

gateIP = 13 # id extracted time window in mV/V (?)
waterRes = 21.23 # water resistivity (Ohm.m)
rec = True
A = 72-1
B = 65-1
Nfix = 71-1


#%% PATHS (see excel ile for correspondance dates/files)

main = 'E:/Padova/Experiments/2020_Rhizotron_Veronika_PRD_vitis_EIT/'
os.chdir(main)

meshPath='E:/Padova/Experiments/2019_Rhizotron_DAFNAE_Lancaster_ERT_MALM/'
icsdPath= main + 'icsd/'

date = '1611'
inputfile = 'MALMIP1116f.bin'
if gateIP:
    icsdPath += date +'_M' + str(gateIP) +'/'
else:
    icsdPath += date +'/'


#%%
## --------- Import mesh --------- ##
fname = meshPath + '1_Mesh/Msh/BaseRhizo_Vrte.msh'
mesh3d=readGmsh(fname, verbose=True)
mesh3d.exportVTK('mesh3d_Markers')
#MR.PyMesh3d() # plot mesh3d_Markers marker in 3d

sensors = []
for node in mesh3d.nodes():
    if node.marker() == -99:
        sensors.append(node.pos())
    elif node.marker() == -999:
        print('-999')
        sensors.append(node.pos())
    elif node.marker() == -1000:
        print('-1000')
        sensors.append(node.pos())

sensors = np.vstack(sensors)


#%%
#abmn=[]
#for nn in range(len(m0)):
#    abmn.append([int(IPcurves.data(t)[nn]+1) for t in ['a', 'b', 'm', 'n']])
#abmn = np.vstack(abmn)


#%%
VRTEpos= MR.VRTEpos(mesh=mesh3d,dim=3) # return VRTE file format icsd MarkerVRTE=991
plt.scatter(VRTEpos[:,0],VRTEpos[:,1],color='b')

Obs_raw = pb.importer.importSyscalPro('./raw_data/' + inputfile) 
dataABMN = [np.array(Obs_raw('a'))+1, np.array(Obs_raw('b'))+1,
            np.array(Obs_raw('m'))+1,np.array(Obs_raw('n'))+1]
dataABMN = np.vstack(dataABMN).T

#%% select observation data

#indices = []
#for i in range(len(abmn)):
#    bool_select = dataABMN == abmn[i, None]
#    for j, x in enumerate(bool_select):
#        if x.all() == True:
#           indices.append(j)

#RmvE_Vec_Bool= np.ones([len(Obs_raw('a'))])
#RmvE_Vec_Bool[indices] = 0

Obs, dataABMN=  MR.PrepareMALMData('./raw_data/' + inputfile, Rec=True, DevErr=1,
                           MinV=1, MaxRc=1, Kfact=1, MinMaxAppRes=1, 
                           SwE=False) # remove reciprocal from original dataset

Obs_IP, dataABMN=  MR.PrepareMALMData('./raw_data/' + inputfile, Rec=True, DevErr=1,
                           MinV=1, MaxRc=1, Kfact=1, MinMaxAppRes=1, 
                           SwE=False, gIP=gateIP) # remove reciprocal from original dataset

#%%
SeqFullMALM= MR.mk_full_malm(dataABMN-1, VRTe = range(len(sensors), len(sensors)+len(VRTEpos)),
                             mesh=mesh3d, R3=True) # output .shm with sensors

rhomap = [[1, waterRes],
          [2, waterRes]]

rho = pg.solver.parseArgToArray(rhomap, mesh3d.cellCount(), mesh3d)
MR.SimulateGreenFcts(mesh_VRTs=mesh3d,rhomapVRTs=rho,schemeVrts=SeqFullMALM, 
                     Name='VRTeSim')

#%%
RemLineNb, Injection, coordE, pointsE= load_geom(main + 'geom/') # geometry file containing electrodes position includinf remotes 

#%%
icsd=i3d(dirName=icsdPath)   
icsd.type='2d'
icsd.obs_err='sqrt' # choose between constant weight and w = 1/sqrt(abs(obs))
icsd.wr=1 #weight regularization
icsd.alphaSxy=True
icsd.x0_prior=False
icsd.x0_ini_guess=False # initial guess
icsd.b=False
icsd.icsd_init()

sol= icsd.invert(x0_prior=False,wr=10, showfig=True)

fig, ax = plt.subplots()
ax.scatter(sensors[:,0],sensors[:,1],color='k',marker='.',label='pot. elec')
ax.scatter(sensors[B,0],sensors[B,1],color='r',marker='v',label='B. elec')
ax.scatter(sensors[A,0],sensors[A,1],color='r',marker='v',label='A. elec')
ax.scatter(sensors[Nfix,0],sensors[Nfix,1],color='g',marker='v',label='Nfix. elec')
for i in range(len(sensors[:,0])):
    ax.annotate(str(i+1), (sensors[i,0],sensors[i,1]))
icsd.showResults(ax=ax)
plt.savefig(icsdPath + 'Csd_gate'+ str(gateIP) + '.png')
plt.show()