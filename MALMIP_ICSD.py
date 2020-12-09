# -*- coding: utf-8 -*-
"""
Created on Tue Nov 10 16:56:48 2020
@author: Benjamin
"""

import os
# MALM Inversion library
from icsd3d_class import iCSD3d as i3d
from plotters.mpl_plot import showObs2d
from importers.read import *
import matplotlib.pyplot as plt
from utils_rhizo import fct_utils as FU


plt.close('all')
#%% 
A = 72-1 # Id electrode A
B = 65-1 # Id electrode B
Nfix = None # None #71-1 # Id electrode N

main = os.getcwd()
os.chdir(main)
meshPath= './mesh/'

icsdPath= './icsd/'
gateIP = 3 #False #3 
date = '0112' # 1310 1611
raw = True

icsdPath += date 
if raw: 
    icsdPath += '_raw' 
if Nfix: 
    icsdPath += '_Nfix' 
if gateIP: 
    icsdPath += '_M' + str(gateIP)
icsdPath += '/'
    
        
mesh3d, sensors = FU.mesh_import(meshPath + 'BaseRhizo_Vrte.msh')

#%%
icsd=i3d(dirName=icsdPath)   
icsd.type='2d'
icsd.obs_err='const' # choose between constant weight and w = 1/sqrt(abs(obs))
icsd.wr=1 #weight regularization
icsd.alphaSxy=True
icsd.x0_prior=False
icsd.x0_ini_guess=False # initial guess
icsd.icsd_init()

sol= icsd.invert(x0_prior=False,wr=0.1, showfig=True)

fig, ax = plt.subplots()
ax.scatter(sensors[:,0],sensors[:,1],color='k',marker='.',label='pot. elec')
ax.scatter(sensors[B,0],sensors[B,1],color='r',marker='v',label='B. elec')
ax.scatter(sensors[A,0],sensors[A,1],color='r',marker='v',label='A. elec')
if Nfix is not None:
    ax.scatter(sensors[Nfix,0],sensors[Nfix,1],color='g',marker='v',label='Nfix. elec')
for i in range(len(sensors[:,0])):
    ax.annotate(str(i+1), (sensors[i,0],sensors[i,1]))
icsd.showResults(ax=ax)
plt.title(icsdPath + 'Csd_gate'+ str(gateIP) + '.png')
plt.savefig(icsdPath + 'Csd_gate'+ str(gateIP) + '.png')
plt.show()