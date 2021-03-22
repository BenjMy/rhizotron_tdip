# coding: utf-8

# 1st step
# In[1]:


# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 11:04:31 2020
@author: Benjamin
"""
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import matplotlib.pyplot as plt

# import libraries
from resipy import R2 # geophysics tools
import numpy as np # matrice operations
import pyvista as pv # if not installed : pip install with conda install pyvista (in a conda terminal)
import shutil
import glob

# Put the current folder
import os 
#MainPath=os.getcwd()
MainPath = 'E:/Padova/Experiments/2020_Rhizotron_Veronika_PRD_vitis_EIT/git_TDIP/'
os.chdir(MainPath)
#https://hkex.gitlab.io/resipy/api.html


# In[2]:


elecs_raw = np.genfromtxt('./mesh/elecsXYZ.csv', delimiter=",",skip_header=1)
elecs = np.copy(elecs_raw)


# In[3]:


# ---------------------------------------------------------#
# create an instance of R3 for each survey
# ---------------------------------------------------------#
csvfiles = glob.glob1(MainPath  + 'raw_data/ERT/',"*.csv")

for i in enumerate([csvfiles[0],csvfiles[6]]):

    k = R2(MainPath + 'ERT_inversion/step' + str(i[0]), typ='R3t')
    k.setTitle('Rhizo_'+str(i[0]))
    # ---------------------------------------------------------#
    # create survey
    # ---------------------------------------------------------#

    k.createSurvey(MainPath + 'raw_data/ERT/'+ str(i[1]) , ftype='Syscal') # read the survey file
    k.filterUnpaired()
    k.filterRecip(percent=10) # in this case this only removes one quadrupoles with reciprocal error bigger than 20 percent
    k.fitErrorPwl()
    # k.setElec(np.c_[elecs[:,0],elecs[:,1],elecs[:,2],elecsFlag])
    k.setElec(np.c_[elecs[:,0],elecs[:,1],elecs[:,2]])
    #k.filterAppResist(vmin=-2000,vmax=2000)
    # ---------------------------------------------------------#
    # import mesh
    # ---------------------------------------------------------#
    elecs= np.genfromtxt('./mesh/elecsXYZ.csv', delimiter=",",skip_header=1)
    
    elecsFlag = np.concatenate([np.zeros(8),np.ones(len(elecs[:,2])-8)])
    k.setElec(np.c_[elecs[:,0],elecs[:,2],elecs[:,1],elecsFlag])
    # k.setElec(np.c_[elecs[:,0],elecs[:,1],elecs[:,2]])
    
    k.importMesh('./mesh/mesh_rhizo_resipy.msh')
    # ---------------------------------------------------------#
    # inversion
    # ---------------------------------------------------------#
    k.param['num_xy_poly'] = 0
    k.param['zmin'] = -np.inf
    k.param['zmax'] = np.inf
    k.param['data_type'] = 1 # using log of resistitivy
    k.err = False # if we want to use the error from the error models fitted before
    k.param['a_wgt'] = 0.001
    k.param['b_wgt'] = 0.05
    k.invert() # this will do the inversion
    # k.saveData(MainPath)

    # pl = pv.Plotter(notebook=False) # init pyvista plotter object
    # fig = k.showResults(ax=pl,index=-1,pvcontour=[2, 2.1], 
    #               pvgrid=True,
    #               attr='difference(percent)',
    #               sens=True,
    #               # vmin=-10, vmax=10,
    #               color_map='bwr',
    #               xlim=[0,100],ylim=[-50,50],zlim=[-100,0],
    #               pvshow=False)
    # pl.show(screenshot=MainPath +  'TLdiffstep_vol' + str(tl) + '.png')
    
    pl = pv.Plotter(notebook=True) # init pyvista plotter object
    k.showResults(ax=pl,
                  attr='Resistivity(log10)', 
                  sens=True, 
                  contour=True, 
                  use_pyvista=True,
                  xlim=[0,0.45],ylim=[-0.03,0],zlim=[0,0.5],
                  pvshow=False)
    pl.show(screenshot=MainPath +  'ERT_inversion/step' + str(i[0]) + '/invdir/'+ i[1] + '_log.png')

    pl = pv.Plotter(notebook=True) # init pyvista plotter object
    k.showResults(ax=pl,
                  attr='Resistivity(ohm.m)', 
                  sens=True, 
                  contour=True, 
                  vmin=8, vmax=10,
                  use_pyvista=True,
                  xlim=[0,0.45],ylim=[-0.03,0],zlim=[0,0.5],
                  pvshow=False)
    pl.show(screenshot=MainPath +  'ERT_inversion/step' + str(i[0]) + '/invdir/'+ i[1] + '.png')
    
    # ---------------------------------------------------------#
    # copy protocol file to TL folder
    # ---------------------------------------------------------#
    src=MainPath +  'ERT_inversion/step' + str(i[0]) + '/invdir/protocol.dat'
    dst_folder=MainPath + 'ERT_inversion/TL/'
    
    if not os.path.exists(dst_folder):
        os.makedirs(dst_folder)
    shutil.copy(src,dst_folder + 'protocol' + str(i[0]) + '.dat')
    


# In[ ]:
k = R2(MainPath + 'ERT_inversion/TL' , typ='R3t')
k.setTitle('Rhizo_inv_TL')

# ---------------------------------------------------------#
# create survey
# ---------------------------------------------------------#
testdir = MainPath + 'ERT_inversion/TL/'
k.createTimeLapseSurvey(testdir, ftype='ProtocolDC')
k.importElec(testdir + '../electrodes.csv')
k.importMesh('./mesh/mesh_rhizo_resipy.msh')

# ---------------------------------------------------------#
# inversion
# ---------------------------------------------------------#
#k.invert(parallel=True) # takes a while because it invert all the surveys together

k.param['num_xy_poly'] = 0
# k.param['z_min'] = -100
# k.param['z_max'] = 0

k.param['zmin'] = -np.inf
k.param['zmax'] = np.inf

k.param['data_type'] = 1 # using log of resistitivy
k.err = False # if we want to use the error from the error models fitted before
k.param['a_wgt'] = 0.01
k.param['b_wgt'] = 0.05
# k.saveData(MainPath)
k.reg_mode = 1
k.invert() # this will do the inversion

# In[ ]:

for tl in range(len(k.surveys)):
    pl = pv.Plotter(notebook=False) # init pyvista plotter object
    fig = k.showResults(ax=pl,index=tl,
                  pvgrid=True,
                  attr='difference(percent)',
                  sens=True,
                  # vmin=-10, vmax=10,
                  color_map='bwr',
                  xlim=[0,100],ylim=[-50,50],zlim=[-100,0],
                  pvshow=False)
    pl.show(screenshot=MainPath + 'ERT_inversion/TL/TLdiffstep' + str(tl) + '.png')
    pl.close()


for tl in range(len(k.surveys)):
    pl = pv.Plotter()
    k.showResults(ax=pl, sens=True, index=tl,  attr= 'Resistivity(ohm.m)',
                  vmin=0, vmax=500,
                  pvslices=([50],[0]), pvgrid=True)
    pl.show(screenshot=MainPath +  'TLstep' + str(tl) + '.png')


# In[ ]:




