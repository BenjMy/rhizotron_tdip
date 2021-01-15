# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 11:33:33 2021
@author: Benjamin
Run icsd TDIP from data collected in the rhizotron
"""

import os 
from icsd3d_class import iCSD3d as i3d
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def icsd_TDIP_plant(icsdPath,inputfile,all_gates):
       
    #mpl_plot.showObs2d(path2files+'/')
    sep = '_'
    NameSave = 'O'+ os.path.basename(inputfile).split(sep, 1)[0] + '.txt'

    icsd=i3d(dirName=icsdPath +'/')   
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