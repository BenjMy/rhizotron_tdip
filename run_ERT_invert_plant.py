# -*- coding: utf-8 -*-
"""
Created on Mon Sep 21 11:04:31 2020

@author: Benjamin
"""

# import libraries
from resipy import R2 # geophysics tools
import numpy as np # matrice operations
import pyvista as pv # if not installed : pip install with conda install pyvista (in a conda terminal)
import os 
from utils_rhizo import fct_utils as FU
import argparse
import pybert as pb 
import pygimli as pg
from pygimli.physics.ert import ERTManager, createGeometricFactors

# Import/read electrode geometry
# put a comment
main = os.getcwd()
os.chdir(main)
    
def invert_Resipy_ERT(inputfileERT):
        

    date = '1217' #
    inputfileERT = 'ERT_0113.csv' #'ERT_1217.csv'
    inputfileMALM = 'MALM_0113.csv' #'MALM_1217.csv'
    # Put the current folder
    #MainPath='E:/Padova/Experiments/2020_Rhizotron_Veronika_PRD_vitis_EIT/rhizo_gui/LAB/Using_code'
    #os.chdir(MainPath)
    # https://hkex.gitlab.io/resipy/api.html
    
    geomPath, meshPath, icsdPath = FU.definePath(main,date)
    
    # ---------------------------------------------------------#
    # create an instance of R2
    # ---------------------------------------------------------#
    k = R2(main, typ='R3t')
    
    # ---------------------------------------------------------#
    # create survey
    # ---------------------------------------------------------#
    # the csv file is exported from prosysII (using spreadsheet export with comma delimitated sep)
    # k.createSurvey(MainPath + '/data/10092020_water_ERT.csv', ftype='Syscal') # read the survey file
    k.createSurvey(main + './raw_data/' + inputfileERT, ftype='Syscal') # read the survey file
    
    k.filterUnpaired()
    k.filterRecip(percent=10) # in this case this only removes one quadrupoles with reciprocal error bigger than 20 percent
    
    k.fitErrorLin()
    
    k.showPseudo()
    k.showError() # plot the reciprocal errors
    
    
    # ---------------------------------------------------------#
    # Import mesh
    # ---------------------------------------------------------#
    # 1st solution: import electrodes geometry and create the mesh using Resipy meshtools
    # 2nd solution : import directly the mesh file (*.msh) --> produced using gmsh
    
    elecs= np.genfromtxt('./mesh/elecsXYZ.csv', delimiter=",",skip_header=1)
    
    elecsFlag = np.concatenate([np.zeros(8),np.ones(len(elecs[:,2])-8)])
    k.setElec(np.c_[elecs[:,0],elecs[:,2],elecs[:,1],elecsFlag])
    # k.setElec(np.c_[elecs[:,0],elecs[:,1],elecs[:,2]])
    
    k.importMesh('./mesh/mesh_rhizo_resipy.msh')
    #k.showMeshInParaview()
    
    # pl = pv.Plotter(notebook=False) # init pyvista plotter object
    # k.showMesh(ax=pl,xlim=[0,0.45],ylim=[-0.03,0],zlim=[0,0.5])
    
    # ---------------------------------------------------------#
    # inversion ERT
    # ---------------------------------------------------------#
    k.param['num_xy_poly'] = 0
    k.param['z_min'] = -100
    k.param['z_max'] = 100
    k.param['data_type'] = 1 # using log of resistitivy
    k.err = True # if we want to use the error from the error models fitted before
    k.invert() # this will do the inversion
    # k.saveData(MainPath)
    
    pl = pv.Plotter(notebook=False) # init pyvista plotter object
    k.showResults(ax=pl,
                  attr='Resistivity(log10)', 
                  sens=True, 
                  contour=True, 
                  use_pyvista=True,
                  xlim=[0,0.45],ylim=[-0.03,0],zlim=[0,0.5])
    
    pl = pv.Plotter(notebook=False) # init pyvista plotter object
    k.showResults(ax=pl,
                  attr='Resistivity', 
                  sens=True, 
                  contour=True, 
                  use_pyvista=True,
                  xlim=[0,0.45],ylim=[-0.03,0],zlim=[0,0.5])
    
    
    pl = pv.Plotter(notebook=False) # init pyvista plotter object
    k.showResults(ax=pl,
                  attr='Sensitivity_map(log10)', 
                  sens=True, 
                  contour=True, 
                  use_pyvista=True,pvshow=True,
                  xlim=[0,0.45],ylim=[-0.03,0],zlim=[0,0.5])
    
    #k.saveInvPlots(MainPath + 'figs')
    # k.showInvError() # all errors should be between -3 and 3
    
    # ---------------------------------------------------------#
    # MALM
    # ---------------------------------------------------------#
    #%% load observation data Robs
    #create a new survey MALM
    
    # k.createSurvey(main + './raw_data/' + inputfileMALM, ftype='Syscal') # read the survey file
    
    # # need to import sequence here? not imported during create survey?
    # seq = k.importSequence('./sequence/seq_10092020_water_MALM.csv')
    # #k.saveSequence(MainPath) # save to check sequence consistency
    # k.saveData(MainPath + '/test/')  # save to check data consistency
    
    
    # #%% Forward modelling of the solution
    # # produce an error: FATAL: Error in electrode definitions - labels defined in protocol.dat and R3t.in are incompatible
    # # problem is that the protocol.dat is written as a 2d protocol while it should be a 3d
    # # overcomed writting the protocol file by hand
    
    # k.forwardCSD(sources=[(0.0,0.0,0.3,1.0)], noise=0.05, iplot=False)
    
    # # question: 
    # # why the R2.in in the fwd folder is not identical to the one in the invdir folder?
    # # the mesh is unchanged from ERT to MALM (in the rhizotron case)
    
    # #%% create a regular grid of possible 
    # grid_vrte = [(0.075,0.35,2),(0,0,0),(0.075,0.45,2)]
    
    # #%% Inversion of CSD
    # k.invertCSD(grid=grid_vrte, x0=None, weightType='const', wreg=1)
    # k.showCSD()

    return k

def invert_pygimli_ERT(inputfileERT,sensors,mesh):
    #%% INVERT ERT data
    
    dataERT = pb.load(main + './raw_data/' + inputfileERT)
    dataERT.setSensorPositions(sensors)
    
    ert = pg.physics.ERTManager(dataERT)  # sr=False, verbose=True, debug=False)
    k = createGeometricFactors(dataERT)
    dataERT.set("k", k)
    dataERT.set('r', dataERT('u')/dataERT('i'))
    dataERT.set('rhoa', dataERT('r')*dataERT('k'))
    dataERT.set('rhoa', dataERT('r')*dataERT('k'))
    dataERT['err'] = ert.estimateError(dataERT, 
                                       absoluteError=0.001, 
                                       relativeError=0.03)
    dataERT.markInvalid(dataERT("rhoa") < 0)
    dataERT.removeInvalid()
    dataERT.save('dataERT.data')
    # ert.setMesh(mesh3d)  # important to do it explicitly as it avoids rebuilding region manager etc.
    # C = pg.matrix.GeostatisticConstraintsMatrix(mesh=ert.fop.paraDomain, I=[25, 5], dip=-25)
    # ert.fop.setConstraints(C)
    model = ert.invert(mesh=mesh,lam=100)
    
    return model


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='inputfileERT')
    parser.add_argument('--ERT', metavar='path', required=True,
                        help='inputfileERT')
    args = parser.parse_args()
    invert_Resipy_ERT(inputfileERT=args.ERT)
