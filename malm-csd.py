# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 16:07:51 2020

@author: Benjamin
"""
#!pip install git+https://github.com/BenjMy/icsd_dev.git
!git clone --branch tdip https://github.com/BenjMy/icsd_dev.git
# first install missing packages


import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append('../src')
from resipy import R2

datadir = './dir_test'
k = R2(typ='R2',dirname=datadir)
elec = np.zeros((10,3))
elec[:,0] = np.arange(elec.shape[0])
k.setElec(elec)
k.createSequence()
k.createMesh('quad', cl=0.1)
k.addRegion(xz=[[2,-1],[6,-1],[6,-2],[2,-2]], res0=20)
k.showMesh(attr='res0')
k.forward()
k.showPseudo()

# create a regular grid of possible sources
gridx, gridy = np.meshgrid(0.5 + np.arange(9), -0.1 - np.arange(4))
gridx = gridx.flatten()
gridy = gridy.flatten()
grid = np.c_[gridx, gridy]
nodes = k.mesh.node

# snap grid points to mesh nodes
from scipy.spatial.distance import cdist
dist = cdist(grid, nodes[:,[0,2]])
imin = np.argmin(dist, axis=1)

fig, ax = plt.subplots()
k.showMesh(ax=ax)
ax.plot(nodes[imin,0], nodes[imin,2], 'ro')


# forward MALM (elec 2 is our source)
seqMALM = np.zeros((len(elec)-3, 4), dtype=int)
seqMALM[:,0] = 1 # elec far away
seqMALM[:,1] = 2 # tip of root (injected in trunk)
seqMALM[:,2] = 3 + np.arange(len(elec)-3)
seqMALM[:,3] = 4 + np.arange(len(elec)-3)
k.sequence = seqMALM

           
# adding multiple_sources
def multiple_sources(idS=[],wS=None):
    RobsSi = np.zeros((seqMALM.shape[0], len(idS)))
    if len(idS)>0:
        for i, s in enumerate(idS):
            k.param['node_elec'][1][1] = imin[s] # the source to be found (CHANGE HERE)
            k.forward()
            RobsSi[:,i] = k.surveys[0].df['resist'].values
        for i, s in enumerate(idS):
            if wS is None:
                RobsSi[:,i] = RobsSi[:,i]/len(idS)
            else:
                RobsSi[:,i] = RobsSi[:,i]*wS[i]
    Robs = np.sum(RobsSi,axis=1)
    return Robs

id_source2find = [13]
Robs = multiple_sources(idS=id_source2find,
                      wS=[0.5,0.5])

        
# adding noise
def addnoise(x, level=0.05):
    return x + np.random.randn(1)*x*level
addnoise = np.vectorize(addnoise)
Robs = addnoise(Robs, 0.01)


# build top part of the matrix A (time consuming part)
# compute response for each possible source
Rsim = np.zeros((seqMALM.shape[0], len(imin)))
Rtest = []
for i, d in enumerate(imin):
    print('\n=====', i, '++++', d)
    k.param['node_elec'][1][1] = d
    k.forward()
    Rsim[:,i] = k.surveys[0].df['resist'].values
    Rtest.append(k.surveys[0].df['resist'].values)
    
# export for iscd dev analysis
np.savetxt('VRTeSim.txt', np.hstack(Rtest), fmt='%1.2e')
np.savetxt('ObsData.txt', Robs,fmt='%4.5f')  
np.savetxt('VRTeCoord.txt', np.c_[nodes[imin,0],nodes[imin,2]],  delimiter='\t', fmt='%0.2f')  


from icsd3d_class import iCSD3d as i3d
from plotters import mpl_plot

#mpl_plot.showObs2d(path2files+'/')

path2files= datadir
icsd=i3d(dirName=path2files +'/')   
icsd.createSurvey(fname_obs='ObsData.txt',fname_sim='VRTeSim.txt')
icsd.regMesh='strc'
icsd.type='2d'
#icsd.mesh='invdir/fwd/forward_model.vtk'
icsd.obs_err='sqrt' # sqrt choose between constant weight and w = 1/sqrt(abs(obs))
icsd.wr=1 #weight regularization
icsd.alphaSxy=False
icsd.x0_prior=False
icsd.x0_ini_guess=False # initial guess
icsd.icsd_init()
# icsd.plotElecs=False
icsd.method_m0='F1'
m0 = icsd.estimateM0(method_m0='F1',show=True)
#m0 = icsd.estimateM0(method_m0='Pearson',show=True)

# icsd.clim=[0,0.1]
# icsd.run_single()
icsd.alphax0=1
sol= icsd.invert(wr=0.001)

fig, ax = plt.subplots()
icsd.showResults(ax=ax)
ax.scatter(nodes[imin[id_source2find],0],nodes[imin[id_source2find],2],
           color='r',marker='v')
plt.show()


# view results
fig, ax = plt.subplots()
k.showResults(ax=ax, color_map='Blues', edge_color='k', sens=False)
cax = ax.scatter(grid[:,0], grid[:,1], 50, x)
cax = ax.tricontourf(grid[:,0], grid[:,1], x)
fig.colorbar(cax, label='Source')
ax.scatter(nodes[imin[id_source2find],0],nodes[imin[id_source2find],2],color='r',marker='v')
