# -*- coding: utf-8 -*-
"""
@author: B. Mary
Synthetic modelling of secondary voltages (decay curves in TDIP)
"""

#%% libraries

import matplotlib.pyplot as plt
import numpy as np
import pybert as pb
import pygimli as pg
import pygimli.meshtools as mt
import pygimli.physics.ert as ert
from pybert import tdip
from pybert.tdip import mipmodelling
import pandas as pd

plt.close('all')
#%% define survey geometry
# Here take a simple 2d line for simplicity - after use rhizotron set-up config

world = mt.createWorld(
    start=[-55, 0], end=[105, -80], worldMarker=True)

conductive_anomaly = mt.createCircle(
    pos=[10, -7], radius=5, marker=2
)

polarizable_anomaly = mt.createCircle(
    pos=[40, -7], radius=5, marker=3
)    
    
geom = world + conductive_anomaly + polarizable_anomaly
pg.show(geom)

scheme = ert.createERTData(
    elecs=np.linspace(start=0, stop=50, num=51),
    schemeName='dd'
)


#%% Create a mesh for the finite element modelling with appropriate mesh quality.

# local refinement of mesh near electrodes
for s in scheme.sensors():
    geom.createNode(s + [0.0, -0.2])
    
# mesh = mt.createMesh(geom, quality=34)
mesh_coarse = mt.createMesh(geom, quality=33)
mesh = mesh_coarse.createH2()
    
#%% # Create a map to set resistivity values in the appropriate regions
rhomap = [[1, 100.],
          [2, 75.],
          [3, 150.]]

# Take a look at the mesh and the resistivity distribution
pg.show(mesh, data=rhomap, label=pg.unit('res'), showMesh=True)

#%% Simulate ERT response

data_ert = ert.simulate(mesh, scheme=scheme, res=rhomap, noiseLevel=1,
                    noiseAbs=1e-6, seed=1337)

data_ert.remove(data_ert['rhoa'] < 0)
#pg.info('Filtered rhoa (min/max)', min(data['rhoa']), max(data['rhoa']))
data_ert.save('simple.dat')

#%% Mesh generation for the inversion
world = mt.createWorld(
    start=[-15, 0], end=[65, -30], worldMarker=False, marker=2)

# local refinement of mesh near electrodes
for s in scheme.sensors():
    world.createNode(s + [0.0, -0.4])

mesh_coarse = mt.createMesh(world, quality=33)
mesh = mesh_coarse.createH2()
for nr, c in enumerate(mesh.cells()):
    c.setMarker(nr)
pg.show(mesh)


#%% Initialize the ERTManager and invert ERT

mgr = ert.ERTManager('simple.dat')
model = mgr.invert(lam=20, verbose=True)

modelPD = mgr.paraModel(model)  # do the mapping
meshPD = pg.Mesh(mgr.paraDomain) # Save copy of para mesh for plotting later

mgr.showResultAndFit()

pg.show(mgr.paraDomain, modelPD, label='Model', cMap='Spectral_r',
        logScale=True, cMin=25, cMax=150)

#%% Regions for chargeability.

charg_ano = 0.005
count = 0
charg = np.zeros(mgr.paraDomain.cellCount())
for c in mgr.paraDomain.cells():
    if (c.center()[0] <=40+5) & (c.center()[0] >=40-5)  & (c.center()[1]<=-7+5) & (c.center()[1] >=-7-5):
        charg[count] = charg_ano
    else:
        charg[count] = .0001
    count = count+1


fig, axes = plt.subplots(1, 2, figsize=(16 / 2.54, 16 / 2.54))
pg.show(mgr.paraDomain, modelPD, label='Model', cMap='Spectral_r',
        logScale=True, cMin=25, cMax=150, ax=axes[0])
pg.show(meshPD, data=charg, ax=axes[1], label=r"$M''$~[?]")
fig.tight_layout()
fig.show()

#%% TDIP modelling

fop=mgr.fop # take the ERT fop
TDmanager_M = mipmodelling.DCIPMModelling(f=fop,mesh=mgr.paraDomain,rho=model) 
# DCIPMModelling only return integral chargeability

#%% Return forward response as function of chargeability m
charPs = TDmanager_M.response(m=charg)*1000
data_mint = data_ert
data_mint.add("mint",charPs)     #add a column of IP(mV/V)  mint for integral, ip for bert
data_mint.save('mint.data')

#%% visu results mint
tdip_mint=tdip.TDIPdata(data=data_mint)
tdip_mint.showRhoa()
tdip_mint.showIntegralChargeability()

tdip_mint.invertRhoa(data=data_mint)
tdip_mint.showResistivity(logScale=True)

#%% Simultaneous IP modelling of several gates using DCIPMSmoothModelling
# define gate times
from pybert import tdip
IPcurves = tdip.TDIPdata('MALMIP1217.bin') # e.g. ABEM or Syscal TXT export

# t = np.array([0.1,0.3,0.6,1.2])
t = IPcurves.t

TDmanager_gates = mipmodelling.DCIPMSmoothModelling(f=fop,mesh=meshPD,
                                                    rho=modelPD,t=t) 

#%% TENTATIVE 1 #
#%% Define M(t) in 1/t^2 to inject into DCIPMSmoothModelling.response

chargMw = np.ones([len(t),1])*charg  # nb of times * nb of cells in the mesh
for i in range(len(t)):
    chargMw[i]=charg/(t[i]*t[i]) # decay in 1/t^2

#%% Return forward response as function of chargeability M(t) arbitrary defined

charPs = TDmanager_gates.response(m=chargMw)*1000
charPsM = np.reshape(charPs,[len(t),len(data_ert['a'])]) # (nb of times * length of sequence)
data_mi = data_ert

for tt in range(len(t)):
    data_mi.add("m"+str(tt),charPsM[i,:])     #add a column of IP(mV/V)  for each gates
    
data_mi.save('all_gates.dat')  #save in this file

#%% visu results for all gates

IPcurves=tdip.TDIPdata(data=data_mi,t=t, MA=charPsM)
IPcurves.MA=charPsM
IPcurves.t = t
IPcurves.showRhoa()
IPcurves.showMa(3)
IPcurves.showIntegralChargeability()

IPcurves.showDecay(nr=np.arange(0,len(IPcurves.data['a'])), showFit=False, 
                   yscale='linear',xscale='linear')

#%% TENTATIVE 2. Synthetic simulation based on Cole-Cole model #
# Return forward response as function of chargeability M(t) 
# using ColeCole modelling of decay curves


# define the model
tauvec = np.array([1e-3, 0.1])  # s
mvec = np.array([0.001, 0.1])  #
cvec = np.array([0.5,  0.5])  # 

# model = (mvec,tauvec)
model = (mvec,tauvec,cvec)

# Define M(t)
CC = mipmodelling.CCTDModelling(t) 

# !CC.Response produce a shape error!
Z = CC.response(model) # return complex resistivity (?)


t = np.array([0.1,0.3,0.6,1.2])
TDmanager_gates_CC = mipmodelling.DCIPMSmoothModelling(f=fop,
                                                       mesh=meshPD,
                                                       rho=modelPD,
                                                       t=t) 
charPs = TDmanager_gates.response(m=Z)*1000

