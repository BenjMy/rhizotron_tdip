import os
# My own library
from utils_rhizo import fct_utils as FU
from MALM4RootsRhizo_class import MALM4RootsRhizo_class as MR
from importers.read import *

#%% define PATHS 
# or set you working directory local path here 
main = os.getcwd()
os.chdir(main)

geomPath= './geom/'
meshPath= './mesh/'


#%% LOAD geometry and mesh
RemLineNb, Injection, coordE, pointsE= load_geom(geomPath) # geometry file containing electrodes position includinf remotes 
mesh3d, sensors = FU.mesh_import(meshPath + 'BaseRhizo_Vrte.msh') # BaseRhizo_Splitted_Vrte_v2.msh'

#%% sequences for full rhizo 

# MALM sequence for full rhizo
SeqMALM = FU.mk_seqMALM(b_return=65, sensors=sensors, savefig=True)

# MALM sequence for full rhizo with moving B and fix M
SeqMALM = FU.mk_seqMALM(b_return=65, sensors=sensors, 
                        mfix=True, mref=71, savefig=True)

# MALM sequence for full rhizo with moving B
# FU.mk_seqMALM(b_return=65, mesh=mesh3d, Ch_return=[1,3])

# # MALM sequence for full rhizo with moving B and fix M
# FU.mk_seqMALM(b_return=65, mesh=mesh3d, Ch_return=[1,3], mfix=71)

#%% sequences for split rhizo 
# ERT

## create sequence ERT36R left and right and merged
# [Left, Right]= MR.ElecLeftRight()
# #schemeERTSplit, ABMNSplit= MR.CreateSchemeERT(mesh3d, Rec=True, Syscal=True,
# #                                       PosXYZid=[Left, Right],
# #                                       savefig=False,SaveName='SchemeSplitERT') 
# #    
            
# # Alternative use :
# #schemeERTSplit= MR.CreateSchemeERTSplit(mesh3d, Rec=True, Syscal=True,
# #                                       PosXYZid=[Left, Right])

# ## --------- Prepare MALM sequences --------- ##
# #SeqMALM= MR.mk_seq(a_injection = 72, b_return = [65],
# #                   mesh=mesh3d, Syscal=True, Rec=False) 
# #len(SeqMALM)
# #
# ## create sequence MALM left and right and merged
# SeqMALMSplit= MR.mk_seqMALM(a_injection = 72, b_return = 71,
#                    mesh=mesh3d, Syscal=True, Rec=True, 
#                    PosXYZid=[Left+1, Right+1],Ch_return=65,
#                    savefig=False, SaveName='SchemeSplitMALM') 
#