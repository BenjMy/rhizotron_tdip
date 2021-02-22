# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 11:36:11 2021

@author: Benjamin
"""

# load observations file
# compute differences

dates = []

for i, d in enumerate(dates):
    obs = np.loadtxt('process_data/' + d + '')
    pg.load('OMALMIP_0218.txt')
    OBS[:,i] = 
    
    
sc=ax.scatter(coordObs[:,1], coordObs[:,2], c=Obs, 
              cmap ='coolwarm',s=5e2, vmin=vmin, vmax=vmax) # norm=matplotlib.colors.Normalize()
cbar = plt.colorbar(sc,ax=ax)
cbar.set_label('V')   
ax.set_ylabel('y [m]',fontsize=15)
ax.set_xlabel('x [m]',fontsize=15)