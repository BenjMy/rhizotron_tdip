# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 11:36:11 2021

@author: Benjamin
"""
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_CC_violin(m0,tau,fit,r):

    #m0,tau,fit = IPcurves_f.fitDecays(show=True)
    #IPcurves_f.showDecay(nr=np.arange(0,len(IPcurves_f.data['a'])), showFit=False, 
    #                   yscale='linear',xscale='linear')
    #m0_02,tau_02,fit_02 = IPcurves_f.fitDecays(show=True, tmin=0.2)
    #IPcurves_f.showDecay(nr=np.arange(0,len(IPcurves_f.data['a'])), showFit=False, 
    #                   yscale='linear',xscale='linear')

    fig, ax = plt.subplots(nrows=1, ncols=4)
    ax[0].violinplot(m0)
    ax[1].violinplot(tau)
    ax[2].violinplot(fit)
    ax[3].violinplot(r) # abs(IPcurves.data('r')))
    # ax.set_title('Default violin plot')
    ax[0].set_ylabel('Observed values')
    ax[0].set_title('m0 (mV/V)')
    ax[1].set_title('Tau (s)')
    ax[2].set_title('fit rms (log)')
    ax[3].set_title('rhoa')
    
