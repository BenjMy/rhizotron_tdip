# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 11:36:11 2021

@author: Benjamin
"""
import matplotlib.pyplot as plt
import numpy as np

def plot_CC_violin(IPcurves,m0_lim=[0,100],tau_lim=[0.1,10],r_lim=[30,3000], 
                   filtCC=False, fit=False, useColeCole=False):
    '''
    Parameters
    ----------
    IPcurves : TYPE
        DESCRIPTION.
    m0_lim : TYPE, optional
        DESCRIPTION. The default is [0,100].
    tau_lim : TYPE, optional
        DESCRIPTION. The default is [0.1,10].
    r_lim : TYPE, optional
        DESCRIPTION. The default is [30,3000].

    Returns
    -------
    valid_r : TYPE
        DESCRIPTION.
    valid_m : TYPE
        DESCRIPTION.
    valid_tau : TYPE
        DESCRIPTION.
    IPcurves_abs : TYPE
        DESCRIPTION.
    plt : TYPE
        DESCRIPTION.

    '''

    #IPcurves
    #IPcurves_abs.MA = abs(IPcurves.MA)
    #IPcurves_abs.fitDecays(show=False)     
    
    if np.unique(IPcurves.data('u').array())==0:
        IPcurves.data.set('r', IPcurves.data('rhoa')/IPcurves.data('k'))
    else:
        IPcurves.data.set('r', IPcurves.data('u')/IPcurves.data('i'))
                
    if fit:
        if useColeCole:
            m0, tau, c, fit = IPcurves.fitModelDecays(useColeCole=useColeCole)
            IPcurves.data.set('m0', m0)
            IPcurves.data.set('tau', tau)
            IPcurves.data.set('c', c)
            IPcurves.data.set('fit', fit)
            IPcurves.data.save('fitted.data')

            fitCC = [IPcurves.data('r'), m0, tau, c, fit]
            IPcurves.showDecay(nr=np.arange(0,len(IPcurves.data['a'])), showFit=True, 
                          yscale='linear',xscale='linear')
            
            if filtCC:
                print('apply fiter on CC values')
                validCC = IPcurves.filter(rmin=r_lim[0], rmax=r_lim[1],
                                           m0min=m0_lim[0], m0max=m0_lim[1],
                                           taumin=tau_lim[0], taumax=tau_lim[1])
                plt_raw = plt_violin(IPcurves,legend='raw CC range values')

                IPcurves.data.set('valid',validCC)
                IPcurves.MA = IPcurves.MA[:, IPcurves.data['valid'].array()==1]
                IPcurves.data.removeInvalid()
                plt_filter = plt_violin(IPcurves,legend='filtered CC range values')

                return validCC, fitCC, plt_raw, plt_filter
            else:
                plt = plt_violin(IPcurves)
                return fitCC, plt
        else:
            m0, tau,fit = IPcurves.fitDecays()
            IPcurves.data.set('m0', m0)
            IPcurves.data.set('tau', tau)
            IPcurves.data.set('fit', fit)
            IPcurves.data.save('fitted.data')
        
            fitCC = [IPcurves.data('r'), m0, tau, fit]          

            if filtCC:
                print('apply fiter on CC values')
                validCC = IPcurves.filter(rmin=r_lim[0], rmax=r_lim[1],
                                           m0min=m0_lim[0], m0max=m0_lim[1],
                                           taumin=tau_lim[0], taumax=tau_lim[1])
                plt_raw = plt_violin(IPcurves,legend='raw CC range values')
                IPcurves.data.set('valid',validCC)
                IPcurves.MA = IPcurves.MA[:, IPcurves.data['valid'].array()==1]
                IPcurves.data.removeInvalid()
                plt_filter = plt_violin(IPcurves,legend='filtered CC range values')
    
                return validCC, fitCC, plt_raw, plt_filter
        
            else:
                plt = plt_violin(IPcurves)
                return fitCC, plt

    else: 
        plt = plt_violin(IPcurves)
        return plt

    

    
def plt_violin(IPcurves,legend='values'):
    '''
    

    Returns
    -------
    None.

    '''
    try:
        IPcurves.data['c']
        print('plot with c')
        fig, ax = plt.subplots(nrows=1, ncols=5)
        ax[0].violinplot(IPcurves.data['m0'])    
        ax[1].violinplot(IPcurves.data['tau'])
        ax[2].violinplot(IPcurves.data['c'])
        ax[3].violinplot(IPcurves.data['fit'])
    
        ax[4].violinplot(IPcurves.data['r'])
   
        ax[0].set_ylabel(legend)
        
        ax[0].set_title('fit m0 [mV/V]')
        ax[1].set_title('fit tau [s]')
        ax[2].set_title('fit $c$ [-]')
        ax[3].set_title('fit rms (log)')
        ax[4].set_title('u/i (V)')
    except: 
        fig, ax = plt.subplots(nrows=1, ncols=4)
        ax[0].violinplot(IPcurves.data['m0'])
        ax[1].violinplot(IPcurves.data['tau'])
        ax[2].violinplot(IPcurves.data['fit'])
    
        ax[3].violinplot(IPcurves.data['r'])
   
        ax[0].set_ylabel(legend)
        
        ax[0].set_title('fit m0 [mV/V]')
        ax[1].set_title('fit tau [s]')
        ax[2].set_title('fit rms (log)')
        ax[3].set_title('u/i (V)')
    
    return plt
        
    