# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 11:36:11 2021

@author: Benjamin
"""
import matplotlib.pyplot as plt
import numpy as np

def plot_CC_violin(IPcurves,m0=[0,100],tau=[0.1,10],r=[30,3000], 
                   filtCC=False, fit=False):
    '''
    Parameters
    ----------
    IPcurves : TYPE
        DESCRIPTION.
    m0 : TYPE, optional
        DESCRIPTION. The default is [0,100].
    tau : TYPE, optional
        DESCRIPTION. The default is [0.1,10].
    r : TYPE, optional
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
    if fit:
        m0, tau, c, fit = IPcurves.fitModelDecays(useColeCole=True)
        IPcurves.data.set('m0', m0)
        IPcurves.data.set('tau', tau)
        IPcurves.data.set('c', c)
        IPcurves.data.set('fit', fit)

    
    IPcurves.showDecay(nr=np.arange(0,len(IPcurves_abs.data['a'])), showFit=True, 
                      yscale='linear',xscale='linear')

    fig, ax = plt.subplots(nrows=1, ncols=4)
    
    ax[0].violinplot(IPcurves.data['m0'])    
    ax[1].violinplot(IPcurves.data['tau'])
    ax[2].violinplot(IPcurves.data['fit'])
    
    
    IPcurves.data.set('r', IPcurves.data('u')/IPcurves.data('i'))
    ax[3].violinplot(IPcurves.data('u')/IPcurves.data('i')) # abs(IPcurves.data('r')))
    # ax.set_title('Default violin plot')
    ax[0].set_ylabel('Values')
    ax[0].set_title('fit m0 (mV/V)')
    ax[1].set_title('fit Tau (s)')
    ax[2].set_title('fit rms (log)')
    ax[3].set_title('u/i (V)')
    
    IPcurves.data.save('fitted.data')
    valid_r = IPcurves.filter(rmin=r[0], rmax=r[1]).array()  # apparent resistivity
    valid_m = IPcurves.filter(m0min=m0[0], m0max=m0[1]).array()  # apparent resistivity
    valid_tau = IPcurves.filter(taumin=tau[0], taumax=tau[1]).array()  # fitted relaxation time
    
    return valid_r,valid_m, valid_tau, IPcurves_abs, plt
    