# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 11:36:11 2021

@author: Benjamin
"""
import matplotlib.pyplot as plt
import numpy as np

def plot_CC_violin(IPcurves,m0_lim=[0,100],tau_lim=[0.1,10],r_lim=[30,3000], 
                   filtCC=False, fit=False):
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
    
    if fit:
        m0, tau, c, fit = IPcurves.fitModelDecays(useColeCole=True)
        IPcurves.data.set('m0', m0)
        IPcurves.data.set('tau', tau)
        IPcurves.data.set('c', c)
        IPcurves.data.set('fit', fit)
        IPcurves.data.save('fitted.data')
    
        IPcurves.data.set('r', IPcurves.data('u')/IPcurves.data('i'))

        fitCC = [IPcurves.data('r'), m0, tau, c, fit]
        IPcurves.showDecay(nr=np.arange(0,len(IPcurves.data['a'])), showFit=True, 
                      yscale='linear',xscale='linear')
        
        if filtCC:
            print('apply fiter on CC values')
            validCC = IPcurves.filter(rmin=r_lim[0], rmax=r_lim[1],
                                       m0min=m0_lim[0], m0max=m0_lim[1],
                                       taumin=tau_lim[0], taumax=tau_lim[1])  
    
            # valid_r = IPcurves.filter(rmin=r_lim[0], rmax=r_lim[1]).array()  # apparent resistivity
            
            # print('m0min:' + str(min(IPcurves.data('m0'))))
            # print('m0max:'+ str(max(IPcurves.data('m0'))))
            # print('m0_lim:'+ str(m0_lim[0]) + '/' + str(m0_lim[1]))
            # valid_m = IPcurves.filter(m0min=m0_lim[0], m0max=m0_lim[1]).array()  # apparent resistivity
            # print('valid_m:'+ str(np.count_nonzero(valid_m==1)))
            # print(np.vstack((IPcurves.data['m0'].array(), valid_m)).T)
    
            # print('taumin:' + str(min(IPcurves.data('tau'))))
            # print('taumax:'+ str(max(IPcurves.data('tau'))))
            # print('tau_lim:'+ str(tau_lim[0]) + '/' + str(tau_lim[1]))
            # valid_tau = IPcurves.filter(taumin=tau_lim[0], taumax=tau_lim[1]).array()  # fitted relaxation time
            # print('valid_tau:'+ str(np.count_nonzero(valid_tau==1)))
            
            # validCC = [valid_r,valid_m,valid_tau]
            #IPcurves.data.set('validCC', fit)
            #IPcurves.data.save('fitted.data')
            plt = plt_violin(IPcurves,legend='filtered CC range values')

            return validCC, fitCC, plt
    
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
    fig, ax = plt.subplots(nrows=1, ncols=5)
    ax[0].violinplot(IPcurves.data['m0'])    
    ax[1].violinplot(IPcurves.data['tau'])
    ax[2].violinplot(IPcurves.data['c'])
    ax[3].violinplot(IPcurves.data['fit'])
    
    ax[4].violinplot(IPcurves.data['r'])
   
    ax[0].set_ylabel(legend)
    
    ax[0].set_title('fit m0 [mV/V]')
    ax[1].set_title('fit $\tau$ [s]')
    ax[2].set_title('fit $c$ [-]')
    ax[3].set_title('fit rms (log)')
    ax[4].set_title('u/i (V)')
    
    return plt
        
    