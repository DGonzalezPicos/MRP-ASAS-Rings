"""
Created on Wed Mar 17 11:10:39 2021
@author: dario

Title: Removing Stellar Variations
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import sys
sys.path.append('/.../application/app/folder')
import util.util as ut
from Code import mcmc
import os, glob
from astropy.stats import sigma_clip
import matplotlib as mpl
from astropy.io import ascii 
import lightkurve as lk
import matplotlib.gridspec as gridspec
import panoptes as pan
mpl.rcParams['lines.linewidth'] = 2
# mpl.rcParams['axes.labelsize'] = 14
# mpl.rcParams['legend.fontsize'] = 14
mpl.rcParams['font.size'] = 14
mpl.rcParams.keys()
plt.ioff() # disable plt plotting automatically to save CPU
# =============================================================================
#                       PREPARE DATA
# =============================================================================
# Load LC

basepath = '/home/dario/AstronomyLeiden/FRP/merged-bands/'
basepath_V = '/home/dario/AstronomyLeiden/FRP/leiden_vband/camfix/'
lcs = ut.LightCurveSet(basepath_V, 'dat', 'V')
#%%
'''
N = len(lcs.files())
print('Number of files in directory ... ', N)

span_list = []
for file in lcs.files():
    data = lcs.data(file)
    span = max(data['HJD']) - min(data['HJD'])
    span_list.append(span)
            
            
   
span = np.mean(span_list)
print('span: {:.2f} years ({:.2f} days) '.format(span/365, span))

# PLOT HISTOGRAM 
fig, ax = plt.subplots(1, figsize=(12,8))

ax.hist(span_list, bins='auto', alpha=0.8, label='Span')
ax.set_xlabel('Data points')

str_total = 'Number of lightcurves: {:}'.format(N)
# ax.text(x=0.15, y=0.75, s=str_total, transform=ax.transAxes)
ax.set_title(str_total)
ax.legend()
plt.show()
plt.savefig('span_hist_merged.png', dpi=200, bbox_inches='tight')
'''

#%%

    
# file = lcs.files()[3829]
# for file in lcs.files()[:100]:
    
def LC_pro(file, show=False):    
    df = lcs.data(file)

    # name = str(int(df['gaia_id'].unique())) --->>> WRONG last 3 digits
    name = file.split('/')[-1][:-4]
    mean_mag = float(df['mean_mag'].unique())
    
    time, flux, errflux = df['HJD'].to_numpy(), df['flux'].to_numpy(), df['flux_err'].to_numpy()
    
    
    
    ## Clip the data until the stdeviation changes less than eps
    time_sv, flux_sv, errflux_sv, count = sv.iterative_clipping(time, flux, errflux, eps=0.01)
    
    # identify the indices of the removed points to plot them too
    intersect, ind1, ind2 = np.intersect1d(time, time_sv, return_indices=True)
    
    time_clip = np.delete(np.array(time), ind1)
    flux_clip = np.delete(np.array(flux), ind1)
    errflux_clip = np.delete(np.array(errflux), ind1)
    
    ## Visualize the distribution of flux error bars
    # _ = sv.hist_errbars(errflux_sv, errflux, count)
    
    ## Plot the lightcurve with the clipped points (red)
    lc_sv = np.array([time_sv, flux_sv, errflux_sv])
    clipped_data = np.array([time_clip, flux_clip, errflux_clip])
    # _ = sv.plot_lc_clipped(lc_sv, clipped_data, file, size=(20,6))
    
    #%%
    # =============================================================================
    #               REMOVE LINEAR TREND
    # # =============================================================================
    # linear_trend = lambda P, time: sv.line(time, *P, time[0])
    # line_curvefit = lambda time, m, c: sv.line(time, m, c, time[0])
    # p_line, cov = curve_fit(line_curvefit, time_sv, flux_sv, sigma=errflux_sv)
    # line_flux = sv.line(time_sv, *p_line, time_sv[0])
    
    # print('---------Linear Fit parameters------ \
    # \n m= {:2.2} +/- {:2.2} \n \
    # b= {:.2} +/- {:.2}'.format(p_line[0], np.sqrt(cov[0][0]),
    #                                   p_line[1], np.sqrt(abs(cov[1][1]))))
        
    # mcmc.plot_models(time_sv, flux_sv, errflux_sv, [linear_trend], [p_line], ['linear'])
    
    #%%
    # =============================================================================
    #                           Lomb-Scargle Periodogram
    # =============================================================================
    
    # transform to numpy arrays for compatibility (pandas KeyError) 
    time_sv = np.array(time_sv)
    flux_sv = np.array(flux_sv)
    errflux_sv = np.array(errflux_sv)
    # line_flux = np.array(line_flux)
    
    ls1_flux = flux_sv
    
    from astropy.timeseries import LombScargle
    
    # this one is to plot the raw periodogram
    ls_p1b, freq, power, _ = sv.lombscargle_periodogram(time, flux-1.0, errflux, dt=time[0], plot=False)
    
    
    
    ls_p1b_sv, freq_sv, power_sv, ls2_flux = sv.lombscargle_periodogram(time_sv, ls1_flux, errflux_sv, dt=time_sv[0], plot=False)
    
    # Perform the LSP "n" times with the residuals of each iteration
    # def mult_lsp(time, flux, errflux, n=2, **kwargs):
    #     parameters = []
    #     for j in range(n):
    #         p, _, _, flux_out = sv.lombscargle_periodogram(time, flux, errflux, plot=True, dt=time[0])
    #         parameters.append(p)
    #         flux = flux_out
    #     return parameters
    
    # CALL THE FUNCTION
    # parameters = mult_lsp(time_sv, ls2_flux, errflux_sv)  
  
    #%%
    # Gather parameters
    t0, amplitude, period, phase = ls_p1b_sv
    # define the model
    # clean_flux = flux -  sv.line(time, *p_line, time[0]) - (amplitude * np.sin(((2 * np.pi * time)/period) + phase))
    
    # _, = mcmc.plot_models(time_sv, flux_sv, errflux_sv, [stellar_variation], [Pb_ls], ['stellar variation'])
    
    
    fig = plt.figure(figsize=(14,8))
    gs = fig.add_gridspec(2, 2)
    
    
    ax2 = fig.add_subplot(gs[0,:])
    label='V+g band'
    # ax2.errorbar(time, clean_flux+1, errflux, fmt='.k', alpha=0.8, label=label)
    time_corr = time - 2450000 # transform HDJ date
    # ax2.errorbar(time_corr, flux, errflux, fmt='.k', alpha=0.6, label=label)
    ax2.plot(time_corr, flux,'.k', alpha=0.6, label=label)
    ax2.set_title('GAIA_ID = ' + name)
    ax2.set_ylabel('Relative flux')
    ax2.set_xlabel('HJD - 2450000')
    mean_mag_str = 'Mean V-mag = {:.1f}'.format(mean_mag)
    
    if mean_mag > 17:
        color_mag = 'red'
    elif (mean_mag < 17) and (mean_mag > 15):
        color_mag = 'darkorange'
    else:
        color_mag = 'darkgreen'
    ax2.text(x=0.02, y=0.89, s=mean_mag_str, transform=ax2.transAxes, color=color_mag)
    
    ymin_auto, ymax_auto = ax2.yaxis.get_data_interval()
    ax2.set_ylim(ymin=np.max([ymin_auto, 0.0]))
    ax2.set_ylim(ymax=np.min([ymax_auto, 2.0]))
    
    
    ax3 = fig.add_subplot(gs[1,0])
    # _ = sv.hist_errbars(errflux_sv, errflux, count, ax=ax3)
    lc = lk.LightCurve(time, flux, errflux)
    lc_folded = lc.fold(period=period)
    _ = lc_folded.plot(ls='', marker='.', markersize=4, color='k', ax=ax3, label='Folded')
    period_str = 'period : {:.2f} d'.format(period)
    ax3.text(x=0.02, y=0.89, s=period_str, transform=ax3.transAxes)
    ax3.legend(loc='upper right')
    ymin_auto, ymax_auto = ax3.yaxis.get_data_interval()
    ax3.set_ylim(ymin=np.max([ymin_auto, 0.0]))
    ax3.set_ylim(ymax=np.min([ymax_auto, 2.0]))
    
    ax4 = fig.add_subplot(gs[1,1])
    data = [freq_sv, power_sv, freq, power]
    _ = sv.plot_periodogram(data, ax=ax4, clipped=False)
    xticks = np.arange(0,np.max(freq_sv),1)
    ax4.set_xticks(xticks)
    
    # plt.tight_layout()
    if show == True:
        print('GAIA_ID = ' + name)
        plt.show()
        # outpath = '/home/dario/AstronomyLeiden/FRP/zoo/' + 'VARIABLE' +'.png'
        # plt.savefig(outpath, dpi=200, bbox_inches='tight')
        # print('Saved {:}'.format(outpath))
        return name

    if show == False:
        outpath = '/home/dario/AstronomyLeiden/FRP/zoo/' + str(name) +'.png'
        
        plt.savefig(outpath, dpi=200, bbox_inches='tight')
        print('Saved {:}'.format(outpath))
        fig.clear()
        plt.close(fig)
        return outpath


#%%
random = np.random.randint(0, len(lcs.files()))
file = lcs.files()[random]
# file = lcs.files()[13057]
# file = '/home/dario/AstronomyLeiden/FRP/merged-bands/941838052582180864.dat'
name = LC_pro(file, show=True)

# filelist = lcs.files()[:100]

# outlist = []
# for file in filelist:
#     outpath = LC_pro(file, show=False)
#     outlist.append(outpath)
    
# _ = pan.panoptes_add(outlist, '100-no_errorbar')

