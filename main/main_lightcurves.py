"""
Created on Tue Jun  8 12:47:36 2021
@author: dario

Title: Main script to read and generate light curve objects
 and upload them to zooniverse
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import os, glob
from astropy.stats import sigma_clip
import matplotlib as mpl
from astropy.io import ascii 
import lightkurve as lk
import matplotlib.gridspec as gridspec
import panoptes as pan
import plot_functions as plot
import support_functions as supp
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
lcs = supp.LightCurveSet(basepath_V, 'dat', 'v')
file = lcs.files()[-4]
df = lcs.data(file)

#%%
outlier_count = []
    
# file = lcs.files()[3829]
# for file in lcs.files()[:100]:
    
def LC_pro(file, show=False):    
    print('Reading {:}'.format(file))
    df = lcs.data(file)
    
    # Fix bad values in 'mag' column
    mag = np.array(df['mag'], dtype='str')
    
    for i,point in enumerate(mag):
        split = point.split('>')
        if len(split) > 1:
            mag[i] = np.nan
    df['mag'] = mag
   
    df['mag'] = df['mag'].astype('float')
    df = df.dropna()
    # name = str(int(df['gaia_id'].unique())) --->>> WRONG last 3 digits
    name = file.split('/')[-1][:-4]
    # mean_mag = float(df['mean_mag'].unique())
    mean_mag = float(df['mag'].mean())
    
    outliers = len(df['mag'][abs(df['mag'] - mean_mag) > 0.5])
    print('Outliers {:}'.format(outliers))
    if outliers > 10:
        outlier_count.append(mean_mag)

    # BY DEFAULT
    skip = False
    
    
    if outliers < 10:
        skip = True
    if mean_mag > 15.95:
        skip = True
    if len(df) < 30:
        skip = True
        
    
    if skip == True:
        print('------------!!!!!!!!!!!!!!----------')
        print('Skipping {:}'.format(file))
        print('------------!!!!!!!!!!!!!!----------')
        
    if skip == False:
        time, flux, errflux = df['HJD'].to_numpy(), df['flux'].to_numpy(), df['flux_err'].to_numpy()
        
        
        
        ## Clip the data until the stdeviation changes less than eps
        time_sv, flux_sv, errflux_sv, count = supp.iterative_clipping(time, flux, errflux, eps=0.01)
        
        # identify the indices of the removed points to plot them too
        intersect, ind1, ind2 = np.intersect1d(time, time_sv, return_indices=True)
        
        time_clip = np.delete(np.array(time), ind1)
        flux_clip = np.delete(np.array(flux), ind1)
        errflux_clip = np.delete(np.array(errflux), ind1)
        
        
        ## Plot the lightcurve with the clipped points (red)
        lc_sv = np.array([time_sv, flux_sv, errflux_sv])
        clipped_data = np.array([time_clip, flux_clip, errflux_clip])
        # _ = plot.plot_lc_clipped(lc_sv, clipped_data, file, size=(20,6))
        
        
        #%%
        # =============================================================================
        #                           Lomb-Scargle Periodogram
        # =============================================================================
        
        # transform to numpy arrays for compatibility (pandas KeyError) 
        time_sv = np.array(time_sv)
        flux_sv = np.array(flux_sv)
        errflux_sv = np.array(errflux_sv)
        
        from astropy.timeseries import LombScargle
        
        # this one is to plot the raw periodogram
        ls_p1b, freq, power, _ = supp.lombscargle_periodogram(time, flux-1.0, errflux, dt=time[0], plot=False)
        
        
        
        ls_p1b_sv, freq_sv, power_sv, ls2_flux = supp.lombscargle_periodogram(time_sv, flux_sv, errflux_sv, dt=time_sv[0], plot=False)
        #%%
        # Gather parameters
        t0, amplitude, period, phase = ls_p1b_sv
        
        
        fig = plt.figure(figsize=(14,8))
        gs = fig.add_gridspec(2, 2)
        
        
        ax1 = fig.add_subplot(gs[0,:])
        ax1.invert_yaxis()
        label='V band'
        # ax1.errorbar(time, clean_flux+1, errflux, fmt='.k', alpha=0.8, label=label)
        time_corr = time - 2450000 # transform HDJ date
        # ax1.errorbar(time_corr, flux, errflux, fmt='.k', alpha=0.6, label=label)
        ax1.plot(time_corr, df['mag'],'.k', alpha=0.6, label=label)

        ax1.set_title('GAIA_ID = ' + name)
        # ax1.set_ylabel('Relative flux')
        ax1.set_ylabel('Magnitude')
        ax1.set_xlabel('HJD - 2450000')
        mean_mag_str = 'Mean V-mag = {:.1f}'.format(mean_mag)
        
        if mean_mag > 17:
            color_mag = 'red'
        elif (mean_mag < 17) and (mean_mag > 15):
            color_mag = 'darkorange'
        else:
            color_mag = 'darkgreen'
        ax1.text(x=0.02, y=0.89, s=mean_mag_str, transform=ax1.transAxes, color=color_mag)
        bars = [0, int(np.round(len(df) / 2.)), -1]
        ax1.errorbar(time_corr[bars], df['mag'].values[bars], df['mag_err'].values[bars], 
                     fmt='.', color=color_mag, alpha=0.3)
        
        ymin_auto, ymax_auto = ax1.yaxis.get_data_interval()
        # ax1.set_ylim(ymin=np.max([ymin_auto, 0.0]))
        # ax1.set_ylim(ymax=np.min([ymax_auto, 2.0]))
        
        
        ax2 = fig.add_subplot(gs[1,0])
        # _ = plot.hist_errbars(errflux_sv, errflux, count, ax=ax2)
        lc = lk.LightCurve(time, flux, errflux)
        lc_folded = lc.fold(period=period)
        _ = lc_folded.plot(ls='', marker='.', markersize=4, color='k', ax=ax2, label='Folded')
        period_str = 'period : {:.2f} d'.format(period)
        ax2.text(x=0.02, y=0.89, s=period_str, transform=ax2.transAxes)
        ax2.legend(loc='upper right')
        ymin_auto, ymax_auto = ax2.yaxis.get_data_interval()
        ax2.set_ylim(ymin=np.max([ymin_auto, 0.0]))
        ax2.set_ylim(ymax=np.min([ymax_auto, 2.0]))
        
        ax3 = fig.add_subplot(gs[1,1])
        data = [freq_sv, power_sv, freq, power]
        # data = [freq_sv, power_sv]
        _ = plot.plot_periodogram(data, ax=ax3, clipped=False)
        xticks = np.arange(0,np.max(freq_sv),1)
        ax3.set_xticks(xticks)
        
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
# file = '/home/dario/AstronomyLeiden/FRP/leiden_vband/camfix/3005888048142801664.dat'
# file =  '/home/dario/AstronomyLeiden/FRP/leiden_vband/camfix/3217523278911100928.dat'
# file = '/home/dario/AstronomyLeiden/FRP/leiden_vband/camfix/274979577422304000.dat'

# test =  LC_pro(file, show=True)
# inspecting individual LCs
#%%
'''
df = lcs.data(file)
df

# Fix bad values in 'mag' column
mag = np.array(df['mag'], dtype='str')

for i,point in enumerate(mag):
    split = point.split('>')
    if len(split) > 1:
        mag[i] = np.nan
df['mag'] = mag
df['mag'] = df['mag'].astype('float')

df = df.dropna()
df.head()
plot.plot_lc(df, y='mag', name=str(int(df['gaia_id'][0])))
plt.show()

df

time, flux, errflux = df['HJD'], df['flux'], df['flux_err']

## Clip the data until the stdeviation changes less than eps
time_sv, flux_sv, errflux_sv, count = supp.iterative_clipping(time, flux, errflux, eps=0.01)
time_sv = np.array(time_sv)
flux_sv = np.array(flux_sv)
errflux_sv = np.array(errflux_sv)

ls_p1b_sv, freq_sv, power_sv, ls2_flux = supp.lombscargle_periodogram(time_sv, flux_sv, errflux_sv, dt=time_sv[0], plot=False)

# def flat_filter(df, eps=0.25):

eps = 0.25
mean_mag = float(df['mag'].mean())
print(len(df['mag'][abs(df['mag'] - mean_mag) > 0.25]))
    
len(df['mag'])
print(int(np.round(len(df) / 2)))

df['mag'][[0,1,2]]
#%%
time_corr = df['HJD'] - 2450000 # transform HDJ date
y='mag'
y_err = 'mag_err'
bars = [0, int(np.round(len(df)/2.)), -1]

fig, ax = plt.subplots()

ax.plot(time_corr, df[y], '.k')
ax.errorbar(time_corr.values[bars], df[y].values[bars], df[y_err].values[bars], fmt='.k', alpha=0.6)

plt.show()
#%%%%
## Individual query 
# gaia_ids = np.array(lcs.files(gaia_id=True))
# filelist = lcs.files()

# gaia_id = 4133345550648167680
# index = int(np.argwhere(gaia_ids == gaia_id))
# file = filelist[index]

# test =  LC_pro(file, show=True)


###################################
'''
# Large query
num = len(lcs.files())
filelist = lcs.files()

outlist = []
for file in filelist:
    outpath = LC_pro(file, show=False)
    outlist.append(outpath)
    
print('{:} LCs out of {:} with 10 or more points beyond 0.5 from mean_mag'.format(len(outlier_count), num))
print('{:} \%'.format(len(outlier_count)*100 / num))


   
outlist = list(filter(None, outlist))

_ = pan.panoptes_add(outlist, 'filter_0.5mag_less_than_16mag')









