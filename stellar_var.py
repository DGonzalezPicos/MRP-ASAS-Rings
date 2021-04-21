"""
Created on Wed Mar 17 11:10:39 2021
@author: dario

Title: Removing Stellar Variations
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
import sv_util as sv
from Code import mcmc
import os, glob
from astropy.stats import sigma_clip
import matplotlib as mpl
from astropy.io import ascii 
mpl.rcParams['lines.linewidth'] = 2
# mpl.rcParams['axes.labelsize'] = 14
# mpl.rcParams['legend.fontsize'] = 14
mpl.rcParams['font.size'] = 14
mpl.rcParams.keys()

# =============================================================================
#                       PREPARE DATA
# =============================================================================
# Load LC

basepath = '/home/dario/AstronomyLeiden/FRP/leiden_vband/lc/camfix/'
basepath_g = '/home/dario/AstronomyLeiden/FRP/gband/'
#%%
basepath_gband = '/home/dario/AstronomyLeiden/FRP/MRP-ASAS-Rings/zv2'
os.chdir(basepath_gband)
file_list_gband = glob.glob(basepath_gband + '/*.csv')
file = file_list_gband[42]

# print('Converting csv to dat...')
# for file in file_list_gband:
#     df = pd.read_csv(file, sep=';')
#     df = df.drop(columns='cam')
#     df = df.drop(columns=['Unnamed: 0','limit'])
#     df = df.dropna()
#     hdr = 'HJD	flux	flux_err mag    mag_err	FWHM ASAS-SN_ID'
#     file_name = basepath_g + file.split('/')[-1][2:-4] + '.dat'
#     try:
#         np.savetxt(file_name, df.values, fmt='%.5f', header=hdr, comments='#')
#     except TypeError:
#         print(file)
# print('Finished!')
# df.columns
# file_list_gband_dat = glob.glob(basepath_g + '*.dat')
# file = file_list_gband_dat[2]
# data = ascii.read(file, delimiter=' ')





#%%
file_list = glob.glob(basepath + '*.dat')
file_list_g = glob.glob(basepath_g + '*.dat')

N = len(file_list)
N_g = len(file_list_g)
print('Number of files in the V-directory ... ', N)
print('Number of files in the g-directory ... ', N_g)
'''
epochs_V = []
epochs_g = []
start = 1.0e12
start_g = 1.0e12
end  = 0.0
end_g = 0.0
file = file_list[3]
data = ascii.read(file, delimiter='\t')
data.colnames
data['HJD']

file = file_list_g[11853]
data = np.loadtxt(file)
# len(data)
# if len(data) < 1:
#     print('empty file... ', file)
    

def span(time, start, end):
    if len(time) > 1:
        if time[0] < start:
            start = time[0]
        if time[-1] > end:
            end = time[-1]
    return start, end

for file, file_g in zip(file_list, file_list_g):
    time, flux, errflux, name = sv.load_dat(file, y='flux')
    time_g, flux_g, errflux_g, name_g = sv.read_gband(file_g)
    epochs_V.append(len(time))
    epochs_g.append(len(time_g))
    
    start, end = span(time, start, end)
    start_g, end_g = span(time_g, start_g, end_g)

            
            
    
print(np.mean(epochs_V))
print(np.mean(epochs_g))
span = end - start
print('span: {:.2f} years ({:.2f} days) '.format(span/365, span))


fig, ax = plt.subplots(1)

def epochs_hist(epochs, band):
    ax.hist(epochs, bins='auto', color='green', alpha=0.8, label='V-band')
    span_str = ('span = {:.2f} years'.format(span/365,1))
    start_str = ('start = {:.2f} HJD'.format(start))
    end_str = ('end = {:.2f} HJD'.format(end))
    ax.text(x=0.5, y=0.75, s=span_str, transform=ax.transAxes)
    ax.text(x=0.5, y=0.65, s=start_str, transform=ax.transAxes)
    ax.text(x=0.5, y=0.55, s=end_str, transform=ax.transAxes)
    
    # ax.set_ylabel('Number of points')
    ax.set_xlabel('Data points')
    ax.legend()
    return ax

# plt.savefig('/home/dario/AstronomyLeiden/FRP/hist_epochs.png', dpi=100)
plt.show()
 
#%%   

set(file_list_g, file_list)
'''
#%%    

band = 'V' 
if band == 'V':
    file = file_list[448]
    time, flux, errflux, name = sv.load_dat(file, y='flux')

    
if band == 'g':
    file = file_list_g[333]
    time, flux, errflux, name = sv.read_gband(file)

# path = '/home/dario/AstronomyLeiden/FRP/MRP-ASAS-Rings/zv2/'
# file_list  = os.listdir(path)
# # rand = np.random.randint(0,2e3)
# # rand = 909


# # Select file
# file = file_list[rand]

# # Load data
# df = pd.read_csv(os.path.join(path,file), sep=';')
# # Sort data by date since it is mixed
# df = df.sort_values('jd')
# # Assign values
# time, flux, errflux = df['jd'], df['flux'], df['flux_err']

# # remove nan's
# mask = ~np.isnan(flux)
# time  = time[mask]
# flux  = flux[mask]

# normalise flux
flux  = flux / np.median(flux)
errflux = errflux / np.median(flux)

## Clip the data until the stdeviation changes less than eps
time_sv, flux_sv, errflux_sv, count = sv.iterative_clipping(time, flux, errflux, eps=0.001)

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
# =============================================================================
linear_trend = lambda P, time: sv.line(time, *P, time[0])
line_curvefit = lambda time, m, c: sv.line(time, m, c, time[0])
p_line, cov = curve_fit(line_curvefit, time_sv, flux_sv, sigma=errflux_sv)
line_flux = sv.line(time_sv, *p_line, time_sv[0])

print('---------Linear Fit parameters------ \
\n m= {:2.2} +/- {:2.2} \n \
b= {:.2} +/- {:.2}'.format(p_line[0], np.sqrt(cov[0][0]),
                                  p_line[1], np.sqrt(abs(cov[1][1]))))
    
# mcmc.plot_models(time_sv, flux_sv, errflux_sv, [linear_trend], [p_line], ['linear'])

#%%
# =============================================================================
#                           Lomb-Scargle Periodogram
# =============================================================================

# transform to numpy arrays for compatibility (pandas KeyError) 
time = np.array(time)
time_sv = np.array(time_sv)
flux_sv = np.array(flux_sv)
errflux_sv = np.array(errflux_sv)
line_flux = np.array(line_flux)

ls1_flux = flux_sv - line_flux

from astropy.timeseries import LombScargle

# this one is just to compare periodograms clipped vs. raw
ls_p1b, freq, power, _ = sv.lombscargle_periodogram(time, flux-1.0, errflux, 
                                                    dt=time[0], plot=False)



ls_p1b_sv, freq_sv, power_sv, ls2_flux = sv.lombscargle_periodogram(time_sv, ls1_flux, errflux_sv, dt=time_sv[0], plot=True)

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
#%% Implementation with the LIGHTKURVE package

# mpl.rcParams['lines.color']='black'
'''
import lightkurve as lk
lc = lk.LightCurve(time, flux, errflux)
np.std(lc.flux)

clc = lc.flatten()
np.std(clc.flux)

pdg = lc.to_periodogram(minimum_period=2.5, maximum_period=20)
pdg.show_properties()
pdg.default_view='period'
max_freq = pdg.frequency_at_max_power
max_period = 1/max_freq
t0, t1, t2 = pdg._LS_object.model_parameters(max_freq)


time_y = np.linspace(time[0], time[-1], 10000)
y = t0 + (t1 * np.sin(2*np.pi * max_freq.value * time_y)) + (t2 * np.cos(2*np.pi * max_freq.value * time_y))


#-----------------------------------------------------------------------------
fig, (ax1, ax2, ax3) = plt.subplots(3,1, figsize=(20,12))

_ = lc.errorbar(ax=ax1, label='Raw')
ax1.plot(time_y, y+1, lw=1)

_ = clc.errorbar(ax=ax2, label='Flatten')

_ = lc.plot_river(max_period.value)
_ = pdg.plot(ax=ax3, label='Periodogram')
plt.show()
'''
#%%
# Gather parameters
t0, amplitude, period, phase = ls_p1b_sv
# define the model
clean_flux = flux -  sv.line(time, *p_line, time[0]) - (amplitude * np.sin(((2 * np.pi * time)/period) + phase))

# _, = mcmc.plot_models(time_sv, flux_sv, errflux_sv, [stellar_variation], [Pb_ls], ['stellar variation'])

#%%
# delta_flux = flux_sv - stellar_variation(Pb_ls, time_sv) + 1
# delta_flux = delta_flux / np.median(delta_flux)
# errflux = errflux / np.median(delta_flux)
# plt.scatter(time, delta_flux)

fig = plt.figure(figsize=(16,8))

ax1 = fig.add_subplot(3,1,1)
_ = sv.plot_lc_clipped(lc_sv, clipped_data, file, ax=ax1)
ax1.set_title('ASAS_ID = ' + name + ' ('+ band+ '-band)')
# ax1 = plt.plot(model_time, model, '-b')


ax2 = fig.add_subplot(3,1,2)
label='SV sub'
ax2.errorbar(time, clean_flux+1, errflux, fmt='.k', alpha=0.8, label=label)


ax3 = fig.add_subplot(3,2,5)
_ = sv.hist_errbars(errflux_sv, errflux, count, ax=ax3)

ax4 = fig.add_subplot(3,2,6)
data = [freq_sv, power_sv, freq, power]
_ = sv.plot_periodogram(data, ax=ax4)
# fig, ax = plt.subplots(2,1,figsize=(12,6))

# ax[0].errorbar(time, flux, errflux, fmt='.k', alpha=0.8)
# ax[0].set_ylabel('Normalised Flux')
# ax[0].set_title('ASAS_ID = ' + file[2:-4])
# ax[0].axhline(y=1)

# plt.setp(ax, ylim=ax[0].get_ylim()) # set the ylim for all other subplots as the first one
# label='SV sub'
# ax[1].set_ylabel('Normalised Flux')
# ax[1].set_xlabel('Time [JD]')
# ax[1].axhline(y=1)
# ax[1].legend()

# plt.show()



# fig, ax = plt.subplots(1,1,figsize=(20,6))
# ax.errorbar(time, clean_flux+1, errflux, fmt='.k', alpha=0.8, label='Clean flux')
# ax.set_ylabel('Normalised Flux')
# ax.set_title('ASAS_ID = ' + file[2:-4])
# ax.legend()
# ax.set_xlabel('Time [JD]')
# # ax.axhspan(np.median(flux_sv)-np.std(flux_sv), np.median(flux_sv)+np.std(flux_sv), alpha=0.4, color='lightgreen')


# plt.show()