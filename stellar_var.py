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
import os
from astropy.stats import sigma_clip

import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 2
# mpl.rcParams['axes.labelsize'] = 14
# mpl.rcParams['legend.fontsize'] = 14
mpl.rcParams['font.size'] = 14
mpl.rcParams.keys()

# =============================================================================
#                       PREPARE DATA
# =============================================================================
# Load LC
path = '/home/dario/AstronomyLeiden/FRP/MRP-ASAS-Rings/zv2/'
file_list  = os.listdir(path)
# rand = np.random.randint(0,2e3)
rand = 765

# Select file
file = file_list[rand]

# Load data
df = pd.read_csv(os.path.join(path,file), sep=';')
# Sort data by date since it is mixed
df = df.sort_values('jd')
# Assign values
time, flux, errflux = df['jd'], df['flux'], df['flux_err']

# # remove nan's
mask = ~np.isnan(flux)
time  = time[mask]
flux  = flux[mask]

# normalise flux
flux  = flux / np.median(flux)
errflux = errflux / np.median(flux)



# sigma_clip 
sigma = 3

mask = sv.sigma_clip_lc(flux, errflux, sigma_lower=3, sigma_upper=1, maxiters=3)
time_sv, flux_sv, errflux_sv = time[~mask], flux[~mask], errflux[~mask]
time_clip, flux_clip, errflux_clip = time[mask], flux[mask], errflux[mask]


## Visualize the distribution of flux error bars
fig, (ax1) = plt.subplots(1,1, figsize=(9,6))
#----------------------------------------------
ax1.hist(errflux_sv, bins='auto', color='green', alpha=0.8, label='Clipped')
ax1.hist(errflux, bins='auto', color='red', alpha=0.4, label='Raw')
ax1.set_ylabel('Number of points')
ax1.set_xlabel('Flux Errorbar [-]')
ax1.legend()
plt.show()
#%%
fig, ax = plt.subplots(1,1,figsize=(20,6))

ax.errorbar(time_sv, flux_sv, errflux_sv, fmt='.k', alpha=0.8, label='Error-clipped')
ax.errorbar(time_clip, flux_clip, errflux_clip, fmt='xr', alpha=0.8, label='Clipped')
ax.set_ylabel('Normalised Flux')
ax.set_title('ASAS_ID = ' + file[2:-4])
ax.legend()
ax.set_xlabel('Time [JD]')
ax.axhspan(np.median(flux_sv)-np.std(flux_sv), np.median(flux_sv)+np.std(flux_sv), alpha=0.4, color='lightgreen')


plt.show()

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
time_sv = np.array(time_sv)
flux_sv = np.array(flux_sv)
errflux_sv = np.array(errflux_sv)
line_flux = np.array(line_flux)

ls1_flux = flux_sv - line_flux

from astropy.timeseries import LombScargle

# freq, power = LombScargle(time, flux, errflux).autopower()
# freq_sc, power_sc = LombScargle(time_sv, flux_sv, errflux_sv).autopower()

ls_p1b, freq, power, _ = sv.lombscargle_periodogram(time, flux, errflux, 
                                                    dt=time[0], plot=False)
ls_p1b_sv, freq_sv, power_sv, _ = sv.lombscargle_periodogram(time_sv, ls1_flux, errflux_sv, 
                                                          dt=time_sv[0], plot=False)


# best_frequency = freq[np.argmax(power)]
# best_frequency_sv = freq_sv[np.argmax(power_sv)]

fig, ax = plt.subplots(1,1,figsize=(9,6))

ax.plot(freq_sv, power_sv, color='green', alpha=0.8, label='Clipped')
ax.axhline(max(power_sv), label='Clipped max(power)', color='green', ls='--', alpha=0.4)
ax.plot(freq, power, alpha=0.5, label='Raw', color='red')
ax.axhline(max(power), label='Raw max(power)', color='red', ls='--', alpha=0.4)

ax.set_ylim(0,None)
ax.set_xlim(freq[0], freq[-1])
ax.set_ylabel('Power')
ax.set_title('ASAS_ID = ' + file[2:-4])
ax.legend()

ax.set_xlabel('Frequency')

plt.show()
# new flux (removed linear trend) for lomb_scargle
# ORIGINAL CODE below
# ls_p1b, ls2_flux = sv.lombscargle_periodogram(time_sv, ls1_flux, errflux_sv, dt=time_sv[0])
# ls_p2b, ls3_flux = sv.lombscargle_periodogram(time_sv, ls2_flux, errflux_sv, dt=time_sv[0])
# ls_p3b, ls4_flux = sv.lombscargle_periodogram(time_sv, ls3_flux, errflux_sv, dt=time_sv[0])
# ls_p4b, ls5_flux = sv.lombscargle_periodogram(time_sv, ls4_flux, errflux_sv, peak_ind=1, dt=time_sv[0])
# ls_p5b, ls6_flux = sv.lombscargle_periodogram(time_sv, ls5_flux, errflux_sv, dt=time_sv[0])


#%%
amplitude, period, phase = ls_p1b_sv
# define the model
clean_flux = flux -  sv.line(time, *p_line, time[0]) - (amplitude * np.sin(2 * np.pi * time / period + phase))

# _, = mcmc.plot_models(time_sv, flux_sv, errflux_sv, [stellar_variation], [Pb_ls], ['stellar variation'])

#%%
# delta_flux = flux_sv - stellar_variation(Pb_ls, time_sv) + 1
# delta_flux = delta_flux / np.median(delta_flux)
# errflux = errflux / np.median(delta_flux)
# plt.scatter(time, delta_flux)

fig, ax = plt.subplots(2,1,figsize=(12,6))

ax[0].errorbar(time, flux, errflux, fmt='.k', alpha=0.8)
ax[0].set_ylabel('Normalised Flux')
ax[0].set_title('ASAS_ID = ' + file[2:-4])
ax[0].axhline(y=1)

plt.setp(ax, ylim=ax[0].get_ylim()) # set the ylim for all other subplots as the first one
label='SV sub'
ax[1].errorbar(time, clean_flux+1, errflux, fmt='.k', alpha=0.8, label=label)
ax[1].set_ylabel('Normalised Flux')
ax[1].set_xlabel('Time [JD]')
ax[1].axhline(y=1)
ax[1].legend()

plt.show()



fig, ax = plt.subplots(1,1,figsize=(20,6))
ax.errorbar(time, clean_flux+1, errflux, fmt='.k', alpha=0.8, label='Clean flux')
ax.set_ylabel('Normalised Flux')
ax.set_title('ASAS_ID = ' + file[2:-4])
ax.legend()
ax.set_xlabel('Time [JD]')
# ax.axhspan(np.median(flux_sv)-np.std(flux_sv), np.median(flux_sv)+np.std(flux_sv), alpha=0.4, color='lightgreen')


plt.show()