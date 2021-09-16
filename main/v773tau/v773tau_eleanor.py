"""
Created on Wed Sep  8 11:47:24 2021
@author: dario

Title: 
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, glob
from astropy.io import ascii 
import sys 
sys.path.append('../')
import support_functions as supp
import plot_functions as plot
import eleanor
from astropy import units as u
from astropy.coordinates import SkyCoord

gaia_id = 163184366130809984
# query the data
# star = eleanor.multi_sectors(gaia=gaia_id, sectors='all')
# no results with id... trying with coordinates
simbad_coords=SkyCoord('04h14m12.9232870093s','+28d12m12.378540314s', frame='icrs')
star = eleanor.multi_sectors(coords=simbad_coords, sectors='all')

def plot_sector_data(sector_data, ax=None):
    ax = ax or plt.gca()
    q = sector_data.quality == 0
    
    # no offset
    # ax.plot(sector_data.time[q], sector_data.raw_flux[q]/np.nanmedian(sector_data.raw_flux[q]), '.k', label='RAW')

    # ax.plot(sector_data.time[q], sector_data.raw_flux[q]/np.nanmedian(sector_data.raw_flux[q])+0.06, '.k', label='RAW')
    ax.plot(sector_data.time[q], sector_data.corr_flux[q]/np.nanmedian(sector_data.corr_flux[q]) + 0.03, '.r', label='CORR')
    ax.plot(sector_data.time[q], sector_data.pca_flux[q]/np.nanmedian(sector_data.pca_flux[q]), '.g', label='PCA')
    # plt.plot(sector_data.time[q], sector_data.psf_flux[q]/np.nanmedian(sector_data.psf_flux[q]) - 0.02, 'b')
    ax.set_ylabel('Normalized Flux')
    ax.set_xlabel('Time [BJD - 2457000]')
    # plt.show()
    return ax

sector_data = eleanor.TargetData(star[0], height=15, width=15, bkg_size=31, do_psf=False, do_pca=True, regressors='corner')

fig, ax = plt.subplots(1, figsize=(9,6))
_ = plot_sector_data(sector_data, ax=ax)
ax.legend()
ax.set_title('V773Tau')
ax.text(s='Sector 19', x=0.05, y=0.05, transform=ax.transAxes, color='orange')
# plt.savefig('eleanor_v773tau.png', dpi=200, bbox_inches='tight')


#%%
from astropy.timeseries import LombScargle

sector_data.corr_flux_norm = sector_data.corr_flux / np.nanmedian(sector_data.corr_flux)

q = sector_data.quality == 0
ls = LombScargle(sector_data.time[q], sector_data.corr_flux_norm[q])
freq, powers = ls.autopower(nyquist_factor=50, minimum_frequency=1/10.,maximum_frequency=1/1.2)

period_days = 1. / freq
period_hours = period_days * 24


best_period = supp.pick_period(freq, powers)
phase = (sector_data.time[q] / best_period) % 1

print("Best period: {:.2f} days".format(best_period))

fig, (ax1, ax2) = plt.subplots(2, figsize=(12,9))
ax1.plot(period_days, powers, '-k')
ax1.axvline(x=best_period, ls='-', lw=8, alpha=0.3, label='{:.2f} d'.format(best_period), color='green')
# ax1.fill_betweenx(x1=best_period-0.2, x2=best_period+0.2, alpha=0.4)
ax1.legend()
ax1.grid()

ax1.set(xlabel='Period (days)', ylabel='Lomb-Scargle Power')
ymin_auto, ymax_auto = ax1.yaxis.get_data_interval()
ax1.set_ylim(ymax=np.max([ymax_auto, 0.10]))



# ax2.errorbar(phase, df.flux, df.flux_err,
#                fmt='.k', ecolor='gray', ms=1, capsize=0, zorder=1)
ax2.plot(phase, sector_data.corr_flux_norm[q], '.k', alpha=0.8, label='Folded')
ax2.legend()
ax2.set(xlabel='Phase [-]',
          ylabel='flux')
# ax2.text(0.02, 0.03, "Period = {0:.2f} days".format(best_period),
#            transform=ax2.transAxes, bbox=dict(boxstyle="square"))
# ymin_auto, ymax_auto = ax2.yaxis.get_data_interval()
# ax2.set_ylim(ymin=np.max([ymin_auto, 0.0]))
# ax2.set_ylim(ymax=np.min([ymax_auto, 2.5]))

phase_model = np.linspace(0.0, 1.0, 100)

model = ls.model(phase_model * best_period, 1/best_period)

ax2.plot(phase_model, model, '-g', lw=5, alpha=0.8, zorder=2)
plt.savefig('eleanor_v773tau_periodogram.png', dpi=200, bbox_inches='tight')

# from scipy.optimize import curve_fit

# def sine(t, a, prd, phs, off):
#     return a*np.sin(2*np.pi*t/prd + phs) + off

# p0, pcov = curve_fit(sine, phase_model, model, 
#                      p0=[0.05, best_period, 0,0.], bounds=(0, [0.1, 2*best_period, np.pi, 1.]))
# ax2.plot(phase_model, sine(phase_model, *p0), '-r', lw=5, alpha=0.8)
# print(p0)
#%% Generate LombScargle periodogram for TESS data
ls = LombScargle(sector_data.time, sector_data.corr_flux)


def find_best_periods(periods, nperiods):
    best_periods = np.array([])
    add = True
    print('Picking the best periods...')
    for period in periods:
        add = True
        # print(period)

        if len(best_periods) >= nperiods:
            print('{:} best periods selected!!'.format(nperiods))
            return best_periods
        if np.min(period - np.array([1., 0.5, 0.33])) < 0.05:
            # print('1-day derived signal')
            add = False
        if len(best_periods) > 1:
            if np.min(period - best_periods) < 0.20:
                # print('Period already in the list')
                add = False
        if add == True:
            print('Adding {:.3f} d'.format(period))
            best_periods = np.append(best_periods, period)
    return best_periods

def generate_periodogram(ls, nbest_periods=1, ax=None, show=True):
    freq, powers = ls.autopower(nyquist_factor=50, minimum_frequency=1/10.,maximum_frequency=1/1.2)
    
    
    periods = np.array(1/freq)
    sort_peaks = np.flip(np.argsort(powers))
    periods  = periods[sort_peaks]
    best_periods = find_best_periods(periods, nperiods=nbest_periods)
    
    if show == True:
        ax = ax or plt.gca()
        ax.plot(1/freq, powers, '-k', label='{:.2f} d'.format(best_periods[0]))
        # ax.legend()
        ax.set(xlabel='Periods [days]', ylabel='LSC power')
        ax.grid()
        for i,best_period in enumerate(best_periods):
            ax.axvspan(best_period-0.05, best_period+0.05, alpha=0.6/(i+1), color='green')
    return ax, best_periods

fig, ax = plt.subplots(1, figsize=(9,5))
_, best_periods = generate_periodogram(ls, ax=ax, nbest_periods=2)
plt.show()        
        

#%% Fit sines to each best_period from the periodogram

fig, ax = plt.subplots(1, len(best_periods), figsize=(16,5), sharey=True)
for i, best_period in enumerate(best_periods):
    phase = (sector_data.time / best_period) % 1
    ax[i].plot(phase, sector_data.corr_flux, '.k', alpha=0.8, label='{:.2f} d'.format(best_period))
    ax[i].legend()
    ax[i].set(xlabel='Phase [-]', ylabel='flux')
    # ax[i].text(0.02, 0.03, "Period = {0:.2f} days".format(best_period),
    #            transform=ax[i].transAxes, bbox=dict(boxstyle="square"))
    # ymin_auto, ymax_auto = ax[i].yaxis.get_data_interval()
    # ax[i].set_ylim(ymin=np.max([ymin_auto, 0.0]))
    # ax[i].set_ylim(ymax=np.min([ymax_auto, 2.5]))
    
    phase_model = np.linspace(0.0, 1.0, 100)
    
    model = ls.model(phase_model * best_period, 1/best_period)
    # theta = ls.model_parameters(1/best_period)
    # print(theta.round(2))
    
    
    # mymodel = offset + theta[0] + ((theta[1]*np.sin(2*np.pi*phase_model / best_period)) + \
    #     (theta[2]*np.cos(2*np.pi*phase_model / best_period)))
    design_matrix = ls.design_matrix(1/best_period, phase_model)
    offset = ls.offset()
    mymodel = offset
    for theta in design_matrix:
        mymodel += theta[0] + ((theta[1]*np.sin(2*np.pi*phase_model / best_period)) + \
        (theta[2]*np.cos(2*np.pi*phase_model / best_period)))
    
    
    # design_matrix = ls.design_matrix(1/best_period, phase_model)
    # np.allclose(model, offset + design_matrix.dot(theta))
    
    ax[i].plot(phase_model, model, '-g', lw=5, alpha=0.8, zorder=2)
    ax[i].plot(phase_model, mymodel, '-r', lw=5, alpha=0.8, zorder=2)
    
#%%% Dirk's notebook
from scipy.optimize import curve_fit

q = sector_data.quality == 0
sector_data.corr_flux[q]/np.nanmedian(sector_data.corr_flux[q])
time, flux, flux_err = sector_data.time[q], sector_data.corr_flux[q]/np.nanmedian(sector_data.corr_flux[q]), \
    sector_data.corr_flux[q]*0.001/np.nanmedian(sector_data.corr_flux[q])


linear_trend = lambda P, time: supp.line(time, *P, time[0])
line_curvefit = lambda time, m, c: supp.line(time, m, c, time[0])
p_line, _ = curve_fit(line_curvefit, time, flux, sigma=flux_err)
line_flux = supp.line(time, *p_line, time[0])

# sanity check with a plot
_, = plot.plot_models(time, flux, flux_err, [linear_trend], [p_line], ['linear'])


print('the best fit parameters are m = %.12f, b = %.12f' % (p_line[0], p_line[1]))

#%%
# new flux (removed linear trend) for lomb_scargle
ls1_flux = flux - line_flux

def lsc_cycle(time, ls_flux, flux_err, ncycles=4, show=True):
    ls_p = []
    if show == True:
        fig, ax = plt.subplots(ncycles, figsize=(9,16))
        fig.tight_layout()
    for i in range(ncycles):
        if show == True:
            ls_p1b,ls2_flux = supp.lombscargle_periodogram_dirk(time, ls_flux, flux_err, ax=ax[i])
        else:
            ls_p1b,ls2_flux = supp.lombscargle_periodogram_dirk(time, ls_flux, flux_err, ax=None) # avoid plotting
        ls_p.append(ls_p1b)
        ls_flux = ls2_flux
    
    # plt.show()
    return np.array(ls_p)


ls_p = lsc_cycle(time, ls1_flux, flux_err, ncycles=6, show=False)

ls_df = pd.DataFrame(ls_p, columns=['amplitude', 'period', 'phase'])
# ls_df.insert(0, 'harmonic', False)

def del_harmonics(ls_df_in):
    ls_df = ls_df_in.copy()
    for i, fund_period in enumerate(ls_df.period):
        test_periods = ls_df.period.iloc[i+1:].values
        ratios = test_periods / fund_period 
        harmonics = np.round(ratios).astype(np.int) 
        discrepancies = 100 * np.abs(harmonics - ratios) 
        for test_period, harmonic, discrepancy in zip(test_periods, harmonics, discrepancies): 
                # print('the discrepancy between %.4f and %.4f is %.2f %% in the %i harmonic' % (fund_period, test_period, discrepancy, harmonic)) 
                if discrepancy < 5.:
                    out_period = min(test_period, fund_period) # remove the lower period (harmonic)
                    print('{:.2f} is a harmonic of {:.2f}...'.format(out_period, max(test_period, fund_period)))
                    ls_df.loc[(ls_df.period == out_period), 'harmonic'] = True
    len_before = len(ls_df)
    try:                
        ls_df = ls_df[ls_df.harmonic != True]
        nharmonics = len_before - len(ls_df)
        print('Removing {:} harmonic signals...'.format(nharmonics))
        del ls_df['harmonic']
    except AttributeError:
        pass
    return ls_df
    
ls_df = del_harmonics(ls_df)    


from scipy.stats import chisquare
# determine best fit parameters from lomb scargle
pb_sines = np.array(ls_p).T.flatten()
Pb_ls = np.concatenate((p_line, pb_sines), 0)

# define the model
print('Building stellar variation model...')
stellar_variation = supp.line(time, *p_line) + supp.sines(time, ls_df.amplitude, ls_df.period, ls_df.phase)


fig, ax = plt.subplots(figsize=(12,7))

print('Plotting model against TESS data...')
ax.plot(time, flux,'.b', alpha=0.7, label='TESS data')

flux_sv = stellar_variation / np.nanmedian(stellar_variation)
# chisq = chisquare(flux, flux_sv)
# ax.text(x=0.1, y =0.85)
ax.plot(time, flux_sv,'-g', lw=4, alpha=0.4, label='SV model {:} sines'.format(len(ls_df)))

ax.set(xlabel='Time', ylabel='Flux')
ax.legend()
plt.show()


#%% Determine the "reasonable" number of sines to fit the data
##### Not sure how to interpret the results (16th Sept.)
# For now use 10-sine terms

chisq = []
for k in range(2,30):
    ls_p = lsc_cycle(time, ls1_flux, flux_err, ncycles=k, show=False)
    ls_df = pd.DataFrame(ls_p, columns=['amplitude', 'period', 'phase'])
    # ls_df = del_harmonics(ls_df)   
    stellar_variation = supp.line(time, *p_line) + supp.sines(time, ls_df.amplitude, ls_df.period, ls_df.phase)
    flux_sv = stellar_variation / np.nanmedian(stellar_variation)
    chisq.append(chisquare(flux, flux_sv)[0])

# _, = plot.plot_models(time, flux, flux_err, [stellar_variation], [Pb_ls], ['stellar variation'])
fig, ax = plt.subplots(figsize=(12,7))

print('Plotting model against TESS data...')
ax.plot(time, flux,'.b', alpha=0.7, label='TESS data')

flux_sv = stellar_variation / np.nanmedian(stellar_variation)
# chisq = chisquare(flux, flux_sv)
# ax.text(x=0.1, y =0.85)
ax.plot(time, flux_sv,'-g', lw=4, alpha=0.4, label='SV model {:} sines'.format(len(ls_df)))

ax.set(xlabel='Time', ylabel='Flux')
ax.legend()
plt.show()
#%%
fig, ax = plt.subplots(1, figsize=(9,6))
ax.plot(range(len(chisq)), chisq, '--', label='{:} terms'.format(len(chisq)))
ax.set(xlabel='Number of sine terms', ylabel='Chi-square')

gradient = abs(np.diff(chisq) / np.array(chisq[1:])) 
best_nterms = np.argwhere(gradient < 0.05)[0][0] # arbitrary 5% change 
best_chisq = chisq[best_nterms-1] # one before its gradient becomes less than 5%
print('Plotting the evolution of Chi-Square with the number of sine-terms...')
print('The optimal number of sine-terms is {:}'.format(best_nterms))
ax.axvline(x=best_nterms, color='green', label='Optimal N-terms = {:}'.format(best_nterms), alpha=.5)
ax.legend()

###
ls_p = lsc_cycle(time, ls1_flux, flux_err, ncycles=best_nterms, show=False)
ls_df = pd.DataFrame(ls_p, columns=['amplitude', 'period', 'phase'])

# determine best fit parameters from lomb scargle
pb_sines = np.array(ls_p).T.flatten()
Pb_ls = np.concatenate((p_line, pb_sines), 0)

# define the model
print('Building stellar variation model...')
stellar_variation = supp.line(time, *p_line) + supp.sines(time, ls_df.amplitude, ls_df.period, ls_df.phase)

########## SAVE MODEL PARAMATERS TO CSV
print('Saving model parameters in .csv...')
df_model = pd.DataFrame([p_line, np.transpose(ls_df.values)])
df_model.to_csv('sv_model.csv')
##########

# _, = plot.plot_models(time, flux, flux_err, [stellar_variation], [Pb_ls], ['stellar variation'])
fig, ax = plt.subplots(figsize=(12,7))

print('Plotting model against TESS data...')
ax.plot(time, flux,'.b', alpha=0.7, label='TESS data')

flux_sv = stellar_variation / np.nanmedian(stellar_variation)
# ax.text(x=0.1, y =0.85)
ax.plot(time, flux_sv,'-g', lw=4, alpha=0.4, label='SV model {:} sines'.format(len(ls_df)))

ax.set(xlabel='Time', ylabel='Flux')
ax.legend()
plt.show()



