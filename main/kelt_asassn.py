"""
Created on Wed Jul 14 17:45:55 2021
@author: dario

Title: Read both KELT and ASAS_SN and merge data
"""
import numpy as np
import matplotlib.pyplot as plt
import support_functions as supp
import plot_functions as plot
import pandas as pd
import glob
import lightkurve as lk
from lmfit import Minimizer, Parameters, report_fit

# Get KELT files 
kelt_path = '/home/dario/AstronomyLeiden/FRP/MRP-ASAS-Rings/KELT_PMS/lc/'
kelt_files = glob.glob(kelt_path+'*.lc')
kelt_ids = []
for file in kelt_files:
    kelt_ids.append(int(file.split('/')[-1][:-3]))
print('{:} KELT files found'.format(len(kelt_ids)))

# Get ASAS_SN gaia_ids
basepath_V = '/home/dario/AstronomyLeiden/FRP/leiden_vband/camfix/'
asas_sn = glob.glob(basepath_V+'*.dat')
asas_sn_lcs = supp.LightCurveSet(basepath_V, 'dat', 'v')
print('{:} ASAS_SN files found'.format(len(asas_sn_lcs.files(gaia_id=True))))

# Find matches between both catalogs
intersect = list(set(asas_sn_lcs.files(gaia_id=True)).intersection(kelt_ids))
print('-----> {:} targets matched'.format(len(intersect)))

df_kelt = supp.load_kelt(intersect[12])
df_asassn = supp.LightCurveObject(intersect[12],'v').data
#
#%%
gaia_id = intersect[12]
outname = '/home/dario/AstronomyLeiden/FRP/MRP-ASAS-Rings/main/merged_lcs/' + str(gaia_id) +'.csv'

fig, ax = plt.subplots(1, figsize=(9,9))
df, ax = supp.kelt_asassn(intersect[12], ax=ax, show=True)
# df.drop(columns='FWHM', inplace=True)
# df.to_csv(outname, index=False)
df = df.round({'flux':5, 'flux_err':5})
plt.show()

#%%
# Load KELT "interesting" targets
files = glob.glob('/home/dario/AstronomyLeiden/FRP/MRP-ASAS-Rings/kelt-interesting/*.png')

kelt_gaia_ids_int = []
# get gaia_ids
for file in files:
    temp = file.split('/')[-1][5:-4]
    kelt_gaia_ids_int.append(int(temp))
    



#%%
plt.ioff() # disable plt plotting automatically to save CPU
outlist = []
'''
for gaia_id in kelt_gaia_ids_int[12:13]:
    try:
        df, _ = supp.kelt_asassn(gaia_id, ax=None, write=False)    
        mean_mag = np.nanmedian(df['mag'])
        mean_mag_str = 'Mean V-mag = {:.1f}'.format(mean_mag)
        
        df_kelt = df[df['survey']=='KELT'] # warning: there is an unwanted /s before "KELT"
        df_kelt_rs = supp.resample(df_kelt, t_bin=1.)
        
        df = pd.concat([df[df['survey']=='ASAS_SN'], df_kelt_rs])
        # Make a 3-fig plot for classification on zooniverse   
        fig = plt.figure(figsize=(14,9))
        gs = fig.add_gridspec(2, 2)
        
        
        ax1 = fig.add_subplot(gs[0,:])
        # ax1.invert_yaxis()
        label='V band'
        # ax1.errorbar(time, clean_flux+1, errflux, fmt='.k', alpha=0.8, label=label)
        time_corr = df['HJD'] - 2450000 # transform HDJ date
        # ax1.errorbar(time_corr, flux, errflux, fmt='.k', alpha=0.6, label=label)
        # ax1.plot(time_corr, df['flux'],'.k', alpha=0.6, label=label)
        
        ax1.set_title('GAIA_ID = ' + str(gaia_id))
        # ax1.set_ylabel('Relative flux')
        ax1.set_ylabel('Flux')
        ax1.set_xlabel('HJD - 2450000')
        
    
        _ = plot.plot_survey(df, ax=ax1, errorbars=True)
        ax1.legend()
        
        ymin_auto, ymax_auto = ax1.yaxis.get_data_interval()
        ax1.set_ylim(ymin=np.max([ymin_auto, 0.0]))
        ax1.set_ylim(ymax=np.min([ymax_auto, 1.5]))
        ax1.grid()
        
        
        ax2 = fig.add_subplot(gs[1,0])
        
        time, flux, flux_err = df['HJD'], df['flux'], df['flux_err']
        # WARNING: flux must be around 0.0 for the periodogram to work
        lc = lk.LightCurve(time, flux, flux_err)
        
        periodogram = lk.periodogram.LombScarglePeriodogram.from_lightcurve(lc, 
                                                                            maximum_period=10,
                                                                            nterms=1)
        # periodogram.show_properties()
        # period = periodogram.period_at_max_power.value
        freq = periodogram.frequency_at_max_power.value
        periods = np.array(periodogram.period)
        sort_peaks = np.flip(np.argsort(periodogram.power))
        periods  = periods[sort_peaks]
        peaks = np.array(periodogram.power[sort_peaks])
        
                
        # select peak and remove 1-day period signals
        peak_ind = 0
        try:
            period = periods[peak_ind]
        except:
            print('No peaks found on Periodogram')
            period = 1.0 # day
        while (abs(period - 1.0) < 0.1) or (abs(period - 0.50) < 0.1) or (abs(period - 0.33) < 0.1)  :
            print('Removed 1-day derived period signal')
            peak_ind += 1 
            period = periods[peak_ind]
                
        lc_folded = lc.fold(period=period)
        
        
        ax3 = fig.add_subplot(gs[1,1])
        periodogram.default_view = 'period'
        _ = periodogram.plot(ax=ax3, color='black')
        # ax3.set_xlim(0,10)
        ax3.set_ylim(bottom=0)
        # xticks = np.arange(0,np.max(freq),1)
        # ax3.set_xticks(xticks)
        
        # create data to be fitted
        x = lc_folded.time
        data = lc_folded.flux
        
        
        # define objective function: returns the array to be minimized
        def fcn2min(params, x, data):
            """Model a sine wave and subtract data."""
            amp = params['amp']
            shift = params['shift']
            period = params['period']
            offset = params['offset']
            model = offset + (amp * np.sin((2*np.pi*x/period) + shift))
            return model - data
        
        
        # create a set of Parameters
        params = Parameters()
        params.add('amp', value=0.1, min=0, max=0.8)
        params.add('shift', value=0.0, min=-np.pi/2., max=np.pi/2.)
        params.add('period', value=period)
        params.add('offset', value=1.)
        
        # do fit, here with the default leastsq algorithm
        minner = Minimizer(fcn2min, params, fcn_args=(x, data))
        result = minner.minimize()
        
        # calculate final result
        final = data + result.residual
        
        # write error report
        report_fit(result)
        
        ax2.plot(x, final, 'r', lw=4, alpha=0.7)
        
        
        
        # plot the folded light curve
        ax2.errorbar(lc_folded.time, lc_folded.flux,fmt='.',
                     label="folded -- {:.2f} d".format(period), 
                     alpha=0.7, ms=3, zorder=1, color='black')
        ax2.legend()
        ax2.set_ylabel("Flux")
        ax2.set_xlabel("Phase")
        # ax2.set_ylim(0.7, 1.2)
        ymin_auto, ymax_auto = ax2.yaxis.get_data_interval()
        y_min = np.max([ymin_auto, 0.5])
        y_max = np.min([ymax_auto, 1.4])
        ax2.set_ylim(y_min, y_max)
        outpath = '/home/dario/AstronomyLeiden/FRP/MRP-ASAS-Rings/main/merged/kelt_interesting/' + str(gaia_id) +'.png'
        outlist.append(outpath)
        plt.savefig(outpath, dpi=200, bbox_inches='tight')
        print('Saved {:}'.format(outpath))
        # fig.clear()
        # plt.close(fig)
        plt.show()
    except:
        print('ParserError')
        pass
'''

#%% Finally, upload the LCs to zooniverse
import panoptes as pan
# _ = pan.panoptes_add(outlist, 'kelt_asassn_interesting2')


#%%
'''
# plot individual target
gaia_id = 6039427503765943936

fig, ax = plt.subplots(figsize=(14,6))
df, _ = supp.kelt_asassn(gaia_id, ax=None, show=False)

_ = plot.plot_survey(df, ax=ax, errorbars=False)
ax.legend()
ax.grid()
ax.set_ylim(top=1.5)

plt.show()
outpath = '/home/dario/AstronomyLeiden/FRP/MRP-ASAS-Rings/' + str(gaia_id) +'.png'
plt.savefig(outpath, dpi=200, bbox_inches='tight')
'''

#%% Updated version from Sept. 13th 2021
# KELT data is rebinned and displayed along ASAS_SN in ax0
# The LombScargle is generated with the NON-REBINNED data
# The "manual" folded light curve is plotted with the LS.MODEL() fit
'''
plt.ioff() # disable plt plotting automatically to save CPU


def three_plot_df(df, gaia_id):
    # Rebin ONLY KELT data
    df_kelt = df[df['survey']=='KELT'] 
    df_kelt_rs = supp.resample(df_kelt, t_bin=1.)
    
    # Merge "resampled" light curve
    df_rs = pd.concat([df[df['survey']=='ASAS_SN'], df_kelt_rs])
    df_rs = supp.pd_sigma(df_rs, y='flux_err', n=2)
    
    # Make a 3-fig plot for classification on zooniverse   
    fig = plt.figure(figsize=(14,9))
    gs = fig.add_gridspec(2, 2)
    
    
    ax0 = fig.add_subplot(gs[0,:])
    
    _ = plot.plot_lc(df_rs, ax=ax0, y='flux', color='navy')
    ax0.grid()
    mean_mag = np.nanmedian(df['mag'])
    mean_mag_str = 'Mean V-mag = {:.1f}'.format(mean_mag)
    ax0.text(0.8, 1.03, s=mean_mag_str, transform=ax0.transAxes, alpha=0.7)
    ax0.set_title('GAIA_ID = {:}'.format(gaia_id))
    
    df = supp.pd_sigma(df, y='flux', n=2)
    ls = LombScargle(df.HJD, df.flux, df.flux_err)
    frequency, power = ls.autopower(nyquist_factor=50,
                                    minimum_frequency=1/100.,
                                    maximum_frequency=1/1.2)
    
    period_days = 1. / frequency
    period_hours = period_days * 24
    

    best_period = supp.pick_period(frequency, power)
    phase = (df.HJD / best_period) % 1
    
    print("Best period: {:.2f} days".format(best_period))

    ax1 = fig.add_subplot(gs[1,0])
    ax1.plot(period_days, power, '-k')
    ax1.axvline(x=best_period, ls='-', lw=8, alpha=0.3, label='{:.2f} d'.format(best_period), color='green')
    # ax1.fill_betweenx(x1=best_period-0.2, x2=best_period+0.2, alpha=0.4)
    ax1.legend()
    ax1.grid()
    
    ax1.set(xlabel='Period (days)', ylabel='Lomb-Scargle Power')
    ymin_auto, ymax_auto = ax1.yaxis.get_data_interval()
    ax1.set_ylim(ymax=np.max([ymax_auto, 0.10]))
    
    
    ax2 = fig.add_subplot(gs[1,1])
    # ax2.errorbar(phase, df.flux, df.flux_err,
    #                fmt='.k', ecolor='gray', ms=1, capsize=0, zorder=1)
    ax2.plot(phase, df.flux, '.k', alpha=0.8, label='Folded')
    ax2.legend()
    ax2.set(xlabel='Phase [-]',
              ylabel='flux')
    # ax2.text(0.02, 0.03, "Period = {0:.2f} days".format(best_period),
    #            transform=ax2.transAxes, bbox=dict(boxstyle="square"))
    ymin_auto, ymax_auto = ax2.yaxis.get_data_interval()
    ax2.set_ylim(ymin=np.max([ymin_auto, 0.0]))
    ax2.set_ylim(ymax=np.min([ymax_auto, 2.5]))
    
    phase_model = np.linspace(0.0, 1.0, 100)

    mag_model = ls.model(phase_model * best_period, 1/best_period)
    
    ax2.plot(phase_model, mag_model, '-g', lw=5, alpha=0.8, zorder=2)
    outpath = '/home/dario/AstronomyLeiden/FRP/MRP-ASAS-Rings/main/merged/zooniverse/' + str(gaia_id) +'.png'
    
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    print('Saved {:}'.format(outpath))
    fig.clear()
    plt.close(fig)
    plt.show()
    return outpath


outlist = [] # to upload to zooniverse later
for gaia_id in intersect:
    try:
        df, _ = supp.kelt_asassn(gaia_id, ax=None, show=False)
        outpath = three_plot_df(df, gaia_id)
        outlist.append(outpath)
    except:
        print('Error')
'''

#%% List of Interesting targets by Dirk (Sept 14th 2021)

gaia_id = 146719488741654784
gaia_id = 149540560701692672

gaia_id_zoomin = [146719488741654784, 149540560701692672, 153049407246313728,
                  164832740220756608, 171858275923130112, 2246820366041389056,
                  3420162069318601728, 3429395596170205056, 3442736005110139264,
                  6245762538722565888]

gaia_id_dips = [[151262700852297728, 5600], [163184366130809984, 5400], 
               [3017010295452313984, 6650], [3210151298261735808, 5550],
               [3216729878197029120, 5500], [3221372600399707520, 8250],
               [3409711868428222848, 6300], [5655070120227528448, 6600],
               [5919033999875261568, 7000], [5905886387735366528, 8100],
               [5955809537093034112, 8300], [6059825467381589888, 8300],
               [6112637935742759680, 8300], [6112762829087987712, 8300]]

other = [1563865079295954560, 3415797729710946688, 4167204942506094208,
         4379892204038009856, 4476241029294484864, 4488729561538989696,
         6072893029952309632]

####

def zoom_in(gaia_id):
    df, _ = supp.kelt_asassn(gaia_id, ax=None, write=False)    
    mean_mag = np.nanmedian(df['mag'])
    mean_mag_str = 'Mean V-mag = {:.1f}'.format(mean_mag)
    
    df_kelt = df[df['survey']=='KELT'] # warning: there is an unwanted /s before "KELT"
    # df_kelt = supp.pd_sigma(df_kelt, y='flux_err', n=1)
    df_kelt_rs = supp.resample(df_kelt, t_bin=1.)
    df_kelt_rs = supp.pd_sigma(df_kelt_rs, y='flux_err', n=2)
    df_kelt_rs = supp.pd_sigma(df_kelt_rs, y='flux', n=2)
    
    df_merged = pd.concat([df_kelt_rs, df[df['survey']=='ASAS_SN']])
    df_merged = supp.pd_sigma(df_merged, y='flux', n=3)
    # Make a 3-fig plot for classification on zooniverse   
    
    #### Get y-limits without the errorbars (this plot is not displayed)
    # fig, ax0 = plt.subplots(1)
    # ax0.plot(df.HJD, df.flux,'.')
    # ymin_auto, ymax_auto = ax0.yaxis.get_data_interval()
    # fig.clear()
    # plt.close(fig)
    ##################
    
    fig, ax = plt.subplots(2, figsize=(12,9), sharex=True)
    
    ax[0].set_title('GAIA_ID = ' + str(gaia_id))
    # ax[0].set_ylabel('Relative flux')
    ax[0].set_ylabel('Flux')
    # ax[0].set_xlabel('HJD - 2450000')
    ax[0].grid()
    ax[0].text(x=0.05, y=0.90, s='Zoomed in', transform=ax[0].transAxes, alpha=0.7)
    
    # ax[0].set_ylim(ymin_auto, ymax_auto) # apply y-lims from hidden plot
    
    
    _ = plot.plot_survey(df_merged, ax=ax[0], errorbars=True)
    ax[0].legend()
    
    _ = plot.plot_survey(df, ax=ax[1], errorbars=True)
    ax[1].text(x=0.05, y=0.90, s='Original view', transform=ax[1].transAxes, alpha=0.7)
    ax[1].legend()
    ax[1].grid()
    plt.show()
    outpath = '/home/dario/AstronomyLeiden/FRP/MRP-ASAS-Rings/main/merged/dirk_list/' + str(gaia_id) +'.png'
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    print('Saved {:}'.format(outpath))
    return None


# for gaia_id in gaia_id_zoomin:
#     _ = zoom_in(gaia_id)
import matplotlib.patches as mpatches

def dip(gaia_id_dip):
    gaia_id, dip = gaia_id_dip 
    print('Dip at {:}'.format(dip))
    df, _ = supp.kelt_asassn(gaia_id, ax=None, write=False)    
    mean_mag = np.nanmedian(df['mag'])
    mean_mag_str = 'Mean V-mag = {:.1f}'.format(mean_mag)
    
    df_kelt = df[df['survey']=='KELT'] # warning: there is an unwanted /s before "KELT"
    # df_kelt = supp.pd_sigma(df_kelt, y='flux_err', n=1)
    df_kelt_rs = supp.resample(df_kelt, t_bin=1.)
    # df_kelt_rs = supp.pd_sigma(df_kelt_rs, y='flux_err', n=2)
    # df_kelt_rs = supp.pd_sigma(df_kelt_rs, y='flux', n=2)
    
    df_merged = pd.concat([df_kelt_rs, df[df['survey']=='ASAS_SN']])
    # df_merged = supp.pd_sigma(df_merged, y='flux', n=3)

    
    fig, ax = plt.subplots(2, figsize=(12,9))
    
    
    _ = plot.plot_survey(df_merged, ax=ax[0], errorbars=True)
    ax[0].set_title('GAIA_ID = ' + str(gaia_id))
    ax[0].text(0.8, 1.03, s=mean_mag_str, transform=ax[0].transAxes, alpha=0.7)
    # ax[0].set_ylabel('Relative flux')
    ax[0].set_ylabel('Flux')
    ax[0].set_xlabel(None)
    sides = 150
    ax[0].set_xlim(dip-sides, dip+sides)
    ax[0].grid()
    dip_patch = mpatches.Patch(color='wheat', label='Dip at {:}'.format(dip))
    handles, labels = ax[0].get_legend_handles_labels()
    handles.append(dip_patch) 
    ax[0].legend(handles=handles)
    # props = dict(boxstyle='square', facecolor='none', alpha=0.2)
    # ax[0].text(x=0.05, y=0.90, s=, 
    #            transform=ax[0].transAxes, alpha=0.7, bbox=props)
    
    # ax[0].set_ylim(ymin_auto, ymax_auto) # apply y-lims from hidden plot
    
    
    
    
    
    
    _ = plot.plot_survey(df_merged, ax=ax[1], errorbars=True)
    ax[1].axvspan(dip-sides, dip+sides, alpha=0.4, color='wheat')
    ax[1].text(x=0.05, y=0.90, s='Original view', transform=ax[1].transAxes, alpha=0.7)
    # ax[1].legend()
    ax[1].grid()
    plt.show()
    outpath = '/home/dario/AstronomyLeiden/FRP/MRP-ASAS-Rings/main/merged/dirk_list/' + str(gaia_id) +'.png'
    plt.savefig(outpath, dpi=200, bbox_inches='tight')
    print('Saved {:}'.format(outpath))
    return None


# for gaia_id_dip in gaia_id_dips:
#     _ = dip(gaia_id_dip)

gaia_id = other[7]
# _ = dip([gaia_id, 6800])
df, _ = supp.kelt_asassn(gaia_id, ax=None, write=False)    
mean_mag = np.nanmedian(df['mag'])
mean_mag_str = 'Mean V-mag = {:.1f}'.format(mean_mag)

df_kelt = df[df['survey']=='KELT'] # warning: there is an unwanted /s before "KELT"
# df_kelt = supp.pd_sigma(df_kelt, y='flux_err', n=1)
df_kelt_rs = supp.resample(df_kelt, t_bin=1.)
df_kelt_rs = supp.pd_sigma(df_kelt_rs, y='flux_err', n=2)
df_kelt_rs = supp.pd_sigma(df_kelt_rs, y='flux', n=2)

df_merged = pd.concat([df_kelt_rs, df[df['survey']=='ASAS_SN']])
# df_merged = supp.pd_sigma(df_merged, y='flux', n=3)

fig, ax = plt.subplots(2, figsize=(12,9), sharex=True)

ax[0].set_title('GAIA_ID = ' + str(gaia_id))
# ax[0].set_ylabel('Relative flux')
ax[0].set_ylabel('Flux')
# ax[0].set_xlabel('HJD - 2450000')
ax[0].grid()
ax[0].text(x=0.05, y=0.90, s='Zoomed in', transform=ax[0].transAxes, alpha=0.7)

# ax[0].set_ylim(ymin_auto, ymax_auto) # apply y-lims from hidden plot


_ = plot.plot_survey(df_merged, ax=ax[0], errorbars=True)
ax[0].legend()

_ = plot.plot_survey(df, ax=ax[1], errorbars=True)
ax[1].text(x=0.05, y=0.90, s='Original view', transform=ax[1].transAxes, alpha=0.7)
ax[1].legend()
ax[1].grid()
plt.show()
outpath = '/home/dario/AstronomyLeiden/FRP/MRP-ASAS-Rings/main/merged/dirk_list/' + str(gaia_id) +'.png'
# plt.savefig(outpath, dpi=200, bbox_inches='tight')
print('Saved {:}'.format(outpath))
    












