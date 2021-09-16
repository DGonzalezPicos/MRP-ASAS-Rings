"""
Created on Tue Jun  8 12:46:27 2021
@author: dario

Title: Support functions for main script
"""
#%% STANDARD MODULES %%#
########################

from astropy.timeseries import LombScargle
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import numpy as np
from astropy.stats import sigma_clip
from astropy.io import ascii 
from scipy.stats import chisquare
import pandas as pd
import warnings
import glob, os
from scipy.optimize import curve_fit
import plot_functions as plot
from scipy.stats import chisquare

##########################
#%% ANALYSIS FUNCTIONS %%#
##########################

class LightCurveSet:
    def __init__(self, path, extension, band):
        self.path = path
        self.extension = '.' + extension
        self.band = band
    
    
    def files(self, gaia_id=False):
        '''
        Return the list of files at a given folder.
        
        Parameters
        ----------
        gaia_id : Boolean (optional)
            If True return only the GAIA_ID for each file. The default is False.

        Returns
        -------
        TYPE
            List of files / names.

        '''
        
        files = glob.glob(os.path.join(self.path, '*' + self.extension))
        
        # Return a list with ONLY the GAIA_ID for each LC
        if gaia_id == True:
            ids = []
            for row in files:
                ids.append(int(row.split('/')[-1][:-4]))
            return ids
        return files

    
    def data(self, file):
        '''
        Return LC data Objects for either band

        Parameters
        ----------
        file : string
            Full file-path.

        Returns
        -------
        df : pd.DataFrame
            LC dataframe with 4 columns [HJD, flux, errflux, gaia_id]
            flux, errflux are relative

        '''
        assert self.extension == '.dat', 'Extension must be .dat'
            
        if self.band == 'v':
            data = ascii.read(file, delimiter='\t')
            df = data.to_pandas()
            df['flux_err'] = df['flux_err'].replace(99.99, np.nan)
            df['mag_err'] = df['mag_err'].replace(99.99, np.nan)
            df = df.dropna()
            name = file.split('/')[-1][:-4]
            df['gaia_id'] = np.full_like(df['HJD'], name)
            ##1 convert to relative flux
            flux_med = np.nanmedian(df['flux'])
            df['flux'] = df['flux'] / flux_med
            df['flux_err'] = df['flux_err'] / flux_med
            ##1
            
        if self.band == 'g':
            data =  ascii.read(file, format='no_header')
            df = data.to_pandas()
            df.columns = ['HJD', 'flux', 'flux_err', 'mag','mag_err', 'limit', 'fwhm', 'quality', 'cam', 'asas_sn_id']
        # if self.band == 'g':
        #    # data is [hjd, flux, errflux, asas_id]
        #     df = pd.DataFrame(np.loadtxt(file,usecols=[0,1,2,-1]), columns=['HJD', 'flux','flux_err','asas_id'])
        #     df.drop(columns=['asas_id']) # remove ASAS_SN_ID 
        #     df = df.dropna()
        #     name = int(file.split('/')[-1][:-4])
        #     df['gaia_id'] = np.full_like(df['HJD'], name)
            
            # convert to relative flux
            flux_med = np.nanmedian(df['flux'])
            df['flux'] = df['flux'] / flux_med
            df['flux_err'] = df['flux_err'] / flux_med
            
        if self.band == 'merged':
           # data is [hjd, flux, flux_err, gaia_id, asas_id, mean_mag]
            df = pd.DataFrame(np.loadtxt(file,usecols=[0,1,2,3,-1]), columns=['HJD', 'flux','flux_err','gaia_id', 'mean_mag'])
            # convert to relative flux
            flux_med = np.nanmedian(df['flux'])
            df['flux'] = df['flux'] / flux_med
            df['flux_err'] = df['flux_err'] / flux_med
            # deal with 'nan'---> remove them
            # df = df.dropna()
            # df['errflux'] = df['errflux'].replace(np.nan, np.nanmedian(df['errflux']))
        
        return df
    
class LightCurveObject:
    def __init__(self, gaia_id, band):
        self.id = gaia_id
        self.filter = band
    
    @property
    def path(self):
        vpath = '/home/dario/AstronomyLeiden/FRP/leiden_vband/camfix/'
        return os.path.join(vpath, str(self.id))
    @property
    def data(self):
        df = get_data(self.id)
        df['Filter'] = np.full_like(df['HJD'], self.filter, dtype=str)
        return get_data(self.id)
    
    def plot(self, ax=None):
        df = self.data
        df['mag'] = df['mag'].astype('float')
        
        ax = ax or plt.gca() # function embedded
        # fig, ax = plt.subplots(1, figsize=(16,7))
        
        cams = df['camera'].unique()
        colors = ['brown','darkorange','chocolate']
        print('Plotting {:} with cameras {:}'.format(self.id, cams))
        for i, cam in enumerate(cams):
            _ = plot.plot_cam(df, cam=cam, y='mag', ax=ax, label=cam, color=colors[i])
            
        ax.invert_yaxis()   
        ax.set_title('GAIA_ID = ' + str(self.id))
        ax.set_ylabel('Magnitude')
        ax.set_xlabel('HJD - 2450000')
        
        ax.legend()
        plt.show()
        return ax
        
def get_data(gaia_id):
    vpath = '/home/dario/AstronomyLeiden/FRP/leiden_vband/camfix/'
    file = os.path.join(vpath, str(gaia_id)+'.dat')
    data = ascii.read(file, delimiter='\t')
    df = data.to_pandas()
    df['flux_err'] = df['flux_err'].replace(99.99, np.nan)
    df['mag_err'] = df['mag_err'].replace(99.99, np.nan)
    df = df.dropna()
    ##1 convert to relative flux
    flux_med = np.nanmedian(df['flux'])
    df['flux'] = df['flux'] / flux_med
    df['flux_err'] = df['flux_err'] / flux_med
    return df    


def read_gband(file):
    data = np.loadtxt(file)

    hjd, errflux = data[:,0], data[:,2]    
        
    flux = []
    for item in data[:,1]:
        a  = str(type(item))
        if a == "<class 'numpy.str_'>":
            flux.append(float(item.split('>')[-1]))
        else:
            flux.append(item)
 
    name = file.split('/')[-1][:-4]
    
    # Fix "infinite values" for the "errflux" data by
    # replacing them with the median
    med = np.median(errflux)
    if med > 90:
        print('Error bars not valid for plot ', name)
    
    errflux = np.array(errflux)
    med = np.median(errflux[errflux<90])    # compute the median clipping values >90
    errflux[errflux>90] = med               # replace conflicting values with the median
    
    # return the relative flux, errflux
    flux_rel = flux / np.median(flux)
    errflux_rel = errflux / np.median(flux)
    
    return hjd, flux_rel, errflux_rel, name
def read_vband(gaia_id):
    path = '/home/dario/AstronomyLeiden/FRP/leiden_vband/camfix/'
    file = path + str(gaia_id) +'.dat'
    data = ascii.read(file, delimiter='\t')
    df = data.to_pandas()
    df['flux_err'] = df['flux_err'].replace(99.99, np.nan)
    # df = df.replace(99.99, np.nan)
    df = df.dropna()
    # name = file.split('/')[-1][:-4]
    
    # convert to relative flux
    flux_med = np.nanmedian(df['flux'])
    df['flux'] = df['flux'] / flux_med
    df['flux_err'] = df['flux_err'] / flux_med
    print(df)
    
    return df

def load_kelt(gaia_id):
    kelt_path = '/home/dario/AstronomyLeiden/FRP/MRP-ASAS-Rings/KELT_PMS/lc/'
    file = kelt_path + str(gaia_id) +'.lc'
    df_kelt = pd.read_csv(file, sep=" ", header=None)
    df_kelt = df_kelt.dropna(axis=1, how='all')
    df_kelt = df_kelt.rename(columns={3:'HJD', 4:'mag',5:'mag_err',6:'camera'})
    df_kelt['HJD'] = df_kelt['HJD'] + 2450000
    
    # generate flux data
    # df_kelt = flux_from_mag(df_kelt)
    return df_kelt

def resample(df_in, t_bin=1.):
    '''
    Given a light curve DataFrame, resample the data into bins of size "t_bin" 
    It applies astropy's sigma_clip to the flux

    Parameters
    ----------
    df_in : DataFrame
        light curve data.
    t_bin : integer
        Bin size in days.
    eclipse: Boolean
        True to resample eclipse or False to skip in-eclipse points

    Returns
    -------
    df_rs : DataFrame
        Resampled light curve.
    '''
    print('------> Rebinning light curve to {:}-day bins'.format(int(t_bin)))
    print('Number of points RAW light curve...', len(df_in))
    df_in = df_in.sort_values('HJD') # important for next step
    
    t_rebin = np.arange(df_in['HJD'].iloc[0], df_in['HJD'].iloc[-1], t_bin)
    time, flux, flux_err, boxes = ([] for i in range(4))
    for t in t_rebin:
        start = t - t_bin/2.
        end = t + t_bin/2.
        boxes.append((start, end)) # only for displaying-purposes
        
        df = df_in[(df_in['HJD'] > start) & (df_in['HJD'] < end)]
        if len(df) > 0:
            flux_clipped = sigma_clip(df['flux'], sigma=3, maxiters=2)
            meanflux = np.mean(flux_clipped)
            
            time.append(t)
            flux.append(meanflux)
            flux_err.append(np.std(flux_clipped))
            
    df_rs = pd.DataFrame(np.transpose([time, flux, flux_err]), columns=('HJD', 'flux', 'flux_err'))
    df_rs.insert(0, 'survey', df_in['survey'].iloc[0])

    print('Number of points RESAMPLED light curve...', len(df_rs))
    print('------> Done!')
    return df_rs

def kelt_asassn(gaia_id, ax=None, save=False, show=False, write=False):
    '''
    Merge data from KELT and ASAS_SN

    Parameters
    ----------
    gaia_id : integer

    Returns
    -------
    df : pd.DataFrame
        merged data
    ax : plot
        merged plot.

    '''
    df_kelt = load_kelt(gaia_id)
    df_kelt.insert(0, 'survey','KELT')
    
    df_asas_sn = LightCurveObject(gaia_id,'v').data
    df_asas_sn.insert(0, 'survey','ASAS_SN')
    mean_mag = np.nanmedian(df_asas_sn['mag'])
    
    # Merge data from KELT AND ASAS_SN by normalising the flux for each camera
    
    df = pd.concat([df_kelt, df_asas_sn])
    for cam in df['camera'].unique():
        df_cam = pd.DataFrame(df[df['camera']==cam])
        df_cam = flux_from_mag(df_cam)
        df.loc[df['camera']==cam,'flux'] = df_cam['flux']
        df.loc[df['camera']==cam,'flux_err'] = df_cam['flux_err']
    
        
        # print(cam)
    
        
    if show == True:
        ax = ax or plt.gca()
        for cam in df['camera'].unique():
            _ = plot.plot_cam(df, cam, 'flux', ax=ax, label=cam)
        
        ax.legend()
        ax.set_xlabel('HJD -2450000')
        ax.set_ylabel('Normalised flux')
        ax.set_ylim(0.25,1.4)
        ax.text(s='V = {:.1f}'.format(mean_mag), x=0.85, y=1.01, 
                transform=ax.transAxes, fontsize=16)
        ax.set_title('GAIA_ID = {:}'.format(gaia_id))
        # plt.show()
    
    if save == True:
        outname = '/home/dario/AstronomyLeiden/FRP/MRP-ASAS-Rings/main/merged/' + str(gaia_id) +'.png'
        outlist.append(outname)
        plt.savefig(outname, dpi=200, bbox_inches='tight')
        print('Saved {:}'.format(outname))
        fig.clear()
        plt.close(fig)
        
    if write == True:
        outname = '/home/dario/AstronomyLeiden/FRP/MRP-ASAS-Rings/main/merged/merged_lcs/' + str(gaia_id) +'.csv'
        print('Writing {:}'.format(outname))
        df.drop(columns='FWHM', inplace=True)
        df = df.round({'flux':6, 'flux_err':6})
        df.to_csv(outname, index=False)
    

    return df, ax

def flux_from_mag(df, ax=None):
    df['flux'] = 10**(-df['mag']/2.5) /  10**(-np.nanmedian(df['mag'])/2.5)
    df['flux_err'] = df['flux'] * df['mag_err'] / 1.09
    return df

def load_asas(file, ax=None):
    
    hdr = ['HJD', 'MAG_0', 'MAG_1',  'MAG_2',  'MAG_3',  'MAG_4', 'MER_0', 'MER_1',
               'MER_2', 'MER_3', 'MER_4','GRADE','FRAME','FLAG']
    data = ascii.read(file, delimiter='\s', names=hdr)
    df = data.to_pandas()
    if df['HJD'][0] < 2450000: # only for the first three files
        df['HJD'] = df['HJD'] + 2450000
    df = df[df['GRADE']=='A']
    df['file_num'] = file.split('.')[-2]
    
    
    # fig, ax = plt.subplots(1, figsize=(12,5))
    ax = ax or plt.gca()
    for i in range(1):
        mag = 'MAG_'+str(i)
        errmag = 'MER_'+str(i)
        ax.errorbar(df['HJD'], df[mag], df[errmag], fmt='.', label=mag)
    return df

def pd_sigma(df, y='mag_err', n=1):
    '''
    Clip values "n"-sigma times beyond the mean on a given column "y"

    Parameters
    ----------
    df : DataFrame
        light curve data
    y : str, optional
        Column to be clipped. The default is 'mag_err'.
    n : int, optional
        Number of sigmas. The default is 1.

    Returns
    -------
    df : DataFrame
        Clipped dataframe.

    '''
    len_in = len(df)
    median = np.nanmedian(df[y])
    std = np.std(df[y])
    sigma = n*std
    df = df[abs(df[y] - median) < sigma]
    len_out = len(df)
    print('Sigma-clipping: {:} points discarded at {:}-sigma'.format(len_in - len_out, n))
    return df
#-----------------------------------------------------------------------------
def iterative_clipping(time, flux, errflux, eps=0.05, maxiter=10):
    '''
    Perform iterative sigma clipping until relative improvement on standard 
    deviation is less than eps (e.g. 5%)
    
    Parameters
    ----------
    flux : array
    errflux: array 
    eps : float
        threshold 

    Returns
    -------
    cleaned: time, flux, errflux

    '''
    time_clip, flux_clip, errflux_clip = [np.array([]) for i in range(3)]
    change = 1.0 
    count = 0
    while change > eps:
        sigmaclip = sigma_clip(errflux, sigma=3)
        mask = sigmaclip.mask
        
        time_clip = np.append(time_clip, time[mask])
        flux_clip = np.append(flux_clip, flux[mask])
        errflux_clip = np.append(errflux_clip, errflux[mask])
        
        change = 1 - np.std(flux[~mask]) / np.std(flux)
        print('{:2.8f}'.format(change))
        time, flux, errflux = time[~mask], flux[~mask], errflux[~mask]
        
        count += 1
        if count > maxiter:
            print('Max number of interations: ', count)
            break
    return time, flux, errflux, count


def sigma_clip_lc(flux, errflux, sigma_lower=2, sigma_upper=4, maxiters=3):
    '''
    Generate combined mask to remove outliers for 'flux value' and for the
    'errorbars'

    Parameters
    ----------
    flux : array
        Values of the LC.
    errflux : array
        Errorbars for the flux.
    sigma_lower : integer
        DESCRIPTION. The default is 2.
    sigma_upper : integer
        DESCRIPTION. The default is 4.
    maxiters : clipping iterations
        DESCRIPTION. The default is 3.

    Returns
    -------
    mask : boolean array
        Array with outliers as TRUE.
    '''
    sigma_flux = sigma_clip(flux, sigma_lower, sigma_upper, maxiters, masked=True)
    sigma_errflux = sigma_clip(errflux, sigma_lower, sigma_upper, maxiters, masked=True)
    mask = sigma_flux.mask & sigma_errflux.mask # combine the two masks
    return mask

def my_sin(time, offset, amplitude, period, phase):
    return np.sin(2*np.pi * time / period + phase) * amplitude + offset

def pick_period(freq, power, peak_num=1):
    
    periods = np.array(1/freq)
    sort_peaks = np.flip(np.argsort(power))
    periods  = periods[sort_peaks]
    
    # select peak and remove 1-day period signals
    peak_ind = 0
    try:
        period = periods[peak_ind]
    except:
        print('No peaks found on Periodogram')
        period = 1.0 # day
    while (abs(period - 1.0) < 0.1) or (abs(period - 0.50) < 0.1) or (abs(period - 0.33) < 0.1)  :
        # print('Removed 1-day derived period signal')
        peak_ind += 1 
        period = periods[peak_ind]
        
    return period


def lombscargle_periodogram(time, flux, error, dt=0, min_period=0.1, 
                            max_period=10, peak_points=10, height=0, 
                            peak_ind=0, show=True, xlim=(0,1)):
    '''
    this function determines the peak period of a given light curve, allows
    one to choose which peak to select, then plots the periodogram and the
    folded light curve
    
    Parameters
    ----------
    time : array of float
        contains time data for the light curve
    flux : array of float
        contains flux data for the light curve
    error : array of float
        contains error data for the light curve
    dt : float
        time shift [default = 0]
    min_period : float
        minimum period to investigate [default = 0.1 days]
    max_period : float
        maximum period to investigate [default = 4 days]
    peak_points : int
        number of points around peaks [default = 10]
    height : float
        minimum height to consider peaks [default = 0]
    peak_ind : ind
        choose a peak number, maximum is default [peak_ind = 0]
    plot : bool
        plot a periodogram and the folded lightcurve with best fit sinusoid
    xlim : tuple
        x-limits of the folded plot [default = (0, 1)]

    Returns
    -------
    Pb : tuple
        contains the best fit parameters for the sine wave (amplitude,
        period, phase)
    residuals : array of float
        contains the residuals between the best fit model and the data
    '''
    time_fixed = time - dt
    # create the periodogram
    model = LombScargle(time_fixed, flux, error)
    
    frequencies, power = model.autopower(minimum_frequency=(1./max_period), 
                                         maximum_frequency=(1/min_period), 
                                         samples_per_peak=peak_points)
     
    # convert to periods
    periods = 1/frequencies
    # identify and extract peaks
    inds, peaks = find_peaks(power, height=height)
    peaks = peaks['peak_heights']
    sort_peaks = np.argsort(peaks)
    inds  = inds[sort_peaks]
    peaks = peaks[sort_peaks]
    
    # select peak and remove 1-day period signals
    try:
        period = periods[inds[-1 - peak_ind]]
    except:
        print('No peaks found on Periodogram')
        period = 1.0 # day
    while (abs(period - 1.0) < 0.1) or (abs(period - 0.50) < 0.1) or (abs(period - 0.33) < 0.1)  :
        print('Removed 1-day derived period signal')
        peak_ind += 1 
        period = periods[inds[-1 - peak_ind]]
        
    # false_alarm = model.false_alarm_probability(peaks[peak_ind]) # always yields 1.0
    # fit the sinusoid
    flux_fit = model.model(time_fixed, 1/period)
    residuals = flux - flux_fit
    t0, t1, t2 = model.model_parameters(1/period)
    # convert theta parameters to amplitude and phase
    amplitude = np.hypot(t1, t2)  * (-1)
    phase = -np.arctan(t1 / t2) + np.pi

    
    Pb = [t0, amplitude, period, phase] # WARNING: change Apr 14 12:44h

        
    # plot the periodogram and folded light curve
    if show == True:
        # periodogram
        # fig = plt.figure(figsize=(16, 8))
        # plt.title('Lomb-Scargle Periodogram of Stellar Variations')
        # plt.xlabel('Period [days]')
        # plt.ylabel('Power [-]')
        # plt.plot(periods, power, 'b-')
        # plt.gca().axvline(x=period, color='k', ls=':')
        # plt.show()
        # folded light curve
        plot.plot_folded(time_fixed, flux, error, Pb, flux_fit)
        
        print('-------------------------------------------------')
        # print('%.6f + %.6f sin(2 pi time / %.4f + %.4f)' % Pb)
        print('t0 \t\t\t amplitude \t period (days) \t phase \n{:.4f} \t {:.3f} \t {:.3f} \t {:.2f}'.format(*Pb))
        print('-------------------------------------------------')      
    return Pb, frequencies, power, residuals

# =============================================================================
#                              Dirk's functions
# =============================================================================

def line(time, slope, y_intercept, dt=0):
    '''
    this is simply the a line function for scipy.optimize
    Parameters
    ----------
    time : array of floats
        contains time data for the light curve
    slope : float
        slope of the line
    y_intercept : float
        y-intercept of the line
    dt : float
        time shift (useful if you for example get the same data in MJD and JD)
    Returns
    -------
    trend : array of floats
        line that follows trend = slope * time_fixed + y_intercept
    '''
    time_fixed = time - dt 
    trend = slope * time_fixed + y_intercept
    return trend


def lombscargle_periodogram_dirk(time, flux, error, ax=None, dt=0, min_period=0.1, 
                            max_period=4, peak_points=10, height=0, 
                            peak_ind=0, show=True, xlim=(0,1)):
    '''
    this function determines the peak period of a given light curve, allows
    one to choose which peak to select, then plots the periodogram and the
    folded light curve
    
    Parameters
    ----------
    time : array of float
        contains time data for the light curve
    flux : array of float
        contains flux data for the light curve
    error : array of float
        contains error data for the light curve
    dt : float
        time shift [default = 0]
    min_period : float
        minimum period to investigate [default = 0.1 days]
    max_period : float
        maximum period to investigate [default = 4 days]
    peak_points : int
        number of points around peaks [default = 10]
    height : float
        minimum height to consider peaks [default = 0]
    peak_ind : ind
        choose a peak number, maximum is default [peak_ind = 0]
    plot : bool
        plot a periodogram and the folded lightcurve with best fit sinusoid
    xlim : tuple
        x-limits of the folded plot [default = (0, 1)]
    Returns
    -------
    Pb : tuple
        contains the best fit parameters for the sine wave (amplitude,
        period, phase)
    residuals : array of float
        contains the residuals between the best fit model and the data
    '''
    time_fixed = time - dt
    # create the periodogram
    model = LombScargle(time_fixed, flux, error)
    frequencies, power = model.autopower(minimum_frequency=(1./max_period), 
                                         maximum_frequency=(1/min_period), 
                                         samples_per_peak=peak_points)
    # convert to periods
    periods = 1/frequencies
    # identify and extract peaks
    inds, peaks = find_peaks(power, height=height)
    peaks = peaks['peak_heights']
    sort_peaks = np.argsort(peaks)
    inds  = inds[sort_peaks]
    peaks = peaks[sort_peaks]
    # select peak
    period = periods[inds[-1 - peak_ind]]
    # fit the sinusoid
    flux_fit = model.model(time_fixed, 1/period)
    residuals = flux - flux_fit
    t0, t1, t2 = model.model_parameters(1/period)
    # convert theta parameters to amplitude and phase
    # amplitude = np.hypot(t1, t2) * np.sign(flux_fit[0])
    amplitude = np.hypot(t1, t2)
    phase = -np.arctan(t1 / t2) + np.pi/2
    Pb = (amplitude, period, phase)
    # plot the periodogram and folded light curve
    if show == True:
     #periodogram
        # fig = plt.figure(figsize=(16, 8))
        # plt.title('Lomb-Scargle Periodogram of Stellar Variations')
        # plt.xlabel('Period [days]')
        # plt.ylabel('Power [-]')
        # plt.plot(periods, power, 'b-')
        # plt.gca().axvline(x=period, color='k', ls=':')
        # plt.show()
        # folded light curve
        # GET THE CORRECT "Pb" values by checking the areas
        Pb = plot_folded(time_fixed, flux, error, Pb, flux_fit, dt, xlim, ax=ax)
        print('{:.6f} sin(2 pi time / {:.4f} + {:.4f})'.format(*Pb))
    return Pb, residuals

def plot_folded(time, flux, error, Pb, flux_fit, dt=0, xlim=(0,1), ax=None):
    '''
    this function makes a phase-folded plot of the provided light curve
    
    Parameters
    ----------
    time : array of float
        contains time data for the light curve
    flux : array of float
        contains flux data for the light curve
    error : array of float
        contains error data for the light curve
    amplitude : float
        best fit amplitude for the sine fit
    Pb : tuple, list, array
        contains best fit parameters for sine curve
            amplitude --> of the sine
            period --> of the sine
            phase --> of the sine
    dt : float
        time shift [default = 0]
    xlim : tuple
        xlims of the plot [default = (0, 1)]
    Returns
    -------
    matplotlib.figure()
    '''
    # correct time
    time_fixed = time - dt
    # extract best fit
    amp, prd, phs = Pb
    phase = (time_fixed % prd) / prd
    # fig = ax.figure(figsize=(16, 8))

    sort = np.argsort(phase)
    ###### dario's fork ######
    fun = sines(time_fixed, [amp], [prd], [phs])
    fun2 = sines(time_fixed, [amp], [prd], [phs+np.pi])
    
    
    area1 = sum((fun-flux_fit)**2)
    area2 = sum((fun2-flux_fit)**2)
    

    if ax != None:
        ax.set(ylabel='Flux [-]')
        ax.text(x=0.1, y=0.85, s='Period = %.4f' % prd,
            transform=ax.transAxes)
       
        ax.errorbar(phase[sort], flux[sort], yerr=error, fmt='.', color='k', alpha=0.7)
        ###############
        ax.plot(phase, sines(time_fixed, [amp], [prd], [phs]),
                 'r.', lw=8, label='A = {:.3f}'.format(area1), alpha=0.7)
        # print('Chisquare for phs: {:.2f}'.format(chisquare(flux_fit,sines(time_fixed, [amp], [prd], [phs]))[0]))
    
    
        ax.plot(phase, sines(time_fixed, [amp], [prd], [phs+np.pi]),
             'b.', lw=8, label='A = {:.3f}'.format(area2), alpha=0.7)
        ax.plot(phase, flux_fit,'g.', lw=8, alpha=0.7)
        ax.set_ylim(-0.1,0.1)
    
        ax.set_xlim(xlim)
        ax.legend(ncol=2)
    
    # Update the value of "phs" is the blue area is smaller
    if  area1 == np.min([area1, area2]):
        pass
    else:
        phs += np.pi
        
    Pb = [amp, prd, phs] # CORRECT VALUES
    return Pb

def sines(time, amplitudes, periods, phases, dt=0):
    '''
    function that returns the sum of several sinusoids
    Parameters
    ----------
    time : array of float
        contains time data for the light curve
    amplitudes : list of floats
        amplitudes of the sines
    periods : list of floats
        periods of the sines
    phases : list of floats
        phases of the sines
    dt : float
        time shift (useful if you for example get the same data in MJD and JD)
    Returns
    -------
    trend : array of float
        sum of sines with given input parameters
    '''
    time_fixed = time - dt
    trend = 0
    for amplitude, period, phase in zip(amplitudes, periods, phases):
        sine = amplitude * np.sin(2 * np.pi * time_fixed / period + phase)
        trend += sine
    return trend

# for each period determine what the minimum
def test_harmonic(fund_period, test_periods, pop): 
    try:
        ratios = test_periods / fund_period 
        harmonics = np.round(ratios).astype(np.int) 
        discrepancies = 100 * np.abs(harmonics - ratios) 
        for test_period, harmonic, discrepancy in zip(test_periods, harmonics, discrepancies): 
            print('the discrepancy between %.4f and %.4f is %.2f %% in the %i harmonic' % (fund_period, test_period, discrepancy, harmonic)) 
            if discrepancy < 5:
                pop.append(test_period)
    except:
        pass
    return pop

