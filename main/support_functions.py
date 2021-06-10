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
            # temp = load_dat(file)
            data = ascii.read(file, delimiter='\t')
            df = data.to_pandas()
            df['flux_err'] = df['flux_err'].replace(99.99, np.nan)
            df = df.dropna()
            name = file.split('/')[-1][:-4]
            # hjd, flux, errflux, name = load_dat(file)
            df['gaia_id'] = np.full_like(df['HJD'], name)
            # convert to relative flux
            flux_med = np.nanmedian(df['flux'])
            df['flux'] = df['flux'] / flux_med
            df['flux_err'] = df['flux_err'] / flux_med
            # data =  np.transpose(np.array([hjd, flux, errflux, name_col]))
            # df = pd.DataFrame(data, columns=['HJD', 'flux','errflux','gaia_id'])
    
        if self.band == 'g':
           # data is [hjd, flux, errflux, asas_id]
            df = pd.DataFrame(np.loadtxt(file,usecols=[0,1,2,-1]), columns=['HJD', 'flux','flux_err','asas_id'])
            df.drop(columns=['asas_id']) # remove ASAS_SN_ID 
            df = df.dropna()
            name = int(file.split('/')[-1][:-4])
            df['gaia_id'] = np.full_like(df['HJD'], name)
            
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
        
    

def load_dat(file, y='mag'):
    '''
    Load LC object from .dat file

    Parameters
    ----------
    file : string
        full path for the file.
    y : string
        WHAT COLUMN TO READ FROMN THE FILE. The default is 'mag'.
        The output is always flux, so when 'mag' is is picked, the flux will be 
        calculated from 'mag' values

    Returns
    -------
    hjd : array
        time.
    flux : array
        Relative flux.
    errflux : array
        Relative error flux.
    name : string
        Gaia ID.

    '''
    
    if y =='flux':
        err_type = 'flux_err'
    else:
        err_type = 'mag_err'
        
    data = ascii.read(file, delimiter='\t')
    hjd = data['HJD']
    y_err = data[err_type]
    
    y_data = [] # either 'flux' or 'mag'
    for item in data[y]:
        a  = str(type(item))
        if a == "<class 'numpy.str_'>":
            y_data.append(float(item.split('>')[-1]))
        else:
            y_data.append(item)

    
    name = file.split('/')[-1][:-4]
    
    # Fix "infinite values" for the "errflux" data by
    # replacing them with the median
    med = np.median(y_err)
    if med > 90:
        print('Error bars not valid for plot ', name)
    
    y_err = np.array(y_err)                 # convert to np.array
    med = np.median(y_err[y_err<90])    # compute the median clipping values >90
    y_err[y_err>90] = med               # replace conflicting values with the median
    
    if y =='mag':     
        mag = np.array(y_data)              # renaming for readability
        errmag = np.array(y_err)            # renaming for readability
        med_mag = np.median(mag)
        flux = 10**(-mag/2.5) 
        flux_rel = flux / np.median(flux)
        errflux_rel = flux_rel * errmag / 1.09 # http://slittlefair.staff.shef.ac.uk/teaching/phy217/lectures/stats/L18/index.html#magnitudes

        
    return hjd, flux_rel, errflux_rel, name

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
def lombscargle_periodogram(time, flux, error, dt=0, min_period=0.1, 
                            max_period=10, peak_points=10, height=0, 
                            peak_ind=0, plot=True, xlim=(0,1)):
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
    while (abs(period - 1.0) < 0.1) or (abs(period - 0.50) < 0.1)  :
        print('Removed 1-day period signal')
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
    if plot == True:
        # periodogram
        # fig = plt.figure(figsize=(16, 8))
        # plt.title('Lomb-Scargle Periodogram of Stellar Variations')
        # plt.xlabel('Period [days]')
        # plt.ylabel('Power [-]')
        # plt.plot(periods, power, 'b-')
        # plt.gca().axvline(x=period, color='k', ls=':')
        # plt.show()
        # folded light curve
        plot_folded(time_fixed, flux, error, Pb, flux_fit)
        
        print('-------------------------------------------------')
        # print('%.6f + %.6f sin(2 pi time / %.4f + %.4f)' % Pb)
        print('t0 \t\t\t amplitude \t period (days) \t phase \n{:.4f} \t {:.3f} \t {:.3f} \t {:.2f}'.format(*Pb))
        print('-------------------------------------------------')      
    return Pb, frequencies, power, residuals