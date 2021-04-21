#!/usr/bin/env python3
# -*- coding: utf-8 -*-

###############################################################################
############################# STELLAR VARIATIONS ##############################
###############################################################################

# This module contains all the functions necessary to model the stellar 
# variations from the K2 light curve of V928 Tau



########################
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
        
        # Return a list with ONLY the ASAS_ID for each LC
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
            hjd, flux, errflux, name = load_dat(file)
            name_col = np.full(shape=len(hjd), fill_value=name, dtype='int')
            data =  np.transpose(np.array([hjd, flux, errflux, name_col]))
            df = pd.DataFrame(data, columns=['HJD', 'flux','errflux','gaia_id'])
    
        if self.band == 'g':
           # data is [hjd, flux, errflux, asas_id]
            df = pd.DataFrame(np.loadtxt(file,usecols=[0,1,2,-1]), columns=['HJD', 'flux','errflux','asas_id'])
            df.drop(columns=['asas_id']) # remove ASAS_SN_ID 
            name = int(file.split('/')[-1][:-4])
            df['gaia_id'] = np.full_like(df['HJD'], name)
            
            # convert to relative flux
            flux_med = np.median(df['flux'])
            df['flux'] = df['flux'] / flux_med
            df['errflux'] = df['errflux'] / flux_med
        return df
        
    
    
    
def csv_to_dat(file):
    name = file.split('/')[-1][2:-4] + '.dat'
    
    txt_file = '/home/dario/AstronomyLeiden/FRP/gband/'+name
    with open(file, 'r') as inp, open(txt_file, 'w') as out:
        for k, line in enumerate(inp):
            if k == 0:
                line = line.replace(';', ' ')
                line = line.replace('jd', 'HJD')
                out.write(line)
            else:
                line = line.replace(';', ' ')
                out.write(line)

    return name

def loadmultdata(rand_files, y='flux', plot=True):
    '''
    Function to visualize N random LCs in files .dat generated from ASAS-API
    It can either show magnitudes or fluxes with the "y" parameter
    
    Parameters
    ----------
    rand_files : array(floats)
        list of files
    y : string
        'flux' or 'mag'. The default is 'flux'.
    plot : Boolean
        Choose whether to plot the LCs. The default is True.
    
    Returns
    -------
    hjd : array(floats)
        
    flux : array(floats)
        
    errflux : array(floats)
      
    names : array(strings)
        names of files
    '''
    
    if y =='flux':
        err_type = 'flux_err'
    else:
        err_type = 'mag_err'
    
    
    hjd, flux, errflux, names = [[] for i in range(4)]
    for i in rand_files:
        data = ascii.read(file_list[i])
        hjd.append(data['HJD'])
        errflux.append(data[err_type])
        
        temp_list = []
        for item in data[y]:
            a  = str(type(item))
            if a == "<class 'numpy.str_'>":
                temp_list.append(float(item.split('>')[-1]))
            else:
                temp_list.append(item)
        
        flux.append(temp_list)
        
        name = file_list[i].split('/')[-1][2:-4]
        names.append(name)
    
    # Fix "infinite values" for the "errflux" data by
    # replacing them with the median
    for j in range(len(errflux)):
        med = np.median(errflux[j])
        if med > 90:
            print(f'Error bars not valid for plot {names[j]}')
            
        for count, k in enumerate(errflux[j]):
            if k > 90: # the wrong values have an assigned "99.99"
                if med < 90:
                    errflux[j][count] = np.median(errflux[j])
                else:
                    errflux[j][count] = 0.0 # no value
    
    if plot == True:
        # Plot LCs
        fig, ax = plt.subplots(len(hjd), 1, figsize=(15,16))
           
        # Make common labels
        fs = 15 # fontsize
        fig.text(0.5, 0.04, 'HJD', ha='center', fontsize=fs)
        ylabel = y + '[]'
        fig.text(0.07, 0.5, s=ylabel, va='center', rotation='vertical', fontsize=fs)
        # Plot all the datasets
        for j in range(len(hjd)):
            ax[j].plot(hjd[j], np.array(flux[j], dtype='float'), '.k')
            ax[j].errorbar(hjd[j], np.array(flux[j], dtype='float'), errflux[j],fmt='k.', ecolor='.7', markersize=4)
            ax[j].text(0.87, 0.05, names[j],transform=ax[j].transAxes,fontsize=15)
             
        # plt.savefig('LC1.png', dpi=100, bbox_inches='tight')
        # plt.savefig('/home/dario/AstronomyLeiden/FRP/MRP-ASAS-Rings/name.png', dpi=100, bbox_inches='tight')
        plt.show()
    return hjd, flux, errflux, names


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
 
    name = file.split('/')[-1][2:-4]
    
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


def lombscargle_periodogram(time, flux, error, dt=0, min_period=0.1, 
                            max_period=20, peak_points=10, height=0, 
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
    period = periods[inds[-1 - peak_ind]]
    while abs(period - 1.0) < 0.1:
        print('Removed 1-day period signal')
        peak_ind += 1 
        period = periods[inds[-1 - peak_ind]]
    
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
        plot_folded(time_fixed, flux, error, Pb, flux_fit, dt, xlim)
        
        print('-------------------------------------------------')
        # print('%.6f + %.6f sin(2 pi time / %.4f + %.4f)' % Pb)
        print('t0 \t\t\t amplitude \t period (days) \t phase \n{:.4f} \t {:.3f} \t {:.3f} \t {:.2f}'.format(*Pb))
        print('-------------------------------------------------')      
    return Pb, frequencies, power, residuals

def plot_folded(time, flux, error, Pb, flux_fit, dt=0, xlim=(0,1)):
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
    t0, amp, prd, phs = Pb
    # t0, t1, t2, prd = Pb
    phase = (time_fixed % prd) / prd
    fig = plt.figure(figsize=(16, 8))
    plt.xlabel('Phase [-]')
    plt.ylabel('Normalised Flux [-]')
    plt.title('Folded Light Curve [Period = %.4f]' % prd)
    sort = np.argsort(phase)
    ###############################
    fun = t0 + amp * np.sin(2*np.pi*time_fixed/prd + phs)
    fun2 = t0 - amp * np.sin(2*np.pi*time_fixed/prd + phs)
    
    area1 = sum((fun[sort]-flux_fit[sort])**2)
    area2 = sum((fun2[sort]-flux_fit[sort])**2)

    print('-------------------------------------------------')
    print('area:  (+) {:.2f} \t (-) {:.2f}'.format(area1, area2))
    print('-------------------------------------------------')
    if abs(area1) > abs(area2):
        print('Flipping amplitude sign!')
        amp = -amp
    ###############################
    plt.errorbar(phase[sort], flux[sort], yerr=error, fmt='.', color='k')
    # plt.plot(phase[sort], sines(time_fixed, [amp], [prd], [phs])[sort],
    #           'bx', label='sines')
    
    plt.plot(phase,t0 + amp * np.sin(2*np.pi*time_fixed/prd + phs),
              'rx', label=r'$t_0 + A \sin{\frac{2\pi}{T}t + \phi}$')
    plt.plot(phase,t0 - amp * np.sin(2*np.pi*time_fixed/prd + phs),
              'bx', label=r'$t_0 - A \sin{\frac{2\pi}{T}t + \phi}$')
    
    # fun = t0 + (t1 * np.sin(2*np.pi*time_fixed/prd)) + (t2 * np.cos(2*np.pi*time_fixed/prd))
    # plt.plot(phase, fun, 'rx', label=r'$C + A \sin{\frac{2\pi}{T}t} + B \cos{\frac{2\pi}{T}t}$')
              

    
    plt.plot(phase[sort], flux_fit[sort],'go', label='Fitted flux')
    plt.legend()
    plt.xlim(xlim)
    plt.show()
    return None

def plot_lc_clipped(clean_data, clipped_data, file, ax=None, **kwargs):
    '''
    Plot the clipped lightcurve (black) and the removed points (red)
    along with a (green) baseline

    Parameters
    ----------
    clean_data : np.array
        array of clipped arrays with time, flux, errflux
    clipped_data : np.array
        array of clipped values with time, flux, errflux
    file : string
        name of file (for title)
    size : tuple
        Figure size. The default is (20,6).

    Returns
    -------
    ax: axes
        axes to be plotted onto fig, ax = plt.subplots()

    '''
    ax = ax or plt.gca()
    ax.errorbar(clean_data[0], clean_data[1], clean_data[2], fmt='.k', alpha=0.8, label='Error-clipped')
    ax.errorbar(clipped_data[0], clipped_data[1],clipped_data[2], fmt='xr', alpha=0.8, label='Clipped')
    ax.set_ylabel('Normalised Flux')
    ax.set_title('ASAS_ID = ' + file[2:-4])
    ax.legend()
    # ax.set_xlabel('Time [JD]')
    flux_sv = clean_data[1]
    ax.axhspan(np.median(flux_sv)-np.std(flux_sv), np.median(flux_sv)+np.std(flux_sv), alpha=0.4, color='lightgreen')
    return(ax)


def hist_errbars(errflux_sv, errflux, count, ax=None):
    '''
    Compare the distribution of error bars among the clipped and raw data

    Parameters
    ----------
    errflux_sv : np.array()
        clipped data
    errflux : np.array()
        raw data
    size : tuple
        Figure size. The default is (9,6).

    Returns
    -------
    fig : plt.fig()
        Figure
    '''
    ax = ax or plt.gca()
    # ax.set_title('Number of clippings: {:}'.format(count))
    ax.hist(errflux_sv, bins='auto', color='black', alpha=0.8, label='Clipped')
    ax.hist(errflux, bins='auto', color='red', alpha=0.4, label='Raw')
    ax.set_ylabel('Number of points')
    ax.set_xlabel('Flux Errorbar [-]')
    ax.legend()
    return(ax)

def plot_periodogram(data, ax=None):
    freq_sv, power_sv, freq, power = data[0], data[1], data[2], data[3]
    
    ax = ax or plt.gca()
    ax.plot(1/freq_sv, power_sv, color='black', alpha=0.8, label='Clipped')
    ax.axhline(max(power_sv), label='Clipped max(power)', color='black', ls='--', alpha=0.4)
    ax.plot(1/freq, power, alpha=0.5, label='Raw', color='red')
    ax.axhline(max(power), label='Raw max(power)', color='red', ls='--', alpha=0.4)
    
    ax.set_ylim(0,None)
    ax.set_xlim(1/freq[-1], 1/freq[0])
    ax.set_ylabel('Power')
    # ax.set_title('ASAS_ID = ' + file[2:-4])
    ax.legend(fontsize=10)
    
    ax.set_xlabel('Period')
    return ax
#######################
#%% MODEL FUNCTIONS %%#
#######################

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

