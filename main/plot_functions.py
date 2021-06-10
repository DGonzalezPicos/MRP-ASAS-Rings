"""
Created on Tue Jun  8 12:44:40 2021
@author: dario

Title: Script with plotting functions
"""
import numpy as np
import matplotlib.pyplot as plt

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
    
    
    # Fit a sin function to get better parameters
    # sci_fit = curve_fit(my_sin, time_fixed, flux_fit, p0=Pb)
    # scipy_fit = my_sin(time_fixed, *sci_fit[0])
    # plt.plot(phase, scipy_fit, '.', label='scipy_fit', alpha=0.7, markersize=14)

    
    plt.plot(phase[sort], flux_fit[sort],'go', label='Fitted flux')
    plt.legend()
    plt.xlim(xlim)
    # plt.show()
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

def plot_periodogram(data, ax=None, clipped=True, raw=True):
    freq_sv, power_sv, freq, power = data[0], data[1], data[2], data[3]
    # freq_sv, power_sv = data[0], data[1]
    
    ax = ax or plt.gca()
    if clipped == True:
        ax.plot(1/freq_sv, power_sv, color='black', alpha=0.8, label='Clipped')
        ax.axhline(max(power_sv), color='black', ls='--', alpha=0.4)
        
    if raw == True:
        ax.plot(1/freq, power, alpha=0.99, label='Periodogram', color='orange')
        ax.axhline(y=0.10, color='red', ls='--', alpha=0.2)

        
        # ax.axhline(max(power), color='red', ls='--', alpha=1.2)
    
    ax.set_ylim(0,None)
    ax.set_xlim(1/freq[-1], 1/freq[0])
    ax.set_ylabel('Power [-]')
    # ax.set_title('ASAS_ID = ' + file[2:-4])
    # ax.legend(fontsize=12, frameon=False)
    ax.set_xlabel('Period [days]')
    return ax