"""
Created on Tue Jun  8 12:44:40 2021
@author: dario

Title: Script with plotting functions
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
mpl.rcParams['font.size'] = 14

def plot_lc(df, name=None, y='mag', ax=None, label=None, errorbars=True, title=False, invert=False, **kwargs):
    
    time_corr = df['HJD'] - 2450000 # transform HDJ date
    y_err = y + '_err'
    
    ax = ax or plt.gca()
    
    if errorbars==True:
        ax.errorbar(time_corr, df[y], df[y_err], fmt='.', alpha=0.6, label=label)
    else:
        ax.plot(time_corr, df[y],'.k', alpha=0.6, label=label, **kwargs)
        
        # bars = [0, int(np.round(len(df) / 2.)), -1]
        # ax.errorbar(time_corr.values[bars], df[y].values[bars], df[y_err].values[bars], fmt='.b', alpha=0.2)
    
    if title == True:
        ax.set_title('GAIA_ID = ' + name)
  
    ax.set_ylabel('Magnitude')
    if invert==True:
        ax.invert_yaxis()
    ax.set_xlabel('HJD - 2450000')
    mean_mag = np.nanmean(df['mag'])
    mean_mag_str = 'Mean V-mag = {:.1f}'.format(mean_mag)
        
    # if mean_mag > 17:
    #     color_mag = 'red'
    # elif (mean_mag < 17) and (mean_mag > 15):
    #     color_mag = 'darkorange'
    # else:
    #     color_mag = 'darkgreen'
    # ax.text(x=0.02, y=0.89, s=mean_mag_str, transform=ax.transAxes, color=color_mag)
        
    if y == 'flux':
        ax.set_ylabel('Flux')
    #     ymin_auto, ymax_auto = ax.yaxis.get_data_interval()
    #     ax.set_ylim(ymin=np.max([ymin_auto, 0.0]))
    #     ax.set_ylim(ymax=np.min([ymax_auto, 2.0]))
    return ax

def plot_full(df, gaia_id, y='mag', ax=None):
    
    y_err = y + '_err'
    ax = ax or plt.gca() # function embedded
    markers = {'V':'.', 'g':'^','B':'v','I':'s','R':'d'}
    colors = {'V':['goldenrod','lightcoral','chocolate','red'], 
              'g':['green', 'teal','aqua', 'indigo','fuchsia'],
              'B':['royalblue','navy','slateblue','indigo'],
              'I':['purple','orchid','violet','pink'],
              'R':['red','tomato','coral','peru']}
    for filter in df['Filter'].unique():
        cams = df[df['Filter']==filter]['camera'].unique()
        # colors = ['brown','darkorange','chocolate']
        print('Plotting {:} in {:}-band with cameras {:}'.format(gaia_id, filter, cams))
        for i, cam in enumerate(cams):
            _ = plot_cam(df, cam=cam, y=y, ax=ax, 
                         # label=cam, color=colors[filter][i],
                         label=cam,
                         marker=markers[filter])
    if y == 'mag':
        ax.invert_yaxis()  
        ax.set_ylabel('Magnitude')
    if y == 'flux':
        ax.set_ylabel('Relative flux')
        
    ax.set_title('GAIA_ID = ' + str(gaia_id))
    
    ax.set_xlabel('HJD - 2450000')
    
    ax.legend(loc=(1.02,0.0))
    plt.show()
    return ax

## Plot light curve DataFrame object with its different cameras
def plot_cam(df, cam, y, ax=None, **kwargs):
            lc = df[df['camera'] == cam] # select camera
            # name = str(int(df['gaia_id'].unique())) 
            time_corr = lc['HJD'] - 2450000 # transform HDJ date
            y_err = y + '_err'
            ax.errorbar(time_corr, lc[y], lc[y_err], fmt='.', alpha=0.4, **kwargs)
            return ax
        
## Embedded function for selected V-band LCs with Gaia_ID (with different cams)
def plot_gaiaID(gaia_id, lco, ax=None, save=False):
    
        file = os.path.join(lco.path, str(gaia_id)+ '.dat')
        df = lco.data(file)
        df['mag'] = df['mag'].astype('float')
        
        ax = ax or plt.gca() # function embedded
        # fig, ax = plt.subplots(1, figsize=(14,7))
        # _ = plot.plot_lc(df, name=str(gaia_id), ax=ax)
    
        
        cams = df['camera'].unique()
        colors = ['brown','darkorange','chocolate']
        print('Plotting {:} with cameras {:}'.format(gaia_id, cams))
        for i, cam in enumerate(cams):
            _ = plot_cam(df, cam=cam, y='mag', ax=ax, label=cam, color=colors[i])
        
        ###########################
        add_gband = False 
        ###########################
        
        if add_gband == True:
            g_path = '/home/dario/AstronomyLeiden/FRP/gband_api_fixed/'
            gband = supp.LightCurveSet(g_path, 'dat', 'v') # read as v-band because new g-data is like old v-band
            file_g = os.path.join(g_path, str(gaia_id)+ '.dat')
            df_g = gband.data(file)
            
            g_cams = df_g['camera'].unique()
            colors_g = ['green', 'teal','navy','plum']
            print(g_cams)   
            for i, cam in enumerate(g_cams):
                _ = plot_cam(df_g, cam=cam, y='mag', ax=ax, label=cam, color=colors_g[i], marker='^')
            
        ax.invert_yaxis()   
        ax.set_title('GAIA_ID = ' + str(gaia_id))
        ax.set_ylabel('Magnitude')
        ax.set_xlabel('HJD - 2450000')
        
        ax.legend()
        plt.show()
        if save == True:
            outpath = '/home/dario/AstronomyLeiden/FRP/ring_candidates/' + str(gaia_id)+'.png'
            plt.savefig(outpath, dpi=200, bbox_inches='tight') ## UNCOMMENT TO SAVE FIG
            print('Saved {:}'.format(outpath))
        return ax
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