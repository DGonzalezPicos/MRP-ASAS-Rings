"""
Created on Fri Jun 25 20:16:55 2021
@author: dario

Title: 
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import support_functions as supp
import plot_functions as plot
import os

info = pd.read_csv('short_list_df.csv')

# =============================================================================
#                   3016935867962755968
# =============================================================================
gaia_id = 3016935867962755968

lc = supp.LightCurveObject(gaia_id, 'V')    

# import manually downloaded g-band from SkyPatrol website
lc_sp = pd.read_csv('3016935867962755968_skypatrol.csv')
lc_sp['mag_err'] = lc_sp['mag_err'].replace(99.99, np.nan)
lc_sp = lc_sp.dropna()
lc_sp['mag'] = lc_sp['mag'].astype('float')
lc_sp = lc_sp.rename(columns={'flux(mJy)':'flux', 'Camera':'camera'})
 # convert to relative flux
flux_med = np.nanmedian(lc_sp['flux'])
lc_sp['flux'] = lc_sp['flux'] / flux_med
lc_sp['flux_err'] = lc_sp['flux_err'] / (1.09*flux_med)
# pick the g-band only
# lc_g = lc_skypatrol[lc_skypatrol['Filter']=='g']



fig, ax = plt.subplots(1, figsize=(16,5))
_ = plot.plot_full(lc_sp, str(gaia_id), y='flux', ax=ax)
# plt.savefig('raw_full_lc.png', dpi=200, bbox_inches='tight')

#%% 
# =============================================================================
#                       Correct offset between bands
# =============================================================================
# take the V-band as the good one and correct the g-band

def band_offset(lc, y='mag'):
    '''
    Given a light curve DataFrame with data in V- and g-band returns a corrected
    light curve for the band-offset.

    Parameters
    ----------
    lc : DataFrame
        full light curve

    Returns
    -------
    lc : DataFrame
        corrected light curve.

    '''
    vband = lc[lc['Filter']=='V'].sort_values('HJD')
    gband = lc[lc['Filter']=='g'].sort_values('HJD')
    
    vband_overlap = vband[vband['HJD']>gband['HJD'].iloc[0]]
    gband_overlap = gband[gband['HJD']<vband['HJD'].iloc[-1]]
    
    assert len(vband_overlap) > 10, 'Error: NO OVERLAP'
    
    offset = np.nanmedian(vband_overlap[y]) - np.nanmedian(gband_overlap[y])
    gband_corr = gband[y] + offset
    lc.loc[lc_sp['Filter']=='g',y] = gband_corr
    
    fig, ax = plt.subplots(1, figsize=(16,5))
    _ = plot.plot_full(lc, str(gaia_id), y=y)
    # plt.savefig('corrected_full_lc.png', dpi=200, bbox_inches='tight')
    return lc
lc_sp_corr = band_offset(lc_sp, y='flux')


#%%
# =============================================================================
#                       TESS DATA
# =============================================================================
tess = np.loadtxt('tess_3016935867962755968.txt', dtype=float, delimiter=',')
df_tess = pd.DataFrame(tess, columns=['BTJD','NORM_SAP_FLUX',
                                      'NORM_SAP_BKG'])
df_tess['HJD'] = df_tess['BTJD'] + 7000

fig, (ax1,ax2) = plt.subplots(2,1, figsize=(16,9))

_ = plot.plot_full(lc_sp_corr, str(gaia_id), y='flux', ax=ax1)
ax1.axvspan(df_tess['HJD'].iloc[0], df_tess['HJD'].iloc[-1], 
            color='blue', alpha=0.2, label='TESS')

ax2.plot(df_tess['HJD'], df_tess['NORM_SAP_FLUX'],'.', color='blue', alpha=0.6, label='TESS')
ax2.set_ylabel('Relative flux')
ax2.set_xlabel('Time HJD - 245000')
ax2.legend()

plt.savefig('full_lc_corrected_TESS.png', dpi=200, bbox_inches='tight')















