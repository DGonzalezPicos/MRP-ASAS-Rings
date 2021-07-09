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
import os, glob
from astropy.io import ascii 

info = pd.read_csv('short_list_df.csv')

# =============================================================================
#                   3016935867962755968
# =============================================================================
gaia_id = 3016935867962755968

lc = supp.LightCurveObject(gaia_id, 'V')    

# import manually downloaded g-band from SkyPatrol website
lc_sp = pd.read_csv('gdr2_30169/3016935867962755968_skypatrol.csv')
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

ax.set_ylim(0,3)

# Plot light curve without correction
# _ = plot.plot_full(lc_sp, str(gaia_id), y='flux', ax=ax)
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
    lc.loc[lc['Filter']=='g',y] = gband_corr
    
    fig, ax = plt.subplots(1, figsize=(16,5))
    _ = plot.plot_full(lc, str(gaia_id), y=y)
    # plt.savefig('corrected_full_lc.png', dpi=200, bbox_inches='tight')
    return lc

# Plot corrected light curve
lc_sp_corr = band_offset(lc_sp, y='flux')


#%%
# =============================================================================
#                       TESS DATA
# =============================================================================
tess = np.loadtxt('gdr2_30169/tess_3016935867962755968.txt', dtype=float, delimiter=',')
df_tess = pd.DataFrame(tess, columns=['BTJD','NORM_SAP_FLUX',
                                      'NORM_SAP_BKG'])

df_tess = df_tess.rename(columns={'NORM_SAP_FLUX':'flux'})
df_tess['HJD'] = df_tess['BTJD'] + 7000

fig, (ax1,ax2) = plt.subplots(2,1, figsize=(16,9))

ax1.set_ylim(0,3)
_ = plot.plot_full(lc_sp_corr, str(gaia_id), y='flux', ax=ax1)
ax1.axvspan(df_tess['HJD'].iloc[0], df_tess['HJD'].iloc[-1], 
            color='blue', alpha=0.2, label='TESS')

ax2.plot(df_tess['HJD'], df_tess['flux'],'.', color='blue', alpha=0.6, label='TESS')
ax2.set_ylabel('Relative flux')
ax2.set_xlabel('Time HJD - 245000')
ax2.legend()

# plt.savefig('full_lc_corrected_TESS.png', dpi=200, bbox_inches='tight')

#%%
# =============================================================================
#                           AAVSO
# =============================================================================

aavso_file = 'gdr2_30169/aavsodata_60d9941262470.txt'

df_aavso = pd.read_csv(aavso_file)
df_aavso = df_aavso.drop(columns='HJD')
df_aavso = df_aavso.rename(columns={'JD':'HJD', 'Magnitude':'mag', 'Uncertainty':'mag_err',
                              'Band':'Filter','Observer Code':'camera'})

# df_aavso['mag'] = df_aavso['mag'].astype('float')

temp = []
for item in df_aavso['mag'].values:
    # a  = str(type(item))
    temp.append(float(item.split('<')[-1]))

df_aavso['mag'] = temp
# df_aavso[['JD', 'Magnitude', 'Uncertainty','Band','Measurement Method']]
df_aavso.keys()
df_aavso['Filter'].unique()
df_aavso['camera'].unique()

df_aavso_V = df_aavso[df_aavso['Filter']=='V']
df_aavso_V['camera'].unique()

fig, ax = plt.subplots(2, figsize=(16,9)) 
# ax.set_ylim(0,3)
_ = plot.plot_full(df_aavso, str(gaia_id), y='mag', ax=ax[0])
ax[0].set_xlabel(None)
ax[0].legend(loc=(1.01, -0.6))
# plt.savefig('raw_full_lc.png', dpi=200, bbox_inches='tight')

# Fix offsets in AAVSO
# offset matrix camera x filter (7x4)
offset_matrix = pd.DataFrame(columns=df_aavso['Filter'].unique())
# offset_matrix['V']['HGUA'] = 0.222

def aavso_corr(df):
    df_corr = df.copy() # new copy
    df_main = df[(df['camera']=='HGUA') & (df['Filter']=='V')] # subset with more points
    baseline = np.nanmedian(df_main['mag']) # baseline as mean mag of main subset
    for f in df['Filter'].unique():
        list_temp = []
        for cam in df['camera'].unique():
            df_temp = df[(df['camera']==cam) & (df['Filter']==f)]
            
            if len(df_temp) < 1: # empty subset
                offset = np.nan
            else: # compute offset to be SUBTRACTED
                offset = np.nanmedian(df_temp['mag']) - baseline
            list_temp.append(offset)
            corr = df_temp['mag'] - offset
            df_corr.loc[(df_corr['Filter']==f) &  (df_corr['camera']==cam),'mag'] = corr
            
            print('{:}-filter and {:}-camera: {:.2f}'.format(f, cam, offset))
        offset_matrix[f] = list_temp
    return df_corr

df_aavso_corr = aavso_corr(df_aavso)


                  
# _ = plot.plot_lc(df_aavso_main, name=str(gaia_id))    

# _ = plot.plot_full(df_aavso, str(gaia_id), y='mag', ax=ax[0])
_ = plot.plot_full(df_aavso_corr, str(gaia_id), y='mag', ax=ax[1])
ax[1].set_title(None)

text_box = dict(facecolor='gray', alpha=0.1)
ax[1].text(x=0.025, y=0.85, s='Corrected', transform=ax[1].transAxes,
           bbox=text_box)
ax[1].get_legend().remove()

plt.show()
# plt.savefig('aavso_corrected', dpi=200, bbox_inches='tight')
#%%
def plot_corrected(df, df_corr,ax=None):

    # ax.set_ylim(0,3)
    _ = plot.plot_full(df, str(gaia_id), y='mag', ax=ax[0])
    ax[0].set_xlabel(None)
    ax[0].legend(loc=(1.01, -0.6))
    
    _ = plot.plot_full(df_corr, str(gaia_id), y='mag', ax=ax[1])
    ax[1].set_title(None)
    
    text_box = dict(facecolor='gray', alpha=0.1)
    ax[1].text(x=0.025, y=0.85, s='Corrected', transform=ax[1].transAxes,
               bbox=text_box)
    ax[1].get_legend().remove()
    
    # plt.show()
    return ax

fig, ax = plt.subplots(2, figsize=(16,9)) 


_ = plot_corrected(df_aavso, df_aavso_corr,ax=ax)
plt.show()    
    
#%%
def flux_from_mag(df):
    flux = 10**(-df['mag']/(2.5*np.nanmedian(df['mag']))) 
    flux_err = flux * df['mag_err'] / 1.09
    flux_norm = flux / np.nanmedian(flux)
    flux_err_norm = flux_err / np.nanmedian(flux)
    df['flux'], df['flux_err'] = flux_norm, flux_err_norm
    return df
#%%
df_aavso = flux_from_mag(df_aavso)
fig, ax = plt.subplots(1, figsize=(16,5))
_ = plot.plot_full(df_aavso, str(gaia_id), y='flux', ax=ax)
plt.show()


#%%
# =============================================================================
#                               ATLAS
# =============================================================================
from astropy.stats import sigma_clip
atlas_file = 'gdr2_30169/atlas_forced_photometry.txt'

df_atlas = pd.read_csv(atlas_file, delimiter = " ")
df_atlas = df_atlas.dropna(axis=1, how='all') # drop useless columns
df_atlas.keys()

# rename to standard notation
df_atlas = df_atlas.rename(columns={'###MJD':'HJD','Unnamed: 2':'mag', 
                                    'Unnamed: 4': 'mag_err', 'm':'camera'})
df_atlas = flux_from_mag(df_atlas) # generate flux and flux_err columns
df_atlas = df_atlas.dropna(subset=['mag_err']) # fix nan values 
df_atlas = df_atlas[df_atlas['mag'] > 0] # fix negative values

df_atlas['mag_err'] = sigma_clip(df_atlas['mag_err'], sigma=2)

def pd_sigma(df, y='mag_err', n=1):
    median = np.nanmedian(df[y])
    std = np.std(df[y])
    sigma = n*std
    df = df[abs(df[y] - median) < sigma]
    return df
df_atlas = pd_sigma(df_atlas)

df_atlas['HJD'] = df_atlas['HJD'] + 2400000 # fix MJD to HJD

# fig, ax = plt.subplots(1, figsize=(16,5))
# _ = plot.plot_lc(df_atlas, str(gaia_id), y='flux', ax=ax, errorbars=True)

# filter data points observed in bands 'c', 'o' 
# link: https://fallingstar-data.com/forcedphot/queue/75486/
options = ['c','o']
df_atlas = df_atlas.loc[df_atlas['camera'].isin(options)]





def corr_offset(lc_in, y, key,bands):
    '''
    Given a light curve DataFrame with camera or band offsets returns a corrected
    light curve.

    Parameters
    ----------
    lc : DataFrame
        full light curve

    Returns
    -------
    lc : DataFrame
        corrected light curve.

    '''
    lc = lc_in.copy()
    vband = lc[lc[key]==bands[0]].sort_values('HJD')
    gband = lc[lc[key]==bands[1]].sort_values('HJD')
    
    vband_overlap = vband[vband['HJD']>gband['HJD'].iloc[0]]
    gband_overlap = gband[gband['HJD']<vband['HJD'].iloc[-1]]
    
    assert len(vband_overlap) > 10, 'Error: NO OVERLAP'
    
    offset = np.nanmedian(vband_overlap[y]) - np.nanmedian(gband_overlap[y])
    gband_corr = gband[y] + offset
    lc.loc[lc[key]==bands[1],y] = gband_corr
    
    # fig, ax = plt.subplots(1, figsize=(16,5))
    # _ = plot.plot_cam(lc, str(gaia_id), y=y)
    # plt.savefig('corrected_full_lc.png', dpi=200, bbox_inches='tight')
    return lc

df_atlas_corr = corr_offset(df_atlas, y='mag', key='camera', bands=['o','c'])

np.nanmedian(df_atlas[df_atlas['camera']=='o']['mag'])
np.nanmedian(df_atlas[df_atlas['camera']=='c']['mag'])



df_atlas[df_atlas['camera']=='o'].sort_values('HJD')

fig, ax = plt.subplots(2, figsize=(16,9))
for cam in df_atlas_corr['camera'].unique():
    _ = plot.plot_cam(df_atlas, cam, y='mag', ax=ax[0], label=cam)
    _ = plot.plot_cam(df_atlas_corr, cam, y='mag', ax=ax[1])
    
ax[0].set_title('ATLAS')
ax[0].set_ylabel('Magnitude'), ax[1].set_ylabel('Magnitude')
ax[0].invert_yaxis(), ax[1].invert_yaxis()
ax[0].legend()
ax[1].set_xlabel('HJD - 2450000')
text_box = dict(facecolor='gray', alpha=0.1)
ax[1].text(x=0.025, y=0.15, s='Corrected', transform=ax[1].transAxes,
           bbox=text_box)
plt.show()
# plt.savefig('atlas_corrected', dpi=200, bbox_inches='tight')

#%%
# =============================================================================
#                       PTF & ZTF
# =============================================================================

def load_tf(file):
    '''
    Function to load light curve files from PTF and ZTF

    Parameters
    ----------
    file : string
        Name of the file (with correct format)
        example: "gdr2_30169/ztf_V0658Ori_1.csv"

    Returns
    -------
    df : DataFrame
        Light curve df.
    '''
    
    df = pd.read_csv(file)
    survey = file.split('/')[-1][:3]
    
    if survey == 'ztf':
        # rename to standard notation
        df = df.rename(columns={'hjd':'HJD', 'magerr':'mag_err'})
        
    if survey == 'ptf':
        df = df.rename(columns={'obsmjd':'HJD','mag_autocorr':'mag', 
                                    'magerr_auto': 'mag_err'})
        df['HJD'] = df['HJD'] + 2400000 # fix MJD to HJD


    df = df[df['mag'] > 0] # fix negative values
    
    # df_ptf = flux_from_mag(df_ptf) # generate flux and flux_err columns

    return df

# Plot all the files from a given survey: [ptf, ztf]
survey = 'ptf'

survey_files = glob.glob('gdr2_30169/'+survey+'*')
nrows = len(survey_files)
fig, ax = plt.subplots(nrows, figsize=(12,12))
ax[0].set_title(survey)
for i in range(1, nrows+1):
    file = 'gdr2_30169/' + survey + '_V0658Ori_' + str(i) +'.csv'
    print('Loading ', file)
    df = load_tf(file)
    df = flux_from_mag(df)
    _ = plot.plot_lc(df, str(gaia_id), y='mag', ax=ax[i-1], errorbars=True, title=False)
    
### Clip points with large error bars?¿

#%%
# Merge PTF files into one 
# 4 out of the 5 contain useful data (ptf2.csv disregard, only one point at ~ 18.6 mag)

def merge_tf(survey, y='mag', exclude=None, ax=None):
    survey_files = np.array(glob.glob('gdr2_30169/'+survey+'*'))
    exclude_file = 'gdr2_30169/' + survey + '_V0658Ori_' + str(exclude) +'.csv'
    survey_files = survey_files[survey_files != exclude_file]

    zero_file = survey_files[0]
    print('Loading ', zero_file)
    df_tf1 = load_tf(zero_file)
    df_tf1['file_num'] = np.full_like(df_tf1['HJD'], zero_file.split('.')[-2][-1], dtype=int)
    df_tf_main = df_tf1.copy()
    for file in survey_files[1:]:
        print('Loading ', file)
        df = load_tf(file)
        file_num = file.split('.')[-2][-1]
        df['file_num'] = np.full_like(df['HJD'], file_num, dtype=int)
        df_tf_main = df_tf_main.append(df, ignore_index=True)
        
    df_tf_main['file_num'].unique()
    
    files = np.array(df_tf_main['file_num'].unique(), dtype=int)
    
    
    #
    ax = ax or plt.gca()
    for file in files:
        df_temp = df_tf_main[df_tf_main['file_num']==file]
        ax.set_title(survey.upper())
        _ = plot.plot_lc(df_temp, name=None, y=y, ax=ax,
                         errorbars=True, title=False, label=str(file))
        

    return df_tf_main, ax

fig, ax = plt.subplots(1, figsize=(16,5))
df_ptf, ax = merge_tf(survey='ptf', exclude=None, ax=ax)
ax.invert_yaxis()
ax.legend(loc=(1.01, 0))

plt.show()
# plt.savefig('ztf_combined_raw', dpi=200, bbox_inches='tight')



#%% merge ZTF data

df_ztf, _ = merge_tf(survey='ztf', exclude=None, ax=ax)
# Divide in two groups: files 1+3 and 2+4

df_ztf_1 = df_ztf[(df_ztf['file_num']==1) | (df_ztf['file_num']==3)]
df_ztf_2 = df_ztf[(df_ztf['file_num']==2) | (df_ztf['file_num']==4)]

np.nanmedian(df_ztf_1['mag'])
np.nanmedian(df_ztf_2['mag'])
offset = np.nanmedian(df_ztf_1['mag']) -  np.nanmedian(df_ztf_2['mag'])
df_ztf.loc[(df_ztf['file_num']==2) | (df_ztf['file_num']==4), 'mag'] = df_ztf_2['mag'] + offset

fig, ax = plt.subplots(1, figsize=(16,5))

files = np.array(df_ztf['file_num'].unique(), dtype=int)

for file in files:
    df_temp = df_ztf[df_ztf['file_num']==file]
    _ = plot.plot_lc(df_temp, y='mag', ax=ax,
                     errorbars=True, title=False, label=str(file))
    
ax.set_title('ZTF')
ax.invert_yaxis()    
ax.legend(loc=(1.01, 0))
text_box = dict(facecolor='gray', alpha=0.1)
ax.text(x=0.025, y=0.85, s='Corrected', transform=ax.transAxes,
           bbox=text_box)
plt.show()
# plt.savefig('ztf_combined_corrected', dpi=200, bbox_inches='tight')
#%%
# =============================================================================
#                               CRTS
# =============================================================================

crts_file = './gdr2_30169/result_web_fileZylnk1.csv'
df_crts = pd.read_csv(crts_file)
df_crts = df_crts.rename(columns={'MJD':'HJD','Mag':'mag', 
                            'Magerr': 'mag_err'})
df_crts['HJD'] = df_crts['HJD'] + 2400000 # fix MJD to HJD
df_crts = flux_from_mag(df_crts)


fig, ax = plt.subplots(1, figsize=(12,5))
_ = plot.plot_lc(df_crts, str(gaia_id), y='mag', ax=ax, errorbars=True, title=False)

#---> Clean data with no offsets

#%% 
# =============================================================================
#                       PLOT ALL DATA
# =============================================================================


fig, ax = plt.subplots(1, figsize=(16,5))



lc_sp_corr = band_offset(lc_sp, y='mag')
_ = plot.plot_lc(lc_sp_corr, y='mag', ax=ax, errorbars=True, label='ASAS_SN', invert=True)
_ = plot.plot_lc(df_atlas_corr, y='mag', ax=ax,errorbars=True, label='ATLAS')
_ = plot.plot_lc(df_crts, y='mag', ax=ax,errorbars=True, label='CRTS')
_ = plot.plot_lc(df_aavso_corr, y='mag', ax=ax,errorbars=True, label='AAVSO')
_ = plot.plot_lc(df_ptf, y='mag', ax=ax,errorbars=True, label='PTF')
_ = plot.plot_lc(df_ztf, y='mag', ax=ax,errorbars=True, label='ZTF')
# ax.plot(df_tess['HJD'], df_tess['flux'], '.', label='TESS')

ax.set_ylim(17.5,13)
ax.legend()
plt.show()

#%% save light curves to csv to read them in the notebook
lc_sp_corr.to_csv('asas_sn_corr.csv', index=False)
df_atlas_corr.to_csv('atlas_corr.csv', index=False)
df_crts.to_csv('crts_corr.csv', index=False)
df_aavso_corr.to_csv('aavso_corr.csv', index=False)
df_ptf.to_csv('ptf_corr.csv', index=False)
df_ztf.to_csv('ztf_corr.csv', index=False)

#%%
# =============================================================================
#                 CORRECT OFFSETS AMONG DIFFERENT SURVEYS
# =============================================================================
# baseline after the eclipse
eclipse = 8000

def baseline(df, eclipse=None):
    if eclipse != None:
        df_baseline = np.nanmedian(df[df['HJD']>eclipse]['mag']) 
    else:
        df_baseline = np.nanmedian(df['mag'])
    return df_baseline

        
        
asas_sn_baseline = baseline(lc_sp_corr, eclipse=eclipse)
atlas_baseline = baseline(df_atlas_corr, eclipse=eclipse)
aavso_baseline = baseline(df_aavso_corr, eclipse=None)
crts_baseline = baseline(df_crts, eclipse=None)
ptf_baseline= baseline(df_ptf, eclipse=None)
ztf_baseline = baseline(df_ztf, eclipse=None)

def blend_asas_sn(df, asas_sn_bl, survey_bl):
    offset = asas_sn_bl - survey_bl
    df_blend = df.copy()
    df_blend['mag'] = df['mag'] + offset
    return df_blend
    
    
# offset = asas_sn_baseline - df_atlas_baseline
# df_atlas_corr_blend = df_atlas_corr.copy()
# df_atlas_corr_blend['mag'] = df_atlas_corr['mag'] + offset


df_atlas_corr_blend = blend_asas_sn(df_atlas_corr, asas_sn_baseline, atlas_baseline)
df_aavso_corr_blend = blend_asas_sn(df_aavso_corr, asas_sn_baseline, aavso_baseline)
df_crts_blend = blend_asas_sn(df_crts, asas_sn_baseline, crts_baseline)
df_ptf_main_blend = blend_asas_sn(df_crts, asas_sn_baseline, ptf_baseline)
df_ztf_blend = blend_asas_sn(df_ztf, asas_sn_baseline, ztf_baseline)

fig, ax = plt.subplots(1, figsize=(16,5))



# lc_sp_corr = band_offset(lc_sp, y='mag')
_ = plot.plot_lc(lc_sp_corr, y='mag', ax=ax, errorbars=True, label='ASAS_SN', invert=True)
_ = plot.plot_lc(df_atlas_corr_blend, y='mag', ax=ax,errorbars=True, label='ATLAS')
_ = plot.plot_lc(df_crts_blend, y='mag', ax=ax,errorbars=True, label='CRTS')
_ = plot.plot_lc(df_aavso_corr_blend, y='mag', ax=ax,errorbars=True, label='AAVSO')
_ = plot.plot_lc(df_ptf_main_blend, y='mag', ax=ax,errorbars=True, label='PTF')
_ = plot.plot_lc(df_ztf_blend, y='mag', ax=ax,errorbars=True, label='ZTF')
# ax.plot(df_tess['HJD'], df_tess['flux'], '.', label='TESS')

ax.set_ylim(17.5,13)
ax.legend(loc=(1.01,0))
plt.show()
# ax.set_xlim(6500, 9000)
# plt.savefig('all_merged_zoomed-in', dpi=200, bbox_inches='tight')
#%% Normalised flux from ASAS_SN by computing the baseline out of the eclipse

lc_sp_corr_flat = lc_sp_corr[lc_sp_corr['HJD']>7900]
mean_flux = np.nanmedian(lc_sp_corr_flat['flux'])
lc_sp_corr['flux'] = lc_sp_corr['flux'] / mean_flux



fig, ax = plt.subplots(1, figsize=(16,5))




_ = plot.plot_lc(lc_sp_corr, str(gaia_id), y='mag', ax=ax, errorbars=True,label='ASAS_SN')
_ = plot.plot_lc(df_atlas_corr, str(gaia_id), y='mag', ax=ax,errorbars=True, label='ATLAS')
_ = plot.plot_lc(df_crts, str(gaia_id), y='mag', ax=ax,errorbars=True, label='CRTS')
_ = plot.plot_lc(df_aavso, str(gaia_id), y='mag', ax=ax,errorbars=True, label='AAVSO')


survey = 'ptf'

survey_files = glob.glob('gdr2_30169/'+survey+'*')
nrows = len(survey_files)
# ax[0].set_title(survey)
for i in range(1, nrows+1):
    file = 'gdr2_30169/' + survey + '_V0658Ori_' + str(i) +'.csv'
    print('Loading ', file)
    df = load_tf(file)
    df = flux_from_mag(df)
    label_str = survey+str(i)
    _ = plot.plot_lc(df, str(gaia_id), y='mag', ax=ax, errorbars=True, title=False, label=label_str)
    
survey = 'ztf'

survey_files = glob.glob('gdr2_30169/'+survey+'*')
nrows = len(survey_files)
# ax[0].set_title(survey)
for i in range(1, nrows+1):
    file = 'gdr2_30169/' + survey + '_V0658Ori_' + str(i) +'.csv'
    print('Loading ', file)
    df = load_tf(file)
    df = flux_from_mag(df)
    label_str = survey+str(i)
    _ = plot.plot_lc(df, str(gaia_id), y='mag', ax=ax, errorbars=True, title=False, label=label_str)
    
    
# ax.plot(df_tess['HJD'], df_tess['flux'], '.', label='TESS')
# ax.invert_yaxis()
ax.set_ylim(17.5,13)
ax.legend(loc=(1.005,0))
plt.show()
# plt.savefig('merged_lightcurve', dpi=200, bbox_inches='tight')

#%%
# =============================================================================
#                               ASAS
# =============================================================================
'''
From: Grzegorz Pojmanski <gp@astrouw.edu.pl>
Subject: Re: New Exciting Star! Weird Dipper with Massive Eclipse

So far I have extracted data from Northern and Southern V archives and
some I data from South. See attached files.
In principle I should be able to get southern data after 2010 - but I
have to reduce this particular field, which will take some time.
'''

'''From ASAS website: 
    for magnitudes fainter then 12 use aperture_number = 0 (MAG_0)
    for magnitudes brighter then 9 use aperture_number = 4 (MAG_4)
    otherwise use:
    aperture_number = 12 - magnitude.
'''

survey = 'asas'
survey_files = glob.glob('gdr2_30169/'+survey+'*')


file = survey_files[0]
# array = np.loadtxt(file)
# df = pd.read_csv(file, sep='\s')

def load_asas(file, ax=None):
    
    hdr = ['HJD', 'MAG_0', 'MAG_1',  'MAG_2',  'MAG_3',  'MAG_4', 'MER_0', 'MER_1',
               'MER_2', 'MER_3', 'MER_4','GRADE','FRAME','FLAG']
    data = ascii.read(file, delimiter='\s', names=hdr)
    df = data.to_pandas()
    df['HJD'] = df['HJD'] + 2450000
    df = df[df['GRADE']=='A']
    
    # fig, ax = plt.subplots(1, figsize=(12,5))
    ax = ax or plt.gca()
    for i in range(1):
        mag = 'MAG_'+str(i)
        errmag = 'MER_'+str(i)
        ax.errorbar(df['HJD'], df[mag], df[errmag], fmt='.', label=mag)
        
    

 
    return df, ax

fig, ax = plt.subplots(3, figsize=(16,12))

for i in range(3):
    _ , _ = load_asas(survey_files[i], ax=ax[i])
    ax[i].invert_yaxis()

ax[2].set_xlabel('HJD')
ax[1].set_ylabel('Magnitude')
# ax.set_title('ASAS')
ax[0].legend(ncol=4, loc=(0.25,1.015))   
plt.show()
# plt.savefig('asas_3-lightcurves', dpi=200, bbox_inches='tight')

#%% Get MAG_0 from ASAS
# pick only file 1 and 3 (file 2 is noisy)
df_asas_0, _ = load_asas(survey_files[0])
temp, _ = load_asas(survey_files[2])
df_asas_0 = df_asas_0.append(temp, ignore_index=True)
df_asas_0 = df_asas_0.rename(columns={'MAG_0':'mag', 'MER_0':'mag_err'})


df_asas_0.to_csv('asas_0.csv', index=False)


