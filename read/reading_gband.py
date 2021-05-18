"""
Created on Sat May  8 18:05:37 2021
@author: dario

Title: Read g-band LCs
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob, os
from astropy.io import ascii 
from pyasassn.client import SkyPatrolClient

path = '/home/dario/AstronomyLeiden/FRP/gband_full/'
filelist = glob.glob(path+'*.dat')


# file = filelist[17]
# gaia_id = 3215778701851108352
# gaia_id = 918534285882642176
# gaia_id = 2270245465569124608
# gaia_id =  450701124878486784
# gaia_id = 3228644190486688512
# gaia_id = 138715731286579200
# gaia_id = 5513727389379855360
gaia_id = 5512322694556142976 # andy 1
# file = os.path.join(path, str(gaia_id)+'.dat')

# data = np.loadtxt(file)
# name = file.split('/')[-1][:-4]

# # df_g = pd.DataFrame(data, columns=['hjd','flux','errflux', 'mag', 'mag_err','cam','asas_sn_id'])
# df_g = pd.DataFrame(data, columns=['hjd','flux','errflux', 'asas_sn_id'])
# df_g = df_g.dropna()

# print(df_g.tail())


#%% V-band data
path = '/home/dario/AstronomyLeiden/FRP/leiden_vband/lc/camfix/'
file_v = os.path.join(path, str(gaia_id)+'.dat')
data = ascii.read(file_v, delimiter='\t')
df = data.to_pandas()
df['flux_err'] = df['flux_err'].replace(99.99, np.nan)
df['mag_err'] = df['mag_err'].replace(99.99, np.nan)
df = df.dropna()
df['flux_norm'] = df['flux'] / np.nanmedian(df['flux'])

# Froms str to float
df['mag'] = pd.to_numeric(df['mag'], errors='coerce')

temp = 10**(-df['mag'] / 2.5)
df['flux_from_mag'] = temp / np.nanmedian(temp)


#%% API data
api = True
if api == True:
    print('Loading SkyPatrolClient...')
    client = SkyPatrolClient("bad_client", "a5a55N_CLIENT")
    client.catalogs.stellar_main.head(12)
    
    lc_api = client.query_list(gaia_id, catalog='stellar_main', id_col='gaia_id', download=True)
    lc_api_info = client.query_list(gaia_id, catalog='stellar_main', id_col='gaia_id', download=False)
    lc_api_info
    lc_api.data.keys()

    flux_norm_api = lc_api.data['flux'] / np.nanmedian(lc_api.data['flux'])
    
    lc_api.data['flux_norm'] = lc_api.data['flux'] / np.nanmedian(lc_api.data['flux'])
    
    print('CAMERAS...', lc_api.data['cam'].unique())
    
    lc_api.data[lc_api.data['cam']=='bD']
    lc_api.data[lc_api.data['cam']=='bt']


    
#%%

# name_web = '2270245465569124608_web.csv'
# df_web = pd.read_csv('/home/dario/AstronomyLeiden/FRP/MRP-ASAS-Rings/' + name_web)
# df_web = df_web.rename(columns={'HJD':'jd', 'Camera':'cam'})

# df_web_clean = df_web.replace(99.99, np.nan)
# df_web_clean = df_web_clean.dropna()
# df_web_clean['mag'] = pd.to_numeric(df_web_clean['mag'])
# df_web_clean['flux_norm'] = df_web['flux(mJy)'] / np.nanmedian(df_web['flux(mJy)'])
#%%
'''
fig, ax = plt.subplots(nrows=2, ncols=2,figsize=(12,9)) 

def plot_gcam(df, cam, y, ax=None, **kwargs):
    lc = df[df['cam'] == cam]
    
    ax = ax or plt.gca()
    ax.plot(lc['jd']-2450000, lc[y], '.', alpha=0.6, **kwargs)
    return ax

for cam in lc_api.data['cam'].unique():
    _ = plot_gcam(lc_api.data, cam=cam, y='flux_norm', ax=ax[0][0], label = 'API Cam: {:}'.format(cam))
    _ = plot_gcam(df_web_clean, cam=cam, y='flux_norm', ax=ax[0][0], label = 'Database Cam: {:}'.format(cam), marker='v', markersize=4)
    _ = plot_gcam(lc_api.data, cam=cam, y='mag', ax=ax[0][1], label = 'API Cam: {:}'.format(cam))
    _ = plot_gcam(df_web_clean, cam=cam, y='mag', ax=ax[0][1], label = 'Database Cam: {:}'.format(cam), marker='v', markersize=4)
    
    _ = plot_gcam(lc_api.data, cam=cam, y='flux_norm', ax=ax[1][0], label = 'API Cam: {:}'.format(cam))
    _ = plot_gcam(df_web_clean, cam=cam, y='flux_norm', ax=ax[1][0], label = 'Database Cam: {:}'.format(cam), marker='v', markersize=4)
    _ = plot_gcam(lc_api.data, cam=cam, y='mag', ax=ax[1][1], label = 'API Cam: {:}'.format(cam))
    _ = plot_gcam(df_web_clean, cam=cam, y='mag', ax=ax[1][1], label = 'Database Cam: {:}'.format(cam), marker='v', markersize=4)


ax[0][0].plot(df['HJD']-2450000, df['flux_norm'],'.', color='gold')
ax[1][0].plot(df['HJD']-2450000, df['flux_norm'],'.', color='gold')
ax[0][0].set_ylim(0,2)
ax[1][0].set_ylim(0,2)
ax[0][0].legend(ncol=4, loc=(0.0, 1.05))
ax[0][0].text(s=('Gaia_ID: {:}'.format(gaia_id)), x=0.8, y=1.20, transform=ax[0][0].transAxes)
ax[0][1].text(s=('Mean V-mag: {:.2f}'.format(np.nanmedian(df['mag']))), x=0.02, y=0.9, transform=ax[0][1].transAxes)
ax[0][1].text(s=('Mean g-mag: {:.2f}'.format(np.nanmedian(lc_api.data['mag']))), x=0.02, y=0.8, transform=ax[0][1].transAxes)

ax[0][0].text(s='V+g flux', x=0.05, y=0.1, transform=ax[0][0].transAxes)
ax[1][0].text(s='g-flux', x=0.05, y=0.1, transform=ax[1][0].transAxes)
ax[0][1].text(s='V+g mag', x=0.05, y=0.1, transform=ax[0][1].transAxes)
ax[1][1].text(s='g-mag', x=0.05, y=0.1, transform=ax[1][1].transAxes)

ax[1][0].set_xlim(df_web['jd'].iloc[0]-2450000,df_web['jd'].iloc[-1]-2450000)
ax[1][1].set_xlim(df_web['jd'].iloc[0]-2450000,df_web['jd'].iloc[-1]-2450000)
ax[1][0].set_xlabel('HJD-2450000')
ax[1][1].set_xlabel('HJD-2450000')

ax[0][1].invert_yaxis()
ax[1][1].invert_yaxis()

ax[0][1].plot(df['HJD']-2450000, df['mag'],'.', color='gold')
ax[1][1].plot(df['HJD']-2450000, df['mag'],'.', color='gold')

plt.show()
# plt.savefig('gband_comp.png', dpi=200)
'''




#%%

def cam_corr(df): 
    cameras = df['cam'].unique()
    print('Cameras...', cameras)
    
    
    cam_medflux = []
    cam_info = []
    print('cam \t cam_points \t cam_flux_median')
    print('-----------------------------------------')
    for cam in cameras:
        cam_medflux.append(np.nanmedian(df[df['cam'] == cam]['flux_norm']))
        print('{:} \t\t {:} \t\t\t{:.4f}'.format(cam, len(df[df['cam'] == cam]), np.nanmedian(df[df['cam'] == cam]['flux_norm'])))
        cam_info.append([cam, len(df[df['cam'] == cam]), np.nanmedian(df[df['cam'] == cam]['flux_norm'])])
        
    cam_df = pd.DataFrame(cam_info, columns=['cam', 'cam_points', 'cam_flux_median'])
    cam_df['cam'].astype(str)
    
    main_cam_obj = cam_df[cam_df['cam_points']==np.max(cam_df['cam_points'])]['cam'] # camera with more points
    main_cam = main_cam_obj.values[0]
    cam_offset_obj = cam_df[cam_df['cam']==main_cam]['cam_flux_median']
    main_flux = cam_offset_obj.values[0]
    
    cameras = list(cameras) 
    cameras.remove(main_cam)
    df_list = [pd.DataFrame(df[df['cam'] == main_cam])]
    for cam in cameras:
            offset = np.nanmedian(df[df['cam'] == cam]['flux_norm']) - main_flux
            
            df_temp = pd.DataFrame(df[df['cam'] == cam])
            df_temp['flux_norm'] =  df[df['cam'] == cam]['flux_norm'] - offset
            df_list.append(df_temp)
    df_corr = pd.concat(df_list)  
    return df_corr

df_corr = cam_corr(lc_api.data)

# add GAIA_ID 
df_corr['gaia_id'] = np.full_like(df_corr['flux'], gaia_id)
df_corr.keys()
#%%

fig, ax = plt.subplots(1,1, figsize=(12,6))


ax.plot(df['HJD']-2450000, df['flux_from_mag'], '.', color='gold', label='V-band: flux from mag')
# ax.plot(df['HJD']-2450000, df['flux_norm'],'.', color='gold', label='V-band')
ax.plot(df_corr['jd']-2450000, df_corr['flux_norm'], '.', color='green', 
        alpha=0.6, label='g-band: Cam-fixed')
ax.plot(lc_api.data['jd']-2450000, lc_api.data['flux_norm'], 'x', color='red', 
        alpha=0.6, label='g-band: API raw')

ax.set_ylabel('Normalised flux')
ax.set_xlabel('HJD-2450000')
ymin_auto, ymax_auto = ax.yaxis.get_data_interval()
ax.set_ylim(ymin=np.max([ymin_auto, 0.0]))
ax.set_ylim(ymax=np.min([ymax_auto, 1.3]))
ax.set_title('GAIA_ID = ' + str(gaia_id))


def mean_mag_label(df,x,y, ax=None):
    
    mean_mag = np.median(df['mag'].to_numpy(dtype=float))
    mean_mag_str = 'Mean V-mag = {:.1f}'.format(mean_mag)
    
    if mean_mag > 17:
        color_mag = 'red'
    elif (mean_mag < 17) and (mean_mag > 15):
        color_mag = 'darkorange'
    else:
        color_mag = 'darkgreen'
    ax.text(x=x, y=y, s=mean_mag_str, transform=ax.transAxes, color=color_mag)
    return None
_ = mean_mag_label(df, x=0.02,y=0.05, ax=ax)
ax.legend(loc='upper left')
plt.show()
plt.savefig('cam_corrected.png', dpi=200, bbox_inches='tight')

