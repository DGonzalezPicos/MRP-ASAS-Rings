"""
Created on Fri Mar  5 14:55:33 2021
@author: dario

Title: Query and download LCs for a given catalog or object coordinates
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time, glob
from lightkurve import LightCurve
from pyasassn.client import SkyPatrolClient
import warnings
from astropy.utils.exceptions import AstropyWarning
warnings.simplefilter('ignore', category=AstropyWarning)
from astropy.io import fits

print('------------------------------')
print('Loading SkyPatrolClient...')

client = SkyPatrolClient("bad_client", "a5a55N_CLIENT")
client.catalogs.stellar_main.head(12)
client.catalogs.catalog_names()
cat_list = client.catalogs
lcs = client.query_list([3207233877299221504], catalog='stellar_main', id_col='gaia_id', download=True)
lcs.data
# lcs.data.loc[10]
query_gaia_id = 5856351499141796736
lcs_info = client.query_list(query_gaia_id, catalog='stellar_main', id_col='gaia_id', download=False)
lcs_info
# lcs = client.query_list([5994772213404409216], catalog='stellar_main', id_col='gaia_id', download=True)
# lcs.data.columns
# plt.plot(lcs.data['jd'], lcs.data['flux'],'.k')
# lcs.data.tail()
# int(lcs['asas_sn_id'])
# int(lcs['gaia_id'])



#%%

# start_tot = time.time()
# gaia_ids = np.loadtxt('zariv2.cat') # OLD 


# Load Zari's catalog
# filename = '/home/dario/AstronomyLeiden/FRP/pms.fits'
# with fits.open(filename) as hdul:
#     print('Reading Zari\'s catalog...')
#     # hdul.info()
#     data = hdul[1].data
    


# Read file list from Tharindu's V-band lightcurves
v_path = '/home/dario/AstronomyLeiden/FRP/leiden_vband/lc/camfix/'
filelist = glob.glob(v_path+'*.dat')
gaia_ids = []

print('Loading gaia_ids from {:}'.format(v_path))
for file in filelist:
    name = file.split('/')[-1][:-4]
    gaia_ids.append(int(name))
    
lcs = client.query_list(gaia_ids, catalog='stellar_main', id_col='gaia_id', download=False)
# lcs
#%%

# f = open('rosetta.out', 'w')
start = time.time()
print('Starting query...')
output = client.query_list(gaia_ids[0], catalog='stellar_main', id_col='gaia_id', download=False)

for gaia_id in gaia_ids[1:10]:
    print('Querying {:}'.format(gaia_id))
    output = output.append(client.query_list(gaia_id, catalog='stellar_main', id_col='gaia_id', download=False), ignore_index=True)
    
print('Finished!!')   

elapsed = time.time() - start
print('Elapsed time {:.2f} s --- {:.2f} min'.format(elapsed, elapsed/60)) 

print(output.tail())

#%%
# Galactic coordinates
# gaia_ids = data['source_id'][3:6]

print('------------------------------')
print('------------------------------')
print('Downloading lightcurves...')
print('------------------------------')
missing = []
start = time.time()
outpath = '/home/dario/AstronomyLeiden/FRP/gband_api/'

step = 400
recall = 0 # MANUAL
for i in range(int(np.round(len(gaia_ids)/step))):
# for i in range(2):    
    print('-------------- {:} --------------'.format(i+recall))
    
    ind = step*(i+recall)
    ind2 = step*(i+1+recall)
    
    if ind > len(gaia_ids):
        # save missing LCs id's
        hdr_missing = 'Number of missing LCs: {:}'.format(len(missing))
        np.savetxt('missing.out', missing, header=hdr_missing)
        print('Done!!!!!!!!!!!!!!!!!')
        end = time.time()
        print('Elapsed time...{:.2f} s'.format(end-start))
        print('------------------------------')
        print('------------------------------')
        break
    if ind2 > len(gaia_ids):
        print('Last query...')
        lcs = client.query_list(gaia_ids[ind:], catalog='stellar_main', id_col='gaia_id', download=True)
    else:
        print('Querying...')
        lcs = client.query_list(gaia_ids[ind:ind2], catalog='stellar_main', id_col='gaia_id', download=True)
      
        
    
    # lcs.data = lcs.data.drop(columns=['cam', 'fwhm'])    # add cam info    
    ids = list(set(lcs.data['asas_sn_id']))
    print('Retrieval ratio of {:} %'.format(100*len(ids)/400))
    for asas_id in ids:
        df = lcs.data.loc[lcs.data['asas_sn_id'] == asas_id]
        array = df.to_numpy()
        hdr = 'HJD flux flux_err mag mag_err limit asas_sn_id'
        filename = outpath + str(asas_id)+'.dat'
        fmt = ['%10.6f', '%6.6f', '%6.6f', '%10.6f', 
               '%10.6f', '%10.6f', '%1.2f', '%s', '%d', '%d']
        
        print('Writing {:} to .dat file...'.format(asas_id))
        # np.savetxt(filename, array, header=hdr, fmt=fmt)
        
        

# save missing LCs id's
hdr_missing = 'Number of missing LCs: {:}'.format(len(missing))
np.savetxt('missing.out', missing, header=hdr_missing)
print('Done!!!!!!!!!!!!!!!!!')
end = time.time()
print('Elapsed time...{:.2f} s'.format(end-start))
print('------------------------------')
print('------------------------------')

# lcs.data.head()
# lcs.data['asas_sn_id']
# lcs.data.loc[lcs.data['asas_sn_id'] == 609885877102]

# lcs.data.columns

# plt.errorbar(lcs.data['jd'], lcs.data['flux'], lcs.data['flux_err'],'.')


# lcs.data.loc[15]
# # lcs.data.plot()
# plt.errorbar(df['hjd'], df['flux'], df['errflux'])
# plt.show()
#%%
# FIX FILE NAMES (from ASAS_ID TO GAIA_ID)

file_list = glob.glob(outpath+'/*.dat')
ids = []


for file in file_list:
    name = file.split('/')[-1][:-4]
    if int(name) < 1e12: 
        ids.append(int(name)) # get all the ASAS_ID
# help(client.query_list)  
     
start = time.time()
out = client.query_list(gaia_ids[:10], catalog='stellar_main', id_col='gaia_id', download=False)

end = time.time()
print('Time {:.1f} minutes'.format((end-start)/60))


#%%
outpath = '/home/dario/AstronomyLeiden/FRP/gband_full/'
inter = set(ids).intersection(out['asas_sn_id'])

for file in list(inter):
    filename = os.path.join(outpath, str(file)+'.dat')
    data = np.loadtxt(filename, usecols=[0,1,2,-1])
    gaia_id = str(out.loc[out['asas_sn_id']==file]['gaia_id'].values[0])
    outname = os.path.join(outpath, gaia_id+'.dat')
    np.savetxt(outname, X=data)
    print('Saving file with GAIA_ID name!')
    

#%% RENAME THE FILES
for i in range(len(out['gaia_id'][:20])):
    filename = os.path.join(outpath, str(out['asas_sn_id'][i])+'.dat')
    # file2 = os.path.join(outpath, str(ids[0])+'.dat')
    try:
        data = np.loadtxt(filename, usecols=[0,1,2,-1])
        
        outname = os.path.join(outpath, out['gaia_id'][i]+'.dat')
        np.savetxt(outname, X=data)
        print('Saving file with GAIA_ID name!')
    except:
        print('No file with this ASAS_SN_ID....', str(out['asas_sn_id'][i]))
    
    
out['gaia_id'][0]
out['asas_sn_id'][0]

out['gaia_id'].loc([out['asas_sn_id']==369368033926])
str(out.loc[out['asas_sn_id']==369368033926]['gaia_id'].values[0])


#%%
# =============================================================================
#               GET SPECIFIC TARGET -----> J0600
# =============================================================================
'''
coords = [90.003155, -31.007731]
lc_j060= client.cone_search(ra_deg=coords[0], dec_deg=coords[1], radius=0.01, 
                   catalog='master_list', download=True)
lc_j060 = lc_j060[lc_j060.catalog_info['asas_sn_id'][0]]
#%%

def overview(lc, name='lc', save=False):
    lk = LightCurve(time=lc.data['jd'], flux=lc.data['flux'], flux_err=lc.data['flux_err'])
    # lk.scatter()
    flat_lk = lk.flatten(window_length=101, polyorder=2, niters=3)
    # flat_lk.scatter()
    
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12,4))
    ax1.set_title('Raw Lightcurve')
    ax2.set_title('Flattened LightCurve')
    lk.errorbar(ax=ax1, color='k')
    flat_lk.errorbar(ax=ax2, color='k')
    plt.show() 
    if save == True:
        plt.savefig('images/'+name+'_raw-flat.png', dpi=200, bbox_inches='tight')

    pg = lk.to_periodogram(oversample_factor=1)
    # pg.plot()
    # pg.plot(view='period', scale='log')
    pg.period_at_max_power
    
    # lk.fold(period=pg.period_at_max_power).scatter()
    
    fig, (ax1, ax2) = plt.subplots(1,2, figsize=(12,4))
    pg.plot(ax=ax1, color='k')
    pg.plot(ax=ax2, color='k',view='period', scale='log')
    ax1.set_title('Periodogram')
    ax2.set_title('Periodogram -- period space')
    # lk.fold(period=pg.period_at_max_power).scatter()
    plt.show() 
    if save == True:
        plt.savefig('images/'+name+'_periodogram.png', dpi=200, bbox_inches='tight')
        
quick_view = overview(lc_j060, name='J060000')
'''