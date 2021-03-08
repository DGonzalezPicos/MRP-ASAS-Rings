"""
Created on Fri Mar  5 14:55:33 2021
@author: dario

Title: 
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import time
from pyasassn.client import SkyPatrolClient
client = SkyPatrolClient("bad_client", "a5a55N_CLIENT")
client.catalogs.stellar_main.head(12)

client.catalogs.master_list

#%%

start = time.time()
gaia_ids = np.loadtxt('zariv2.cat')
query = client.query_list(list(gaia_ids[:20]), catalog='stellar_main',
                          id_col='gaia_id', download=True)
end = time.time()
print('Elapsed time...{:.2f}'.format(end-start))

stats = query.stats()
asas_id = query.catalog_info['asas_sn_id']

lc = query[asas_id[2]]
lc.data['flux_err']
lc.plot()
plt.savefig('images/API-LC.png', dpi=100, bbox_inches='tight')

np.mean(stats['epochs'].values)
#%%
# =============================================================================
#               GET SPECIFIC TARGET 
# =============================================================================
coords = [90.003155, -31.007731]
lc_j060= client.cone_search(ra_deg=coords[0], dec_deg=coords[1], radius=0.01, 
                   catalog='master_list', download=True)
lc_j060 = lc_j060[lc_j060.catalog_info['asas_sn_id'][0]]
#%%
from lightkurve import LightCurve

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
#%%
from astropy.timeseries import LombScargle
frequency, power = LombScargle(t=lk.time, y=lk.flux, dy=lk.flux_err).autopower()
plt.plot(frequency, power)

best_frequency = frequency[np.argmax(power)]
folded_lk = lk.fold(period=1/best_frequency)
folded_lk.scatter()