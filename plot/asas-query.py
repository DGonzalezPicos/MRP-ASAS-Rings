"""
Created on Wed Mar 10 19:14:34 2021
@author: dario

Title: Reading LCs from ASAS-SN query
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
import pandas as pd
import os
import lightkurve as lightk

data = ascii.read("stats.txt") 

data['mean_mag'].info.format = '2.3f'
data['std_mag'].info.format = '2.3f'


np.max(data['mean_mag'])


# Plot statistics from the dataset (two histograms)
'''
fig, (ax1,ax2) = plt.subplots(1,2, figsize=(12,6))
#######################
ax1.hist(data['mean_mag'], bins='auto', color='darkcyan', alpha=0.8)
ax1.set_xlim(8,22)
ax1.set_ylabel('Number of targets')
ax1.set_xlabel('Magnitude')
#######################
ax2.hist(data['epochs'], bins='auto', color='seagreen', alpha=0.8)
# ax1.set_xlim(8,22)
ax2.set_ylabel('Number of targets')
ax2.set_xlabel('Epochs')

plt.show()
# plt.savefig('./images/lc-statistics.png', dpi=150, bbox_inches='tight')
'''
#%%
# =============================================================================
#                           READ LCs
# =============================================================================
size = 10

cwd = os.getcwd()
os.chdir(os.path.join(cwd,'zv2'))
file_list  = os.listdir()
rand = np.random.randint(0,2e3)
df = pd.read_csv(file_list[rand], sep=';')

# df.columns
# df['asas_sn_id'][0]
def overview(lc, name='lc', save=False):
    lk = lightk.LightCurve(time=lc['jd'].values, flux=lc['flux'].values,
                           flux_err=lc['flux_err'].values)
    
    fig, ax1 = plt.subplots(1,2, figsize=(12,4))
    ax1.set_title('Raw Lightcurve')
    ax2.set_title('Flattened LightCurve')
    lk.errorbar(ax=ax1, color='k', label=lc['asas_sn_id'][0])
    flat_lk.errorbar(ax=ax2, color='k', label=lc['asas_sn_id'][0])
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
        
quick_view = overview(df, name='J060000')




