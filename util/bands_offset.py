"""
Created on Mon May  3 09:28:05 2021
@author: dario

Title: bands-offset
Compute offset between V and g band data with same HJD 
"""
import numpy as np
import matplotlib.pyplot as plt
import sv_util as sv
import glob, os
import pandas as pd
import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 2
# mpl.rcParams['axes.labelsize'] = 14
# mpl.rcParams['legend.fontsize'] = 14
mpl.rcParams['font.size'] = 14
mpl.rcParams.keys()


# Define Class Objects
vband = sv.LightCurveSet('/home/dario/AstronomyLeiden/FRP/leiden_vband/lc/camfix', 'dat', 'v')

gband = sv.LightCurveSet('/home/dario/AstronomyLeiden/FRP/gband_full', 'dat', 'g')

vband.files(gaia_id=True)
gband.files(gaia_id=True)

# Find files with data on both bands (crossmatch Gaia_ID)
intersect = set(gband.files(gaia_id=True)).intersection(vband.files(gaia_id=True))

#%% Merge Bands into new files saved at 'outpath' with filename GAIA_ID
print('-------------------------------------')
print('Merging bands...')
for name in list(intersect):
    file_v = os.path.join(vband.path, str(name) + '.dat')
    file_g = os.path.join(gband.path, str(name) + '.dat')
    # print(file_v, file_g)
    
    lc_v = vband.data(file_v)
    lc_g = gband.data(file_g)
    
    # lc_v['HJD'] = lc_v['HJD'] - 2450000
    # lc_g['HJD'] = lc_g['HJD'] - 2450000
    
    # sort values by date
    lc_g = lc_g.sort_values(by=['HJD'])
    
    try:
        print('End V-band: {:.0f} and start g-band {:.0f}'.format(lc_v['HJD'].iloc[-1], lc_g['HJD'].iloc[0]))

        # create mask for overlapping region to compute the median
        v_overlap_mask = lc_v['HJD'] > np.amin(lc_g['HJD'].iloc[0])
        g_overlap_mask = lc_g['HJD'] < np.amax(lc_v['HJD'].iloc[-1])
        
        flux_v_med = np.nanmedian(lc_v['flux'][v_overlap_mask])
        flux_g_med = np.nanmedian(lc_g['flux'][g_overlap_mask])
        
        med_diff = 0.5 * (flux_v_med - flux_g_med)
        
        # fix the flux values by adding or subtracting the median
        lc_v['flux'] -= med_diff
        lc_g['flux'] += med_diff
         # convert to pandas dataframe to manipulate
        df_merge = pd.concat([lc_v, lc_g])
        asas_sn_id = int(df_merge['asas_id'].iloc[-1])
        df_merge['asas_id'] = asas_sn_id
        
        array = df_merge.to_numpy()
        outpath = os.path.join('/home/dario/AstronomyLeiden/FRP/merged-bands-offset/',str(name)+'.dat')
        hdr = 'HJD \t flux \t errflux \t gaia_id \t asas_id'
        np.savetxt(outpath, array, header=hdr)
    
    except:
        print('Empty file...')
        continue
print('Finished!!!')
print('-------------------------------------')

output_filelist = glob.glob('/home/dario/AstronomyLeiden/FRP/merged-bands-offset/*.dat')
len(output_filelist )


#%% Compare the offset fix for an individual file

name = list(intersect)[319]

file_v = os.path.join(vband.path, str(name) + '.dat')
file_g = os.path.join(gband.path, str(name) + '.dat')
# print(file_v, file_g)

lc_v = vband.data(file_v)
lc_g = gband.data(file_g)

lc_v['HJD'] = lc_v['HJD'] - 2450000
lc_g['HJD'] = lc_g['HJD'] - 2450000

# sort values by date
lc_g = lc_g.sort_values(by=['HJD'])
print('End V-band: {:.0f} and start g-band {:.0f}'.format(lc_v['HJD'].iloc[-1], lc_g['HJD'].iloc[0]))
overlap = lc_g['HJD'].iloc[0] - lc_v['HJD'].iloc[-1]
if overlap < 0:
    print('Overlap')
    
v_overlap_mask = lc_v['HJD'] > np.amin(lc_g['HJD'].iloc[0])
g_overlap_mask = lc_g['HJD'] < np.amax(lc_v['HJD'].iloc[-1])

flux_v_med = np.nanmedian(lc_v['flux'][v_overlap_mask])
flux_g_med = np.nanmedian(lc_g['flux'][g_overlap_mask])

med_diff = 0.5 * (flux_v_med - flux_g_med)

new_flux_v = lc_v['flux'] - med_diff
new_flux_g = lc_g['flux'] + med_diff

new_hjd = np.append(lc_v['HJD'], lc_g['HJD'])
new_flux = np.append(new_flux_v, new_flux_g)
# Date overlap

fig, (ax1, ax2) = plt.subplots(2,1, figsize=(16,9), sharex=True, sharey=True)

ax1.scatter(lc_v['HJD'], lc_v['flux'], label='V-band', color='orange', alpha=0.6, s=8)
ax1.scatter(lc_g['HJD'], lc_g['flux'], label='g-band', color='green', alpha=0.6, s=8)

ax2.scatter(lc_v['HJD'], new_flux_v, label='V-band', color='orange', alpha=0.6, s=8)
ax2.scatter(lc_g['HJD'], new_flux_g, label='g-band', color='green', alpha=0.6, s=8)
# ax2.scatter(new_hjd, new_flux)
# ax2.set_xlabel('HJD - 2450000')

ax1.legend()
ax2.legend()

ymin_auto, ymax_auto = ax1.yaxis.get_data_interval()
ax1.set_ylim(ymin=np.max([ymin_auto, 0.0]))
ax1.set_ylim(ymax=np.min([ymax_auto, 2.0]))

plt.show()






