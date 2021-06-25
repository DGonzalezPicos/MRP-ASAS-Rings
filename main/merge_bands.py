"""
Created on Tue Apr 20 00:36:56 2021
@author: dario

Title: Manage and merge LCs
"""
import numpy as np
import matplotlib.pyplot as plt
import support_functions as supp
import plot_functions as plot
import glob, os
import pandas as pd
import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 2
# mpl.rcParams['axes.labelsize'] = 14
# mpl.rcParams['legend.fontsize'] = 14
mpl.rcParams['font.size'] = 14
mpl.rcParams.keys()


# Define Class Objects
print('-------Generating LC objects-------')
vband = supp.LightCurveSet('/home/dario/AstronomyLeiden/FRP/leiden_vband/camfix/', 'dat', 'v')
gband = supp.LightCurveSet('/home/dario/AstronomyLeiden/FRP/gband_api_fixed/', 'dat', 'g')

vband.path

vband.files(gaia_id=True)
gband.files(gaia_id=True)

# Find files with data on both bands (crossmatch Gaia_ID)
intersect = list(set(gband.files(gaia_id=True)).intersection(vband.files(gaia_id=True)))
print('Intersection completed: {:} objects in both bands'.format(len(intersect)))

#%% Merge Bands into new files saved at 'outpath' with filename GAIA_ID

'''
print('-------------------------------------')
print('Merging bands...')
for name in intersect:
    file_v = os.path.join(vband.path, str(name) + '.dat')
    file_g = os.path.join(gband.path, str(name) + '.dat')
    # print(file_v, file_g)
    
    lc_v = vband.data(file_v)
    ##### added to solve new issue with non-numerical values '>17.33'
    lc_v['mag_err'] = lc_v['mag_err'].replace(99.99, np.nan)
    lc_v = lc_v.dropna()
    ###########
    lc_g = gband.data(file_g)
    
    
    # compute the median magnitude
    lc_v['mag'] = pd.to_numeric(lc_v['mag'])
    mean_mag = np.nanmedian(lc_v['mag'])
    
    
    # convert to pandas dataframe to manipulate
    df_merge = pd.concat([lc_v, lc_g])
    asas_sn_id = int(df_merge['asas_id'].iloc[-1])
    df_merge['asas_id'] = asas_sn_id
    
    
    df_merge = df_merge.drop(columns=['camera', 'mag', 'mag_err', 'FWHM'])
    df_merge['mean_mag'] = mean_mag
    
    df_merge.dtypes
    array = df_merge.to_numpy()
    outpath = os.path.join('/home/dario/AstronomyLeiden/FRP/merged-bands/', str(name)+'.dat')
    hdr = 'HJD \t flux \t flux_err \t gaia_id \t asas_id \t mean_mag'
    
    # type_list = ['%16.8e', '%s', '%16.8e', '%16.8e', '%16.8e', '%16.8e', '%16.8e', '%16.8e', '%16.8e']
    print('Saving {:}'.format(outpath))
    # np.savetxt(outpath, array, header=hdr)
print('Finished!!!')
print('-------------------------------------')
# float(lc_g['HJD'][0])
'''


#%% PLot an example 
# gaia_id = 239547059291163008

gaia_id = intersect[2228]


file_v = os.path.join(vband.path, str(gaia_id) + '.dat')
file_g = os.path.join(gband.path, str(gaia_id) + '.dat')
# print(file_v, file_g)

lc_v = vband.data(file_v)
##### added to solve new issue with non-numerical values '>17.33'
lc_v['mag_err'] = lc_v['mag_err'].replace(99.99, np.nan)
lc_v = lc_v.dropna()
###########
lc_g = gband.data(file_g)


fig, (ax, ax1) = plt.subplots(2,1, figsize=(12,6))


ax.errorbar(lc_v['HJD']- 2450000, lc_v['flux'], lc_v['flux_err'], 
                fmt='.', label='V-band', color='orange', alpha=0.8)

ax.errorbar(lc_g['HJD']- 2450000, lc_g['flux'], lc_g['flux_err'],
                fmt='.', label='g-band', color='green', alpha=0.8)
########
ax1.errorbar(lc_v['HJD']- 2450000, lc_v['mag'], lc_v['mag_err'], 
                fmt='.', label='V-band', color='orange', alpha=0.8)

ax1.errorbar(lc_g['HJD']- 2450000, lc_g['mag'], lc_g['mag_err'],
                fmt='.', label='g-band', color='green', alpha=0.8)
############

ax.set_title('GAIA_ID: {:}'.format(str(gaia_id)))

ax.set_ylabel('Relative flux')
ax1.set_ylabel('Magnitude')
ax1.invert_yaxis()

ax1.set_xlabel('HJD - 2450000')
ax.legend()
# plt.savefig('testing_band_merging.png', bbox_inches='tight', dpi=200)
plt.show()

