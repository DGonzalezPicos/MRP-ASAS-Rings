"""
Created on Sun Feb 14 23:24:46 2021
@author: dario

Title: Plotting lightcurves from 2 ASAS-SN sources 
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os, glob

basepath = '/home/dario/AstronomyLeiden/FRP/MRP-ASAS-Rings/data' # path to images
os.chdir(basepath)


    
# Read data from multiple files
def loadmultdata(files):
    hjd, mag, errmag = [[] for i in range(3)]
    for f in files:
        data = pd.read_csv(f)
        hjd.append(data['HJD'])
        mag.append(data['mag'])
        errmag.append(data['mag_err'])
    return hjd, mag, errmag


# Get all csv files from the folder
files = glob.glob('*csv')
# pd = pd.read_csv(files[0])
# pd.columns
hjd, mag_str, errmag = loadmultdata(files)
len(mag_str[0])
len(mag_str[1])

type(mag_str[0][0])
type(mag_str[1][0])

# Fix the data format for magnitude data
mag = []
for j in range(len(mag_str)):
    if type(mag_str[j][0]) == str:
        temp = []
        for k in mag_str[j]:
            temp.append(float(k.split('>')[-1]))
        mag.append(temp)
    else:
        mag.append(mag_str[j])

# Fix "infinite values" for the "errmag" data by
# replacing them with the median
for j in range(len(errmag)):
    for count, k in enumerate(errmag[j]):
        if k > 90: # the wrong values have an assigned "99.99"
            errmag[j][count] = np.median(errmag[j])

#%%
# Plot LCs
fig, ax = plt.subplots(len(hjd), 1, figsize=(9,12))
   
# Make common labels
fs = 15 # fontsize
fig.text(0.5, 0.04, 'HJD', ha='center', fontsize=fs)
fig.text(0.02, 0.5, 'Magnitude (V)', va='center', rotation='vertical', fontsize=fs)
# Plot all the datasets
for j in range(len(hjd)):
    # if j == 0:
    #     ax[j].get_shared_x_axes().join(ax[j], ax[j-1])
    #     ax[j].set_xticklabels([])
    if len(hjd) > 1:  
        med = np.median(mag[j])
        sigma = np.sqrt(np.std(mag[j]))/(3.5)
        ax[j].set_ylim(med-(0.6*sigma), med+sigma)
        ax[j].invert_yaxis()
        
        
        ax[j].plot(hjd[j], mag[j], '.k')
        ax[j].errorbar(hjd[j], mag[j], errmag[j],fmt='k.', ecolor='.7', markersize=2)
        ax[j].text(0.87, 0.05, files[j][:-4],transform=ax[j].transAxes,fontsize=15)
    # if len(hjd) == 1:
    #     ax.plot(hjd[j], mag[j], '.k')
    #     ax.errorbar(hjd[j], mag[j], errmag[j],fmt='k.', ecolor='.7', markersize=4)
    #     ax.text(0.87, 0.05, files[j][:-4],transform=ax.transAxes,fontsize=15)
     
plt.savefig('/home/dario/AstronomyLeiden/FRP/MRP-ASAS-Rings/LC1.png', dpi=100, bbox_inches='tight')
plt.show()
#%% 
df = pd.read_csv(files[0])
df.columns
files[0][:-4]
