"""
Created on Sun Feb 14 23:24:46 2021
@author: dario

Title: Plotting lightcurves from 2 ASAS-SN sources 
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os

basepath = '/home/dario/AstronomyLeiden/FRP/ASAS' # path to images
os.chdir(basepath)

# Read data from multiple files
def loadmultdata(files):
    hjd, mag, errmag = [[] for i in range(3)]
    for f in files:
        data = pd.read_csv(f)
        hjd.append(data['hjd'])
        mag.append(data['mag'])
        errmag.append(data['mag err'])
    return hjd, mag, errmag

# Probably better to put the csv files on a folder and get
# all the files from the folder 

# Call function
files = ['AP18326566.csv', 'AP18326567.csv', 'AP28742359.csv']
hjd, mag, errmag = loadmultdata(files)

# Plot LCs
fig, ax = plt.subplots(len(hjd), 1, figsize=(15,16))
   
# Make common labels
fs = 15 # fontsize
fig.text(0.5, 0.04, 'HJD', ha='center', fontsize=fs)
fig.text(0.07, 0.5, 'Magnitude (V)', va='center', rotation='vertical', fontsize=fs)
# Plot all the datasets
for j in range(len(hjd)):
    if j == 0:
        ax[j].get_shared_x_axes().join(ax[j], ax[j-1])
        ax[j].set_xticklabels([])
        
    ax[j].plot(hjd[j], mag[j], '.k')
    ax[j].errorbar(hjd[j], mag[j], errmag[j],fmt='k.', ecolor='.7', markersize=4)
    ax[j].text(0.87, 0.05, files[j][:-4],transform=ax[j].transAxes,fontsize=15)
     
plt.savefig('/home/dario/AstronomyLeiden/FRP/MRP-ASAS-Rings/LC1.png', dpi=100, bbox_inches='tight')
plt.show()
#%% 
df = pd.read_csv(files[0])
df.columns
files[0][:-4]
