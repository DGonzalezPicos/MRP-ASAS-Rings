"""
Created on Sun Feb 14 23:24:46 2021
@author: dario

Title: Plotting lightcurves from 2 ASAS-SN sources 
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
import glob
from astropy.io import ascii

basepath = '/home/dario/AstronomyLeiden/FRP/ASAS' # path to images
os.chdir(basepath)

file_list = glob.glob(basepath + '/*.dat')
file_list[4]
N = len(file_list)
print('Number of files in the directory ... ', N)

#%%
print(pd.read_csv('AP18326566.csv').columns)

# Read data from multiple files

# float(mag[1][0].replace('>', ''))

def loadmultdata(rand_files):
    # converters = {'col3': [ascii.convert_numpy(np.float)]}
    hjd, mag, errmag, names, temp_list = [[] for i in range(5)]
    for i in rand_files:
        data = ascii.read(file_list[i])
        hjd.append(data['col1'])
        
        # for string in data['col3']:
        #     temp_list.append(float(str(string).replace('>', '')))
                
        # mag.append(temp_list)
        mag.append(data['col3'])
        errmag.append(data['col4'])
        
        name = file_list[i].split('/')[-1][2:-4]
        names.append(name)
    return hjd, mag, errmag, names

# Probably better to put the csv files on a folder and get
# all the files from the folder 

# Call function
# files = ['AP18326566.csv', 'AP18326567.csv', 'AP28742359.csv']

# choose XX random files to analyze
rand_files = np.random.randint(0, N, 3)
hjd, mag, errmag, names = loadmultdata(rand_files)

hjd
mag
# Plot LCs
fig, ax = plt.subplots(len(hjd), 1, figsize=(15,16))
   
# Make common labels
fs = 15 # fontsize
fig.text(0.5, 0.04, 'HJD', ha='center', fontsize=fs)
fig.text(0.07, 0.5, 'Magnitude (V)', va='center', rotation='vertical', fontsize=fs)
# Plot all the datasets
for j in range(len(hjd)):
    ax[j].plot(hjd[j], np.array(mag[j], dtype='float'), '.k')
    ax[j].errorbar(hjd[j], np.array(mag[j], dtype='float'), errmag[j],fmt='k.', ecolor='.7', markersize=4)
    ax[j].text(0.87, 0.05, names[j],transform=ax[j].transAxes,fontsize=15)
     
# plt.savefig('LC1.png', dpi=100, bbox_inches='tight')
plt.savefig('/home/dario/AstronomyLeiden/FRP/MRP-ASAS-Rings/LC3.png', dpi=100, bbox_inches='tight')
plt.show()
#%% 
df = pd.read_csv(files[0])
df.columns
files[0][:-4]
