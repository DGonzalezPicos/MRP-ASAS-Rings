"""
Created on Sun Feb 14 23:24:46 2021
@author: dario

Title: Plotting lightcurves from rnadom LCs from the Full Catalog of
 Variable ASAS-SN sources 
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
    hjd, mag, errmag, names = [[] for i in range(4)]
    for i in rand_files:
        data = ascii.read(file_list[i])
        hjd.append(data['col1'])
        errmag.append(data['col4'])
        
        temp_list = []
        for item in data['col3']:
            a  = str(type(item))
            if a == "<class 'numpy.str_'>":
                temp_list.append(float(item.split('>')[-1]))
            else:
                temp_list.append(item)
        
        mag.append(temp_list)
        
        
        
        name = file_list[i].split('/')[-1][2:-4]
        names.append(name)
    return hjd, mag, errmag, names



# Call function
# choose XX random files to analyze
rand_files = np.random.randint(0, N, 4)
hjd, mag, errmag, names = loadmultdata(rand_files)

errmag 
# Fix "infinite values" for the "errmag" data by
# replacing them with the median

for j in range(len(errmag)):
    med = np.median(errmag[j])
    if med > 90:
        print(f'Error bars not valid for plot {names[j]}')
        
    for count, k in enumerate(errmag[j]):
        if k > 90: # the wrong values have an assigned "99.99"
            if med < 90:
                errmag[j][count] = np.median(errmag[j])
            else:
                errmag[j][count] = 0.0 # no value
            
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

