"""
Created on Thu Apr 22 14:13:47 2021
@author: dario

Title: File management
"""
import os, glob

path = '/home/dario/AstronomyLeiden/FRP/gband_full/'

filelist = glob.glob(path+'*.dat')
#%%
for file in filelist[:200]:
    name = file.split('/')[-1][:-4]
    print(name, len(name))
    if len(name) < 16:
        os.remove(file)
