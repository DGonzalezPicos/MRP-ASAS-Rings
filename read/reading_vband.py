"""
Created on Sat May  8 18:05:37 2021
@author: dario

Title: Read V-band LCs
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import glob, os
from astropy.io import ascii 

path = '/home/dario/AstronomyLeiden/FRP/leiden_vband/lc/camfix/'
filelist = glob.glob(path+'*.dat')


file = filelist[252]




data = ascii.read(file, delimiter='\t')
df = data.to_pandas()
df['flux_err'] = df['flux_err'].replace(99.99, np.nan)
# df = df.replace(99.99, np.nan)
df = df.dropna()
name = file.split('/')[-1][:-4]

df['HJD']

# df = pd.DataFrame(data, columns=['hjd','flux','errflux','asas_sn_id'])

# # df = df.replace('>99.99', np.nan)
# df = df.dropna()

print(df.head())
