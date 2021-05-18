"""
Created on Tue Mar 16 16:34:11 2021
@author: dario

Title: Transform CSV lightcurves to PNG plots for the GUI
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
plt.ioff() # disable plt plotting automatically to save CPU
plt.rc('lines', linewidth=2, color='r')  # Syntax 1

path = '/home/dario/AstronomyLeiden/FRP/MRP-ASAS-Rings/zv2/'
files = os.listdir(path)[0:6] # INPUT: Number of LCs to transform

outpath = '/home/dario/AstronomyLeiden/FRP/MRP-ASAS-Rings/lc_plots/'
for file in files:
    df = pd.read_csv(path+file, delimiter=';')
    name = file[2:-4]
    fig, ax = plt.subplots(figsize=(12,6))
        
    norm = np.median(df['flux'])
    ax.errorbar(df['jd'], df['flux']/norm, df['flux_err']/norm,
                fmt='.k', elinewidth=1, alpha=0.8,
                label=name)
    # plt.tight_layout()
    ax.legend()
    ax.set_xlabel('Time [HJD]')
    ax.set_ylabel('Normalised Flux')
    # plt.show()
    plt.savefig(outpath+name+'.png', dpi=90, bbox_inches='tight')
    fig.clear()
    plt.close(fig)
    
    