"""
Created on Tue Feb 16 22:21:50 2021
@author: dario

Title: Plot target coordinates from Zari's catalog
"""
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy import units as u
from astropy.coordinates import SkyCoord

filename = '/home/dario/AstronomyLeiden/FRP/pms.fits'
with fits.open(filename) as hdul:
    hdul.info()
    data = hdul[1].data
 
# Galactic coordinates
lon = data['l']
lat = data['b']

# Define SkyCoord object
c_gal = SkyCoord(lon*u.deg,lat*u.deg, frame='galactic')

# Transform to ICRS frame
c_gal_icrs = c_gal.icrs

# Save coordinates as .txt
np.savetxt('eqcoords.txt', np.transpose([c_gal_icrs.ra, c_gal_icrs.dec]), fmt='%.8f')
len(np.loadtxt('eqcoords.txt'))


# Transform to radians (plotting purposes)
ra_rad = c_gal_icrs.ra.wrap_at(180 * u.deg).radian
dec_rad = c_gal_icrs.dec.radian


# %%
# Plot targets' coordinates
def plot_coord(data, coord_frame, save=False):
    fig = plt.figure(figsize=(12,6))
    ax = fig.add_subplot(111, projection="mollweide")
    plt.title(coord_frame+" Coordinates")
    plt.grid(True)
    ax.plot(data[0], data[1], 'o', markersize=0.7, alpha=0.3)
    # ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'])
    # plt.subplots_adjust(top=0.95,bottom=0.0)
    if save == True:
        plt.savefig('cat_' + coord_frame +'.png', dpi=100)
    plt.show()

# Call function to plot both frames    
galactic = plot_coord([c_gal.l.wrap_at('180d').radian, c_gal.b.radian], 'Galactic') 
eq = plot_coord([ra_rad, dec_rad], 'Equatorial')  



