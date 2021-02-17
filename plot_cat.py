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

ra_rad = c_gal_icrs.ra.wrap_at(180 * u.deg).radian
dec_rad = c_gal_icrs.dec.radian

#%%

# Plot targets in galactic coordinates
plt.figure(figsize=(12,6))
plt.subplot(111, projection="aitoff")
plt.title("Galactic Coordinates")
plt.grid(True)
plt.plot(c_gal.l.wrap_at('180d').radian, c_gal.b.radian, 'o', markersize=0.7, alpha=0.3)

plt.subplots_adjust(top=0.95,bottom=0.0)
plt.savefig('cat_gal.png', dpi=100)
plt.show()

# %%
# Plot targets in equatorial coordinates
fig = plt.figure(figsize=(12,6))
ax = fig.add_subplot(111, projection="mollweide")
plt.title("Equatorial Coordinates")
plt.grid(True)
ax.plot(ra_rad, dec_rad, 'o', markersize=0.7, alpha=0.3)
ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'])
# plt.subplots_adjust(top=0.95,bottom=0.0)
plt.savefig('cat_eq.png', dpi=100)
plt.show()