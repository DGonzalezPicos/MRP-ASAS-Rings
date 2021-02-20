"""
Created on Fri Feb 19 14:13:40 2021
@author: dario

Title: Cross-matching targets from two catalogs (and coordinate plot)
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units as u
# Load ASAS catalog
file = '/home/dario/AstronomyLeiden/FRP/ASAS/asassn-catalog.csv'
print('Reading ASAS-SN catalog...')
df = pd.read_csv(file, low_memory=False)

name, gaia_id= df['id'], df['gdr2_id']
ra, dec = df['raj2000'], df['dej2000']
l, b = df['l'], df['b']

print('---> {:} objects loaded'.format(len(df['id'])))
# Get only the objects with Gaia ID (remove NaN)
asas_gaia_id = gaia_id.dropna()
print('-------> Only {:} objects ({:.2f}%) have Gaia ID'.format(len(asas_gaia_id), 100*len(asas_gaia_id)/len(gaia_id)))
print('-----------------------------------')

# Load Zari's catalog
filename = '/home/dario/AstronomyLeiden/FRP/pms.fits'
with fits.open(filename) as hdul:
    print('Reading Zari\'s catalog...')
    # hdul.info()
    data = hdul[1].data
 
# Galactic coordinates
name_id = data['source_id']
l, b = data['l'], data['b']
print('---> {:} objects loaded'.format(len(name_id)))
print('-----------------------------------')

zari_id = pd.Series(name_id,dtype='float64') # unify format
print('Looking for targets in both catalogs...')
intersect = pd.Series(list(set(asas_gaia_id).intersection(set(zari_id))))
print('----->> Found {:} targets in common'.format(len(intersect)))

#%%
# Get the coordinates of the targets

gal_c, eq_c = ([] for i in range(2))
for obj in asas_gaia_id.unique():
    if obj in intersect.values:
        temp = intersect[intersect==obj].index
        index = int(temp.values)
        gal_c.append([l[index], b[index]])
        eq_c.append([ra[index], dec[index]])
        
gal_c = np.transpose(gal_c)        
gal_coord = SkyCoord(gal_c[0]*u.deg, gal_c[1]*u.deg, frame = 'galactic')
#%%
# Plot the targets from the catalog intersection
def plot_coord(data, coord_frame, save=False):
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111, projection="mollweide")
    plt.title(coord_frame+" Coordinates")
    plt.grid(True)
    ax.plot(data[0], data[1], 'o', markersize=1, alpha=0.3)
    # ax.set_xticklabels(['14h','16h','18h','20h','22h','0h','2h','4h','6h','8h','10h'])
    # plt.subplots_adjust(top=0.95,bottom=0.0)
    if save == True:
        plt.savefig('/home/dario/AstronomyLeiden/FRP/MRP-ASAS-Rings/images/cat_' + coord_frame +'.png', 
                    dpi=100, bbox_inches='tight')
    plt.show()


galactic = plot_coord([gal_coord.l.wrap_at('180d').radian, gal_coord.b.radian], 'CrossMatchGalactic', save=False) 


#%%
# Save the target's name and coordinates
# from astropy.io import ascii
# from astropy.table import Table
# data_to_write = Table({'id':intersect,
#                        'l': gal_coord.l,
#                        'b': gal_coord.b})

# ascii.write(data_to_write, 'variable_targets.txt', overwrite=True)

# csv format is prefered over ascii
df_var = pd.DataFrame({'id':intersect,
                        'l': gal_coord.l,
                        'b': gal_coord.b})
output='var_targets.csv'
print('Writing resulting catalog to "{:}"'.format(output))
df_var.to_csv(output)

