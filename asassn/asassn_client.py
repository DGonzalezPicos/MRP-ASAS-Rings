"""
Created on Tue May 18 09:17:04 2021
@author: dario

Title: Basic ASAS-SN client query
"""
import numpy as np
from pyasassn.client import SkyPatrolClient
import glob


print('------------------------------')
print('Loading SkyPatrolClient...')
print('------------------------------')

client = SkyPatrolClient("bad_client", "a5a55N_CLIENT")
client.catalogs.stellar_main.head(12)

# print('Available catalogs {:}'.format(client.catalogs.catalog_names()))

# Gaia_IDs from already downloaded objects (worked in the past)
# gaia_id = 3207233877299221504 # ValueError: No objects to concatenate

# Additional IDs to test
# gaia_id = 5789045856890265216 # this one works 
# gaia_id = 2032226624418331136 # ValueError: No objects to concatenate
# gaia_id = 5994772213404409216 # this one works 
gaia_id = 3342413609255935872 # ValueError: No objects to concatenate

# If "download=False"  only retrieve [ASAS_SN_ID, RA, DEC, GAIA_ID]
lc_info = client.query_list([gaia_id], catalog='stellar_main', id_col='gaia_id', download=False)
print('Info downloaded for {:}'.format(gaia_id))

# If "download=True" retrieve all light curve data 
lc_full = client.query_list([gaia_id], catalog='stellar_main', id_col='gaia_id', download=True)
print('Data downloaded for {:}'.format(gaia_id))

#%%
"""
Created on Tue May 18 09:17:04 2021
@author: dario

Title: Basic ASAS-SN client query
"""
import numpy as np
from pyasassn.client import SkyPatrolClient
import glob


print('------------------------------')
print('Loading SkyPatrolClient...')
print('------------------------------')

client = SkyPatrolClient("bad_client", "a5a55N_CLIENT")
client.catalogs.stellar_main.head(15)
# Read file list from Tharindu's V-band lightcurves (only first time)
try:
    _ = len(gaia_ids)
    
except NameError:
    v_path = '/home/dario/AstronomyLeiden/FRP/leiden_vband/lc/camfix/'
    filelist = glob.glob(v_path+'*.dat')
    gaia_ids = []
    
    print('Loading gaia_ids from {:}'.format(v_path))
    for file in filelist:
        name = file.split('/')[-1][:-4]
        gaia_ids.append(int(name))


id_error = []
query = gaia_ids[10:14]
for gaia_id in query:
    lc_info = client.query_list([gaia_id], catalog='stellar_main', id_col='gaia_id', download=False)
    print('Info downloaded for {:}'.format(gaia_id))
    
    
    try:
        lc_full = client.query_list([gaia_id], catalog='stellar_main', id_col='gaia_id', download=True)
        print('Data downloaded for {:}'.format(gaia_id))
        print('Number of epochs {:}'.format(len(lc_full.data)))
        lc_full.data['gaia_id'] = np.full_like(lc_full.data['jd'], lc_info['gaia_id'])
        print(lc_full.data.tail())
    except ValueError:
        id_error.append(gaia_id)
        print('ValueError')

print('Retrieval succes {:.2f} %'.format(100*len(id_error)/len(query)))
