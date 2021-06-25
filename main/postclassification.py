"""
Created on Tue Jun 22 20:10:10 2021
@author: dario

Title: Post classification analysis
Read the CSV classification stats from zooniverse

"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from panoptes_aggregation.csv_utils import unjson_dataframe
import os
import plot_functions as plot
import support_functions as supp
import matplotlib as mpl
mpl.rcParams['font.size'] = 14

zoo_path = '/home/dario/AstronomyLeiden/FRP/zoo/'
# the `data.*` columns are read in as strings instead of arrays

file_stats = 'classifying-classifications.csv'
print('Reading classification stats... {:}'.format(file_stats))
data = pd.read_csv(os.path.join(zoo_path, file_stats))

class_types = [] # extract the classification types
for i in range(len(data)):
    class_types.append(str(data['annotations'][i].split("\"")[-6]))
    

    
data['types'] = class_types # add column to original dataframe

    
candidates = data[data['types'] == 'RINGCANDIDATE'] # select only RINGCANDIDATES 


# use unjson_dataframe to convert them to lists
# all values are updated in place leaving null values untouched
# unjson_dataframe(data)
# data['subject_ids']

gaia_ids_rings = []
for i in range(len(candidates)):
    # class_types.append(str(data['annotations'][i]).split("\"")[-6]) # subject ID
    # if len(data['subject_data'].iloc[i].split('\"')[-2]) > 15:
    try: # pick only the ones with gaia_id
        gaia_ids_rings.append(int(candidates['subject_data'].iloc[i].split('\"')[-2]))
    except ValueError:
        continue
    
gaia_ids_rings = list(dict.fromkeys(gaia_ids_rings)) # REMOVE DUPLICATES
#%% Plot and Save list of 'gaia_ids_rings'

v_path = '/home/dario/AstronomyLeiden/FRP/leiden_vband/camfix/'
lco = supp.LightCurveSet(v_path, 'dat', 'v')

rand = np.random.randint(0, len(gaia_ids_rings))
gaia_id = gaia_ids_rings[rand]


for gaia_id in gaia_ids_rings[0:2]:
        fig, ax = plt.subplots(1, figsize=(12,6))
        _ = plot.plot_gaiaID(gaia_id, lco, ax=ax)
        
        # fig.clear()
        # plt.close(fig) #### ENABLE WHEN PLOTTING SEVERAL 
    
#%% SINGLE GAIA_ID

rand = np.random.randint(0, len(gaia_ids_rings)) # plot random lcs from list
gaia_id = gaia_id = gaia_ids_rings[rand] # manual input
gaia_id = 6034203483517407616

df = lco.data(os.path.join(v_path, str(gaia_id)+'.dat'))

fig, ax = plt.subplots(1, figsize=(16,6))
_ = plot.plot_gaiaID(gaia_id, lco, ax=ax)

    
#%% Simbad query
from astroquery.simbad import Simbad
from astropy.io import ascii

id_list = []
for gaia_id in gaia_ids_rings:
    id_list.append('GAIA DR2 '+str(gaia_id))
    
query = Simbad.query_objects(id_list)
query['GAIA_ID'] = gaia_ids_rings
#%% just print without querying again
query[['MAIN_ID', 'RA','DEC', 'GAIA_ID']].pprint_all()
# ascii.write(query, 'simbad_query_candidates.csv', format='csv')
ascii.write(query, 'candidate_list.txt')
query_df = query.to_pandas()
#%% Load rosetta.out for ASAS_SN_ID <<----->> GAIA_ID

ros = np.loadtxt('rosetta.out')
df_ros = pd.DataFrame(ros_array, columns=['asas_sn_id', 'ra', 'dec', 'gaia_id'])

# convert to int
df_ros['gaia_id'] = df_ros['gaia_id'].astype('int')
df_ros['asas_sn_id'] = df_ros['asas_sn_id'].astype('int')

df_ros[df_ros['gaia_id']==3016935867962755968]

#%%
# Read short list of candidates and upload it to zooniverse
short_file = '/home/dario/AstronomyLeiden/FRP/MRP-ASAS-Rings/ring_candidates/short_candidate_list.txt'
short = ascii.read(short_file, format='no_header', delimiter='\t')

ids, ras, decs, comments = ([] for i in range(4))
for row in range(len(short)):

    comments.append(str(short[row]).split(' ')[-1])
    ids.append(str(short[row]).split(' ')[-7])
    
    if row == 2: # weird row formatting
        ra_temp = str(short[row]).split(' ')[-16:-13]
        dec_temp = str(short[row]).split(' ')[-12:-9]
        
    else:
        ra_temp = str(short[row]).split(' ')[-15:-12]
        dec_temp = str(short[row]).split(' ')[-11:-8]
        
    ra_merged = [' '.join(ra_temp)]
    ras.append(ra_merged[0])
    dec_merged = [' '.join(dec_temp)]
    decs.append(dec_merged[0])

df_short = pd.DataFrame(np.transpose([ids,ras,decs,comments]), columns=['GAIA_ID', 'RA', 'DEC','COMMENTS'])
df_short.to_csv('short_list_df.csv', index=False)

df = lco.data(os.path.join(v_path, str(gaia_id)+'.dat'))

#%% Generate new plots for the short-list candidates
outlist = [] # save to upload to zooniverse
for target in df_short.values[:-1]:
    gaia_id, ra, dec, comment = target
    
    fig, ax = plt.subplots(1, figsize=(16,6))
    _ = plot.plot_gaiaID(gaia_id, lco, ax=ax)
    ax.text(x=0.06, y=1.06,s='RA', transform=ax.transAxes)
    ax.text(x=0.19, y=1.06,s='DEC', transform=ax.transAxes)
    ax.text(x=0.02, y=1.01,s=ra, transform=ax.transAxes)
    ax.text(x=0.15, y=1.01,s=dec, transform=ax.transAxes)
    plt.show()
    
    outname = '/home/dario/AstronomyLeiden/FRP/MRP-ASAS-Rings/ring_candidates/short/' + gaia_id +'.png'
    outlist.append(outname)
    plt.savefig(outname, dpi=200, bbox_inches='tight')
    print('Saved {:}'.format(outname))


#%% Finally, upload the LCs to zooniverse
import panoptes as pan
_ = pan.panoptes_add(outlist, 'short_list')

#%%
# =============================================================================
#                       TESS QUERY 
# =============================================================================
