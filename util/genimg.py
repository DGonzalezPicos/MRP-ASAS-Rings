"""
Created on Thu Feb 25 16:26:47 2021
@author: dario

Title: Generate images from LCs
"""
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os
from plotLC import get_image, get_lightcurve
plt.ioff()
basepath = '/home/dario/AstronomyLeiden/FRP/ASAS' # path to images
os.chdir(basepath)
file = 'asassn-catalog.csv'
print('Reading ASAS-SN catalog...')
df = pd.read_csv(file, low_memory=False, nrows=2000)
len(df)
df.columns
df['asassn_name']
df['id']
set(df['Type'])

df['Classified']

file = 'lc' + str(df['id'][0]) + '.dat'
hjd, mag, errmag, _ = get_lightcurve(file)
lc = [hjd, mag, errmag]
img = img_from_lc(lc)
plt.imshow(img, cmap='gray')
plt.axis('off')
plt.show()
#%%
objs = list(zip(df['id'], df['Type'], df['class_probability']))

# Get objects that are well classified
# i.e "class_probability=1"
X, y = ([] for _ in range(2))
for name, obj_type, value in objs:
    if value > 0.99:
        X.append(100000+np.where(name)[0][0])
        y.append(obj_type)
#%%
from sklearn.model_selection import train_test_split

# X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.33, random_state=42)

# out = '/home/dario/AstronomyLeiden/FRP/MRP-ASAS-Rings/images/nn'
# print('Generating images...')
# for file in X_train:
#     file = 'lc' + str(file) + '.dat'
#     fig = get_image(file, outpath=out)

## TODO: Convert image to np.array instead of saving as png

#%%
# import imageio

# img = imageio.imread('./img_lc_test.png')

# plt.imshow(img)
# plt.axis('off')
import io
x = np.linspace(1,100,100)
y = np.random.normal(size=100)

fig = plt.figure(figsize=(4,3))
plt.plot(x,y,'.k')
plt.axis('off')
plt.show()
plt.tight_layout()
fig.bbox.bounds

io_buf = io.BytesIO()
fig.savefig(io_buf, format='raw', dpi=200, pad_inches=0, transparent=True)
io_buf.seek(0)
img_arr = np.reshape(np.frombuffer(io_buf.getvalue(), dtype=np.uint8),
                     newshape=(int(fig.bbox.bounds[3]), int(fig.bbox.bounds[2]), -1))

#%%
img_arr.shape
gray_img_arr = np.mean(img_arr, -1)
gray_img_arr.shape
plt.imshow(gray_img_arr, cmap='gray')
plt.axis('off')
plt.show()
io_buf.close()

#%%

def img_from_lc(lc, size=(3,5), dpi=200):
    fig, ax = plt.subplots(figsize=size)
    ax.errorbar(lc[0],lc[1],lc[2],fmt='.k')
    ax.axis('off')
    # plt.show()
    plt.tight_layout()
    fig.bbox.bounds
    io_buf = io.BytesIO()
    fig.savefig(io_buf, format='raw', dpi=dpi, pad_inches=0, transparent=True)
    fig.clear()
    io_buf.seek(0)
    img_arr = np.reshape(np.frombuffer(io_buf.getvalue(), dtype=np.uint8),
                         newshape=(int(fig.bbox.bounds[2]), int(fig.bbox.bounds[3]), -1))
    gray_img = np.mean(img_arr, -1)
    return gray_img





