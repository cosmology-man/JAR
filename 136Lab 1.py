#!/usr/bin/env python
# coding: utf-8

# In[1]:


#imports
from __future__ import print_function, absolute_import, division, unicode_literals
import numpy as np
import glob, os, sys
import pdb
from matplotlib import pyplot as plt
from scipy.stats import norm
import statistics 
from astropy.io import fits
from itertools import chain


# In[2]:


#testing
bias0 = fits.open('cleanbias01.fits')
#bias0.info()
bias0.verify('fix')
data0 = bias0[0].data
print(data0.shape)

print(data0[0][0])
print(data0[0])
print(data0)
plt.figure(figsize=(10,10))
plt.imshow(data0, vmax= 1100)
plt.colorbar()


print(np.min(data0))
print(np.max(data0))


# In[3]:


#importing all bias files and organizing data
bi1 = fits.open('cleanbias01.fits')
bi2 = fits.open('cleanbias02.fits')
bi3 = fits.open('cleanbias03.fits')
bi4 = fits.open('cleanbias04.fits')
bi5 = fits.open('cleanbias05.fits')
bi6 = fits.open('cleanbias06.fits')
bi7 = fits.open('cleanbias07.fits')
bi8 = fits.open('cleanbias08.fits')
bi9 = fits.open('cleanbias09.fits')
bi10 = fits.open('cleanbias10.fits')

bi1.verify('fix')
bi2.verify('fix')
bi3.verify('fix')
bi4.verify('fix')
bi5.verify('fix')
bi6.verify('fix')
bi7.verify('fix')
bi8.verify('fix')
bi9.verify('fix')
bi10.verify('fix')
#bi2.info()

a = bi1[0].data
b = bi2[0].data
c = bi3[0].data
d = bi4[0].data
e = bi5[0].data
f = bi6[0].data
g = bi7[0].data
h = bi8[0].data
m = bi9[0].data
n = bi10[0].data

print(a[0][0])
print(a[0])
print(a)
print(len(a[0]))
print(len(a))


# In[4]:


#median of bias fits
masterbias = np.ndarray((len(a), len(a[0])))
for j in range(len(a[0])):
    for i in range(len(a)):
        temp = []
        temp = np.append(temp, a[i][j])
        temp = np.append(temp, b[i][j])
        temp = np.append(temp, c[i][j])
        temp = np.append(temp, d[i][j])
        temp = np.append(temp, e[i][j])
        temp = np.append(temp, f[i][j])
        temp = np.append(temp, g[i][j])
        temp = np.append(temp, h[i][j])
        temp = np.append(temp, m[i][j])
        temp = np.append(temp, n[i][j])
        t = np.median(temp)
        masterbias[i][j] = t
    
print(masterbias)
print(masterbias.shape)

plt.figure(figsize= (10,10))
plt.imshow(masterbias)
plt.colorbar()
plt.title("Master bias")
plt.savefig("master bias.pdf")
#master bias of cleaned bias images - overscan already removed


# In[5]:


#import B Twilight Flats
fb1 = fits.open('cleanTwilightFlatB01.fits')
fb2 = fits.open('cleanTwilightFlatB02.fits')
fb3 = fits.open('cleanTwilightFlatB03.fits')
fb4 = fits.open('cleanTwilightFlatB04.fits')
fb5 = fits.open('cleanTwilightFlatB05.fits')

#import V Twilight Flats
fv1 = fits.open('cleanTwilightFlatV06.fits')
fv2 = fits.open('cleanTwilightFlatV07.fits')
fv3 = fits.open('cleanTwilightFlatV08.fits')
fv4 = fits.open('cleanTwilightFlatV09.fits')
fv5 = fits.open('cleanTwilightFlatV10.fits')

#import R Twilight Flats
fr1 = fits.open('cleanTwilightFlatR11.fits')
fr2 = fits.open('cleanTwilightFlatR12.fits')
fr3 = fits.open('cleanTwilightFlatR13.fits')
fr4 = fits.open('cleanTwilightFlatR14.fits')
fr5 = fits.open('cleanTwilightFlatR15.fits')

fb1.verify('fix')
fb2.verify('fix')
fb3.verify('fix')
fb4.verify('fix')
fb5.verify('fix')
fv1.verify('fix')
fv2.verify('fix')
fv3.verify('fix')
fv4.verify('fix')
fv5.verify('fix')
fr1.verify('fix')
fr2.verify('fix')
fr3.verify('fix')
fr4.verify('fix')
fr5.verify('fix')
#fv2.info()

bflat1 = fb1[0].data
bflat2 = fb2[0].data
bflat3 = fb3[0].data
bflat4 = fb4[0].data
bflat5 = fb5[0].data
vflat1 = fv1[0].data
vflat2 = fv2[0].data
vflat3 = fv3[0].data
vflat4 = fv4[0].data
vflat5 = fv5[0].data
rflat1 = fr1[0].data
rflat2 = fr2[0].data
rflat3 = fr3[0].data
rflat4 = fr4[0].data
rflat5 = fr5[0].data
print(len(bflat1[0]))
print(len(bflat1))


# In[6]:


#creating median array
bflatmed = np.ndarray((len(bflat1), len(bflat1[0])))
for j in range(len(bflat1[0])):
    for i in range(len(bflat1)):
        temp = []
        temp = np.append(temp, bflat1[i][j])
        temp = np.append(temp, bflat2[i][j])
        temp = np.append(temp, bflat3[i][j])
        temp = np.append(temp, bflat4[i][j])
        temp = np.append(temp, bflat5[i][j])
        t = np.median(temp)
        bflatmed[i][j] = t
    
print(bflatmed)
print(bflatmed.shape)

plt.figure(figsize= (10,10))
plt.imshow(bflatmed)
plt.colorbar()
plt.title("Median of B Twilight Flat fits")
plt.savefig("master b flat.pdf")


# In[7]:


#creating median array
vflatmed = np.ndarray((len(vflat1), len(vflat1[0])))
for j in range(len(vflat1[0])):
    for i in range(len(vflat1)):
        temp = []
        temp = np.append(temp, vflat1[i][j])
        temp = np.append(temp, vflat2[i][j])
        temp = np.append(temp, vflat3[i][j])
        temp = np.append(temp, vflat4[i][j])
        temp = np.append(temp, vflat5[i][j])
        t = np.median(temp)
        vflatmed[i][j] = t
    
print(vflatmed)
print(vflatmed.shape)

plt.figure(figsize= (10,10))
plt.imshow(vflatmed)
plt.colorbar()
plt.title("Median of V Twilight Flat fits")
plt.savefig("master v flat.pdf")


# In[8]:


#creating median array
rflatmed = np.ndarray((len(rflat1), len(rflat1[0])))
for j in range(len(rflat1[0])):
    for i in range(len(rflat1)):
        temp = []
        temp = np.append(temp, rflat1[i][j])
        temp = np.append(temp, rflat2[i][j])
        temp = np.append(temp, rflat3[i][j])
        temp = np.append(temp, rflat4[i][j])
        temp = np.append(temp, rflat5[i][j])
        t = np.median(temp)
        rflatmed[i][j] = t
    
print(rflatmed)
print(rflatmed.shape)

plt.figure(figsize= (10,10))
plt.imshow(rflatmed)
plt.colorbar()
plt.title("Median of R Twilight Flat fits")
plt.savefig("master r flat.pdf")


# In[9]:


#create unbiased flats: remove master bias from all master flats
unbiasedbflat= np.subtract(bflatmed, masterbias)
unbiasedvflat= np.subtract(vflatmed, masterbias)
unbiasedrflat= np.subtract(rflatmed, masterbias)
#print(unbiasedbflat)
#print(unbiasedbflat[0])
#print(unbiasedvflat)
#print(unbiasedvflat[0])
#print(unbiasedrflat)
#print(unbiasedrflat[0])


# In[10]:


#plot unbiased b flat
plt.figure(figsize= (10,10))
plt.imshow(unbiasedbflat)
plt.colorbar()
plt.title("Unbiased B Twilight Flat fits")
plt.savefig("unbiased master b flat.pdf")


# In[11]:


#plot unbiased v flat
plt.figure(figsize= (10,10))
plt.imshow(unbiasedvflat)
plt.colorbar()
plt.title("Unbiased V Twilight Flat fits")
plt.savefig("unbiased master v flat.pdf")


# In[12]:


#plot unbiased r flat
plt.figure(figsize= (10,10))
plt.imshow(unbiasedrflat)
plt.colorbar()
plt.title("Unbiased R Twilight Flat fits")
plt.savefig("unbiased master r flat.pdf")

print(len(unbiasedrflat))
print(len(unbiasedrflat[0]))
print(unbiasedrflat)
print(unbiasedrflat[0])


# In[13]:


#normalize the master flat to its mean (which means dividing each pixel by the pixel array mean value) ->plot and save
bnormflat = np.ndarray((len(unbiasedbflat), len(unbiasedbflat[0])))

for i in range(len(unbiasedbflat)):
    fmean = np.mean(unbiasedbflat[i])
    for j in range(len(unbiasedbflat[0])):
        bnormflat[i][j]= (unbiasedbflat[i][j])/fmean

print(bnormflat)
print(bnormflat[0])


plt.figure(figsize= (10,10))
plt.imshow(bnormflat)
plt.colorbar()
plt.title("Normalization of B Twilight Flat fits")
plt.savefig("normalized unbiased master b flat.pdf")


# In[14]:


#normalize the master flat to its mean (which means dividing each pixel by the pixel array mean value) ->plot and save
vnormflat = np.ndarray((len(unbiasedvflat), len(unbiasedvflat[0])))

for i in range(len(unbiasedvflat)):
    fmean = np.mean(unbiasedvflat[i])
    for j in range(len(unbiasedvflat[0])):
        vnormflat[i][j]= (unbiasedvflat[i][j])/fmean

print(vnormflat)
print(vnormflat[0])


plt.figure(figsize= (10,10))
plt.imshow(vnormflat)
plt.colorbar()
plt.title("Normalization of V Twilight Flat fits")
plt.savefig("normalized unbiased master v flat.pdf")


# In[15]:


#normalize the master flat to its mean (which means dividing each pixel by the pixel array mean value) ->plot and save
rnormflat = np.ndarray((len(unbiasedrflat), len(unbiasedrflat[0])))

for i in range(len(unbiasedrflat)):
    fmean = np.mean(unbiasedrflat[i])
    for j in range(len(unbiasedrflat[0])):
        rnormflat[i][j]= (unbiasedrflat[i][j])/fmean

print(rnormflat)
print(rnormflat[0])


plt.figure(figsize= (10,10))
plt.imshow(rnormflat)
plt.colorbar()
plt.title("Normalization of R Twilight Flat fits")
plt.savefig("normalized unbiased master r flat.pdf")


# In[16]:


#import science fits, organize data
b_szher = fits.open('cleanHerB.fits')
b_m15 = fits.open('cleanM15B.fits')
b_m29 = fits.open('cleanM29B.fits')
v_szher = fits.open('cleanHerV.fits')
v_m15 = fits.open('cleanM15V.fits')
v_m29 = fits.open('cleanM29V.fits')
r_szher = fits.open('cleanHerR.fits')
r_m15 = fits.open('cleanM15R.fits')
r_m29 = fits.open('cleanM29R.fits')

b_szher.verify('fix')
b_m15.verify('fix')
b_m29.verify('fix')
v_szher.verify('fix')
v_m15.verify('fix')
v_m29.verify('fix')
r_szher.verify('fix')
r_m15.verify('fix')
r_m29.verify('fix')

bher = b_szher[0].data
vher = v_szher[0].data
rher = r_szher[0].data
bm15 = b_m15[0].data
vm15 = v_m15[0].data
rm15 = r_m15[0].data
bm29 = b_m29[0].data
vm29 = v_m29[0].data
rm29 = r_m29[0].data

print(bher)
print(bher.shape)


# In[17]:


#subtract bias from science image then divide science image by normalized unbiased flat for final image
unbiasedbher = np.subtract(bher,masterbias)
unbiasedvher = np.subtract(vher,masterbias)
unbiasedrher = np.subtract(rher,masterbias)
unbiasedbm15 = np.subtract(bm15,masterbias)
unbiasedvm15 = np.subtract(vm15,masterbias)
unbiasedrm15 = np.subtract(rm15,masterbias)
unbiasedbm29 = np.subtract(bm29,masterbias)
unbiasedvm29 = np.subtract(vm29,masterbias)
unbiasedrm29 = np.subtract(rm29,masterbias)


finalbher = np.ndarray((len(bher), len(bher[0])))
finalvher = np.ndarray((len(vher), len(vher[0])))
finalrher = np.ndarray((len(rher), len(rher[0])))
finalbm15 = np.ndarray((len(bm15), len(bm15[0])))
finalvm15 = np.ndarray((len(vm15), len(vm15[0])))
finalrm15 = np.ndarray((len(rm15), len(rm15[0])))
finalbm29 = np.ndarray((len(bm29), len(bm29[0])))
finalvm29 = np.ndarray((len(vm29), len(vm29[0])))
finalrm29 = np.ndarray((len(rm29), len(rm29[0])))

for i in range(len(finalbher)):
    finalbher[i] = np.true_divide(unbiasedbher[i], bnormflat[i])
    finalvher[i] = np.true_divide(unbiasedvher[i], vnormflat[i])
    finalrher[i] = np.true_divide(unbiasedrher[i], rnormflat[i])
    finalbm15[i] = np.true_divide(unbiasedbm15[i], bnormflat[i])
    finalvm15[i] = np.true_divide(unbiasedvm15[i], vnormflat[i])
    finalrm15[i] = np.true_divide(unbiasedrm15[i], rnormflat[i])
    finalbm29[i] = np.true_divide(unbiasedbm29[i], bnormflat[i])
    finalvm29[i] = np.true_divide(unbiasedvm29[i], vnormflat[i])
    finalrm29[i] = np.true_divide(unbiasedrm29[i], rnormflat[i])
    
print(finalbher[0])
print(finalbher)
print(finalbher.shape)


# In[18]:


plt.figure(figsize= (10,10))
plt.imshow(finalbher, vmin=0, vmax=1000)
plt.colorbar()
plt.title("Final B SZ Her Image")
plt.savefig("Final B SZ Her Image")


# In[19]:


plt.figure(figsize= (10,10))
plt.imshow(finalvher, vmin=0, vmax=1000)
plt.colorbar()
plt.title("Final V SZ Her Image")
plt.savefig("Final V SZ Her Image")


# In[20]:


plt.figure(figsize= (10,10))
plt.imshow(finalrher, vmin=0, vmax=1000)
plt.colorbar()
plt.title("Final R SZ Her Image")
plt.savefig("Final R SZ Her Image")


# In[21]:


plt.figure(figsize= (10,10))
plt.imshow(finalbm15, vmin=0, vmax=30000)
plt.colorbar()
plt.title("Final B M15 Image")
plt.savefig("Final B M15 Image")


# In[22]:


plt.figure(figsize= (10,10))
plt.imshow(finalvm15, vmin=0, vmax=50000)
plt.colorbar()
plt.title("Final V M15 Image")
plt.savefig("Final V M15 Image")


# In[23]:


plt.figure(figsize= (10,10))
plt.imshow(finalrm15, vmin=0, vmax=50000)
plt.colorbar()
plt.title("Final R M15 Image")
plt.savefig("Final R M15 Image")


# In[24]:


plt.figure(figsize= (10,10))
plt.imshow(finalbm29, vmin=0, vmax=5000)
plt.colorbar()
plt.title("Final B M29 Image")
plt.savefig("Final B M29 Image")


# In[25]:


plt.figure(figsize= (10,10))
plt.imshow(finalvm29, vmin=0, vmax=6000)
plt.colorbar()
plt.title("Final V M29 Image")
plt.savefig("Final V M29 Image")


# In[26]:


plt.figure(figsize= (10,10))
plt.imshow(finalrm29, vmin=0, vmax=7000)
plt.colorbar()
plt.title("Final R M29 Image")
plt.savefig("Final R M29 Image")


# In[27]:


pip install sep


# In[ ]:





# In[28]:


import numpy as np
import sep

# additional setup for reading the test image and displaying plots
import matplotlib.pyplot as plt
from matplotlib import rcParams

get_ipython().run_line_magic('matplotlib', 'inline')

rcParams['figure.figsize'] = [10., 8.]


# In[29]:


#M29B

# measure a spatially varying background on the image
bm29bkg = sep.Background(finalbm29)
# get a "global" mean and noise of the image background:
print(bm29bkg.globalback)
print(bm29bkg.globalrms)
# evaluate background as 2-d array, same size as original image
bm29bkg_image = bm29bkg.back()
# bm29bkg_image = np.array(bm29bkg) # equivalent to above
# show the background
plt.figure(figsize= (10,10))
plt.imshow(bm29bkg_image, interpolation='nearest', cmap='gray', origin='lower')
plt.colorbar();
plt.title("M29 B Background")


# In[30]:


# evaluate the background noise as 2-d array, same size as original image
bm29bkg_rms = bm29bkg.rms()
# show the background noise
plt.figure(figsize= (10,10))
plt.imshow(bm29bkg_rms, interpolation='nearest', cmap='gray', origin='lower')
plt.colorbar();
plt.title("M29 B Background RMS")


# In[31]:


# subtract the background
bm29data_sub = finalbm29 - bm29bkg

bm29obj = sep.extract(bm29data_sub, 20, err=bm29bkg.globalrms)
# how many objects were detected
print(len(bm29obj))
#print(bm29obj)

from matplotlib.patches import Ellipse

# plot background-subtracted image
fig, ax = plt.subplots()
m, s = np.mean(bm29data_sub), np.std(bm29data_sub)
im = ax.imshow(bm29data_sub, interpolation='nearest', cmap='gray',
               vmin=m-s, vmax=m+s, origin='lower')

# plot an ellipse for each object
for i in range(len(bm29obj)):
    e = Ellipse(xy=(bm29obj['x'][i], bm29obj['y'][i]),
                width=6*bm29obj['a'][i],
                height=6*bm29obj['b'][i],
                angle=bm29obj['theta'][i] * 180. / np.pi)
    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax.add_artist(e)
    
print(np.max(bm29data_sub))


# In[32]:


#SZHERB

# measure a spatially varying background on the image
bherbkg = sep.Background(finalbher)
# get a "global" mean and noise of the image background:
print(bherbkg.globalback)
print(bherbkg.globalrms)
# evaluate background as 2-d array, same size as original image
bherbkg_image = bherbkg.back()
# bherbkg_image = np.array(bherbkg) # equivalent to above
# show the background
plt.figure(figsize= (10,10))
plt.imshow(bherbkg_image, interpolation='nearest', cmap='gray', origin='lower')
plt.colorbar();
plt.title("SZHER B Background")


# In[33]:


# evaluate the background noise as 2-d array, same size as original image
bherbkg_rms = bherbkg.rms()
# show the background noise
plt.figure(figsize= (10,10))
plt.imshow(bherbkg_rms, interpolation='nearest', cmap='gray', origin='lower')
plt.colorbar();
plt.title("SZHER B Background RMS")


# In[34]:


# subtract the background
bherdata_sub = finalbher - bherbkg

bherobj = sep.extract(bherdata_sub, 10, err=bherbkg.globalrms)
# how many objects were detected
print(len(bherobj))
#print(bherobj)

from matplotlib.patches import Ellipse

# plot background-subtracted image
fig, ax = plt.subplots()
m, s = np.mean(bherdata_sub), np.std(bherdata_sub)
im = ax.imshow(bherdata_sub, interpolation='nearest', cmap='gray',
               vmin=m-s, vmax=m+s, origin='lower')

# plot an ellipse for each object
for i in range(len(bherobj)):
    e = Ellipse(xy=(bherobj['x'][i], bherobj['y'][i]),
                width=6*bherobj['a'][i],
                height=6*bherobj['b'][i],
                angle=bherobj['theta'][i] * 180. / np.pi)
    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax.add_artist(e)


# In[35]:


#M15B

# measure a spatially varying background on the image
bm15bkg = sep.Background(finalbm15)
# get a "global" mean and noise of the image background:
print(bm15bkg.globalback)
print(bm15bkg.globalrms)
# evaluate background as 2-d array, same size as original image
bm15bkg_image = bm15bkg.back()
# bm15bkg_image = np.array(bm15bkg) # equivalent to above
# show the background
plt.figure(figsize= (10,10))
plt.imshow(bm15bkg_image, interpolation='nearest', cmap='gray', origin='lower')
plt.colorbar();
plt.title("M15 B Background")


# In[36]:


# evaluate the background noise as 2-d array, same size as original image
bm15bkg_rms = bm15bkg.rms()
# show the background noise
plt.figure(figsize= (10,10))
plt.imshow(bm15bkg_rms, interpolation='nearest', cmap='gray', origin='lower')
plt.colorbar();
plt.title("M15 B Background RMS")


# In[37]:


# subtract the background
bm15data_sub = finalbm15 - bm15bkg

bm15obj = sep.extract(bm15data_sub, 15, err=bm15bkg.globalrms)
# how many objects were detected
print(len(bm15obj))
#print(bm15obj)

from matplotlib.patches import Ellipse

# plot background-subtracted image
fig, ax = plt.subplots()
m, s = np.mean(bm15data_sub), np.std(bm15data_sub)
im = ax.imshow(bm15data_sub, interpolation='nearest', cmap='gray',
               vmin=m-s, vmax=m+s, origin='lower')

# plot an ellipse for each object
for i in range(len(bm15obj)):
    e = Ellipse(xy=(bm15obj['x'][i], bm15obj['y'][i]),
                width=6*bm15obj['a'][i],
                height=6*bm15obj['b'][i],
                angle=bm15obj['theta'][i] * 180. / np.pi)
    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax.add_artist(e)
    
print(bm15data_sub.shape)


# In[38]:


#M29V

# measure a spatially varying background on the image
vm29bkg = sep.Background(finalvm29)
# get a "global" mean and noise of the image background:
print(vm29bkg.globalback)
print(vm29bkg.globalrms)
# evaluate background as 2-d array, same size as original image
vm29bkg_image = vm29bkg.back()
# vm29bkg_image = np.array(vm29bkg) # equivalent to above
# show the background
plt.figure(figsize= (10,10))
plt.imshow(vm29bkg_image, interpolation='nearest', cmap='gray', origin='lower')
plt.colorbar();
plt.title("M29 V Background")


# In[39]:


# evaluate the background noise as 2-d array, same size as original image
vm29bkg_rms = vm29bkg.rms()
# show the background noise
plt.figure(figsize= (10,10))
plt.imshow(vm29bkg_rms, interpolation='nearest', cmap='gray', origin='lower')
plt.colorbar();
plt.title("M29 V Background RMS")


# In[40]:


# subtract the background
vm29data_sub = finalvm29 - vm29bkg

vm29obj = sep.extract(vm29data_sub, 20, err=vm29bkg.globalrms)
# how many objects were detected
print(len(vm29obj))
#print(vm29obj)

from matplotlib.patches import Ellipse

# plot background-subtracted image
fig, ax = plt.subplots()
m, s = np.mean(vm29data_sub), np.std(vm29data_sub)
im = ax.imshow(vm29data_sub, interpolation='nearest', cmap='gray',
               vmin=m-s, vmax=m+s, origin='lower')

# plot an ellipse for each object
for i in range(len(vm29obj)):
    e = Ellipse(xy=(vm29obj['x'][i], vm29obj['y'][i]),
                width=6*vm29obj['a'][i],
                height=6*vm29obj['b'][i],
                angle=vm29obj['theta'][i] * 180. / np.pi)
    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax.add_artist(e)


# In[41]:


#M15V

# measure a spatially varying background on the image
vm15bkg = sep.Background(finalvm15)
# get a "global" mean and noise of the image background:
print(vm15bkg.globalback)
print(vm15bkg.globalrms)
# evaluate background as 2-d array, same size as original image
vm15bkg_image = vm15bkg.back()
# vm15bkg_image = np.array(vm15bkg) # equivalent to above
# show the background
plt.figure(figsize= (10,10))
plt.imshow(vm15bkg_image, interpolation='nearest', cmap='gray', origin='lower')
plt.colorbar();
plt.title("M15 V Background")


# In[42]:


# evaluate the background noise as 2-d array, same size as original image
vm15bkg_rms = vm15bkg.rms()
# show the background noise
plt.figure(figsize= (10,10))
plt.imshow(vm15bkg_rms, interpolation='nearest', cmap='gray', origin='lower')
plt.colorbar();
plt.title("M15 V Background RMS")


# In[43]:


# subtract the background
vm15data_sub = finalvm15 - vm15bkg

vm15obj = sep.extract(vm15data_sub, 20, err=vm15bkg.globalrms)
# how many objects were detected
print(len(vm15obj))
#print(vm15obj)

from matplotlib.patches import Ellipse

# plot background-subtracted image
fig, ax = plt.subplots()
m, s = np.mean(vm15data_sub), np.std(vm15data_sub)
im = ax.imshow(vm15data_sub, interpolation='nearest', cmap='gray',
               vmin=m-s, vmax=m+s, origin='lower')

# plot an ellipse for each object
for i in range(len(vm15obj)):
    e = Ellipse(xy=(vm15obj['x'][i], vm15obj['y'][i]),
                width=6*vm15obj['a'][i],
                height=6*vm15obj['b'][i],
                angle=vm15obj['theta'][i] * 180. / np.pi)
    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax.add_artist(e)


# In[44]:


#SZHERR

# measure a spatially varying background on the image
rherbkg = sep.Background(finalrher)
# get a "global" mean and noise of the image background:
print(rherbkg.globalback)
print(rherbkg.globalrms)
# evaluate background as 2-d array, same size as original image
rherbkg_image = rherbkg.back()
# rherbkg_image = np.array(rherbkg) # equivalent to above
# show the background
plt.figure(figsize= (10,10))
plt.imshow(rherbkg_image, interpolation='nearest', cmap='gray', origin='lower')
plt.colorbar();
plt.title("SZHER R Background")

# evaluate the background noise as 2-d array, same size as original image
rherbkg_rms = rherbkg.rms()
# show the background noise
plt.figure(figsize= (10,10))
plt.imshow(rherbkg_rms, interpolation='nearest', cmap='gray', origin='lower')
plt.colorbar();
plt.title("SZHER R Background RMS")

# subtract the background
rherdata_sub = finalrher - rherbkg

rherobj = sep.extract(rherdata_sub, 10, err=rherbkg.globalrms)
# how many objects were detected
print(len(rherobj))
#print(bherobj)

from matplotlib.patches import Ellipse

# plot background-subtracted image
fig, ax = plt.subplots()
m, s = np.mean(rherdata_sub), np.std(rherdata_sub)
im = ax.imshow(rherdata_sub, interpolation='nearest', cmap='gray',
               vmin=m-s, vmax=m+s, origin='lower')

# plot an ellipse for each object
for i in range(len(rherobj)):
    e = Ellipse(xy=(rherobj['x'][i], rherobj['y'][i]),
                width=6*rherobj['a'][i],
                height=6*rherobj['b'][i],
                angle=rherobj['theta'][i] * 180. / np.pi)
    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax.add_artist(e)


# In[45]:


#SZHERV

# measure a spatially varying background on the image
vherbkg = sep.Background(finalvher)
# get a "global" mean and noise of the image background:
print(vherbkg.globalback)
print(vherbkg.globalrms)
# evaluate background as 2-d array, same size as original image
vherbkg_image = vherbkg.back()
# vherbkg_image = np.array(vherbkg) # equivalent to above
# show the background
plt.figure(figsize= (10,10))
plt.imshow(vherbkg_image, interpolation='nearest', cmap='gray', origin='lower')
plt.colorbar();
plt.title("SZHER V Background")

# evaluate the background noise as 2-d array, same size as original image
vherbkg_rms = vherbkg.rms()
# show the background noise
plt.figure(figsize= (10,10))
plt.imshow(vherbkg_rms, interpolation='nearest', cmap='gray', origin='lower')
plt.colorbar();
plt.title("SZHER V Background RMS")

# subtract the background
vherdata_sub = finalvher - vherbkg

vherobj = sep.extract(vherdata_sub, 10, err=vherbkg.globalrms)
# how many objects were detected
print(len(vherobj))
#print(vherobj)

from matplotlib.patches import Ellipse

# plot background-subtracted image
fig, ax = plt.subplots()
m, s = np.mean(vherdata_sub), np.std(vherdata_sub)
im = ax.imshow(vherdata_sub, interpolation='nearest', cmap='gray',
               vmin=m-s, vmax=m+s, origin='lower')

# plot an ellipse for each object
for i in range(len(vherobj)):
    e = Ellipse(xy=(vherobj['x'][i], vherobj['y'][i]),
                width=6*vherobj['a'][i],
                height=6*vherobj['b'][i],
                angle=vherobj['theta'][i] * 180. / np.pi)
    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax.add_artist(e)


# In[46]:


#M29R

# measure a spatially varying background on the image
rm29bkg = sep.Background(finalrm29)
# get a "global" mean and noise of the image background:
print(rm29bkg.globalback)
print(rm29bkg.globalrms)
# evaluate background as 2-d array, same size as original image
rm29bkg_image = rm29bkg.back()
# rm29bkg_image = np.array(rm29bkg) # equivalent to above
# show the background
plt.figure(figsize= (10,10))
plt.imshow(rm29bkg_image, interpolation='nearest', cmap='gray', origin='lower')
plt.colorbar();
plt.title("M29 R Background")

# evaluate the background noise as 2-d array, same size as original image
rm29bkg_rms = rm29bkg.rms()
# show the background noise
plt.figure(figsize= (10,10))
plt.imshow(rm29bkg_rms, interpolation='nearest', cmap='gray', origin='lower')
plt.colorbar();
plt.title("M29 R Background RMS")

# subtract the background
rm29data_sub = finalrm29 - rm29bkg

rm29obj = sep.extract(rm29data_sub, 10, err=rm29bkg.globalrms)
# how many objects were detected
print(len(rm29obj))
#print(rm29obj)

from matplotlib.patches import Ellipse

# plot background-subtracted image
fig, ax = plt.subplots()
m, s = np.mean(rm29data_sub), np.std(rm29data_sub)
im = ax.imshow(rm29data_sub, interpolation='nearest', cmap='gray',
               vmin=m-s, vmax=m+s, origin='lower')

# plot an ellipse for each object
for i in range(len(rm29obj)):
    e = Ellipse(xy=(rm29obj['x'][i], rm29obj['y'][i]),
                width=6*rm29obj['a'][i],
                height=6*rm29obj['b'][i],
                angle=rm29obj['theta'][i] * 180. / np.pi)
    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax.add_artist(e)


# In[47]:


#M15R

# measure a spatially varying background on the image
rm15bkg = sep.Background(finalrm15)
# get a "global" mean and noise of the image background:
print(rm15bkg.globalback)
print(rm15bkg.globalrms)
# evaluate background as 2-d array, same size as original image
rm15bkg_image = rm15bkg.back()
# rm15bkg_image = np.array(rm15bkg) # equivalent to above
# show the background
plt.figure(figsize= (10,10))
plt.imshow(rm15bkg_image, interpolation='nearest', cmap='gray', origin='lower')
plt.colorbar();
plt.title("M15 R Background")

# evaluate the background noise as 2-d array, same size as original image
rm15bkg_rms = rm15bkg.rms()
# show the background noise
plt.figure(figsize= (10,10))
plt.imshow(rm15bkg_rms, interpolation='nearest', cmap='gray', origin='lower')
plt.colorbar();
plt.title("M15 R Background RMS")

# subtract the background
rm15data_sub = finalrm15 - rm15bkg

rm15obj = sep.extract(rm15data_sub, 10, err=rm15bkg.globalrms)
# how many objects were detected
print(len(rm15obj))
#print(rm15obj)

from matplotlib.patches import Ellipse

# plot background-subtracted image
fig, ax = plt.subplots()
m, s = np.mean(rm15data_sub), np.std(rm15data_sub)
im = ax.imshow(rm15data_sub, interpolation='nearest', cmap='gray',
               vmin=m-s, vmax=m+s, origin='lower')

# plot an ellipse for each object
for i in range(len(rm15obj)):
    e = Ellipse(xy=(rm15obj['x'][i], rm15obj['y'][i]),
                width=6*rm15obj['a'][i],
                height=6*rm15obj['b'][i],
                angle=rm15obj['theta'][i] * 180. / np.pi)
    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax.add_artist(e)


# In[ ]:




