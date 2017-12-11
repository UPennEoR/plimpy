from astropy.io import fits
import numpy as np, matplotlib.pyplot as plt
import sys, glob
from matplotlib import gridspec

norm = True

hdulist = fits.open('4pol2457548.45923.fits')
hdr = hdulist[0].header

racen,deccen = hdr['CRVAL1'],hdr['CRVAL2']
racenpix,deccenpix = int(hdr['CRPIX1'])-1, int(hdr['CRPIX2'])-1
radelt,decdelt = hdr['CDELT1'],hdr['CDELT2']

ramin = racen + (radelt*racenpix)
ramax = racen - (radelt*racenpix)
decmin = deccen - (decdelt*deccenpix)
decmax = deccen + (decdelt*deccenpix)
extent = (ramin,ramax,decmin,decmax)

nrow,ncol = 2,2

top=1.-0.5/(nrow+1)
bottom=0.5/(nrow+1)
left=0.5/(ncol+1)
right=1-0.5/(ncol+1)

f,axarr = plt.subplots(nrow,ncol,sharex=True,sharey=True)
gs = gridspec.GridSpec(nrow,ncol, width_ratios=[1, 1], wspace=0.01, hspace=0.1, top=top,
                       bottom=bottom, left=left, right=right) 
iarr = [0,0,1,1]
jarr = [0,1,0,1]
for i in range(4):
    ax = plt.subplot(gs[iarr[i],jarr[i]])
    if not norm:
        if i == 0:
            vmax = 0.6
            vmin = -0.025
            cmap = 'viridis'
        else:
            vmax = 0.008
            vmin = -0.008
            cmap = 'RdYlGn'
    else:
        if i==0:
            vmax = 1.
            vmin = -0.04
            cmap = 'viridis'
        else:
            vmax = 0.01
            vmin = -0.01
            cmap = 'RdYlGn'
    if i<=1:
        ax.set_xticklabels([])
    else:
        ax.set_xlabel('R.A. [deg.]',size=12)
        
    if i in [1,3]:
        #ax.set_yticklabels([])
        print('yay')
    else:
        ax.set_ylabel('Dec. [deg.]',size=12)
    if norm:
        N = 0.6
    else:
        N = 1.
    img = hdulist[0].data[i,0,:,:][::-1]/N #np.arcsinh(hdulist[0].data[i,0,:,:])[::-1]
    im = ax.imshow(img,cmap=cmap,vmax=vmax,vmin=vmin,extent=extent)
    f.colorbar(im,ax=ax)
#plt.savefig('./image_7548.pdf')
plt.show()
