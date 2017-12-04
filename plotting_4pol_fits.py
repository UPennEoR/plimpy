from astropy.io import fits
import numpy as np, matplotlib.pyplot as plt
import sys
from astropy import wcs

hdulist = fits.open(sys.argv[0])
w = wcs.WCS(hdulist[0].header)

# here: read in CRVAL and CRPIX - 1 (fortran indexing) and CDELT to get ra_{min,max} and dec_{min,max}

f,axarr = plt.subplots(2,2,sharex=True,sharey=True)
for i in range
    if i == 0:
        vmax = 0.6
        vmin = 0.
        cmap = 'viridis'
    else:
        vmax = 0.008
        vmin = -0.008
        cmap = 'RdYlGn'
    ax = axarr.ravel()[i]
    # extent = (ramin,ramax,decmin,decmax)
    im = ax.imshow(np.arcsinh(test.data[i,0,:,:])[::-1],cmap=cmap,vmax=vmax,vmin=vmin) #,extent=extent)
    f.colorbar(im,ax=ax)
# somehow do axis labels (dones't  reallly matter which way round)
# TODO: get rid of white space around images??
f.tight_layout()
plt.show()

