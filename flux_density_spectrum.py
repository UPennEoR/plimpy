#!/bin/sh

#  flux_density_spectrum.py
#  
#
#  Created by Tashalee Billings on 7/25/18.
#  

"""
   The purpose of this document is to look at how the flux density of the visibility data behaves at different channels. These values have not been scaled but a good wat to check if scaling is needed is to look at this plot. This is to be ran in Python.
"""

import numpy as np, matplotlib.pyplot as plt
import glob
import argparse

from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename

parser = argparse.ArgumentParser(description='Creates Unscaled Flux Density Spectrum. Run with Python as: python flux_desntiy_spectrum.py --<args>=')
parser.add_argument("--fits_image", type=str, default=None, help="Fits image file name. eg '01_zen.2457755.89889.conj.HH.uv.TrueVis.image.fits' ")
parser.add_argument("--ra_pixel", type=int, default=None, help="Open CASA image file, pick a channel and identify the RA pixel number associated with the largest Flux Density of the GC.")
parser.add_argument("--dec_pixel", type=int, default=None, help="Open CASA image file, pick a channel and identify the DEC pixel number associated with the largest Flux Density of the GC.")

args = parser.parse_args()

if __name__ == __main__: #This only runs is we run directly but not if you try to import it.

    filename = args.fits_image
    data=fits.open(filename)[0].data # (4, 150, 512, 512)=(stokes,nimage,??,??)

    chan = range(150,900,5) # 1D array of length 900-150/5 = 150
    nimage = np.shape(data)[1]
    freq = np.linspace(115,188,nimage) #(start freq, end freq, num of points between)
    ra_pixel,dec_pixel = args.ra_pixel,args.dec_pixel #268,261 # Pixel number for GC

    for ind in range(nimage):
        
        plt.plot(chan[ind],data[0,ind,ra_pixel,dec_pixel],'.')
        #plt.plot(chan[ind],data[0,ind,:,:].max(),'.') # ploting Stokes I
        #plt.plot(chan[ind],data[0,ind,ra_pixel,dec_pixel]-data[0,ind,:,:].max(),'.')
        plt.xlabel("Channel [Arb]")
        plt.ylabel("Flux Density [Arb]")

    plt.title('{} Unscaled Flux Density Spectrum'.format(filename))
    plt.savefig('{}_Unscaled_FluxDensitySpectrum.png'.format(filename))
    #plt.show()
    plt.close()

