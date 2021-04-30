from astropy.wcs import WCS
import healpy as hp
from pyuvdata import UVData
from astropy.io import fits
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import numpy as np
import glob as glob
import casatasks
from casatasks import importuvfits, exportfits, tclean,importfits
import re
import cv2

## A wrapper for converting uvh5(uvfits) files to fits files and 
## producing images (casa) and videos 
def uvh5_2_fits(filename, outputpath, 
                good_ants=list(set([9,10,20,22,31,43,53,64,65,72,80,81,88,89,96,97,104,105,112]) - set([22,43,80,81])), ra=0, dec=0, niter=0, \
           weighting='briggs', robust=0, imsize=[512,512], pbcor=False, cell=['500 arcsec'], \
           specmode='mfs', nterms=1, stokes='IQUV', interactive=False, pblimit=-1, **kwargs):
    uvd = UVData()
    uvd.read(filename, antenna_nums=good_ants)
    if uvd.phase_type != 'phased':
        uvd.phase_to_time(np.unique(uvd.time_array)[uvd.Ntimes//2])
    dec, ra = uvd.phase_center_dec_degrees, uvd.phase_center_ra_degrees
    if filename.endswith(".uvh5"):
        uvfitsfile = outputpath + filename.split('/')[-1][:-5] + ".uvfits"
        uvd.write_uvfits(uvfitsfile,spoof_nonessential=True,force_phase=True)
        msfile = uvfitsfile[:-7] + ".ms"
        imagefilestem = msfile[:-3] + ".no_deconvolution"
    elif filename.endswith(".uvfits"):
        uvfitsfile = filename
        msfile = outputpath + uvfitsfile.split('/')[-1][:-7] + ".ms"
        imagefilestem = msfile[:-3] + ".no_deconvolution"
    else:
        print("The filename must end with '.uvh5' or '.uvfits'!")
    
    importuvfits(uvfitsfile, msfile)
    tclean(vis=msfile, imagename=imagefilestem, niter=niter, \
           weighting=weighting, robust=robust, imsize=imsize, pbcor=pbcor, cell=cell, \
           specmode=specmode, nterms=nterms, stokes=stokes, interactive=interactive, pblimit=pblimit, 
           phasecenter='J2000 %sdeg %sdeg' % (ra, dec), **kwargs)
    exportfits(imagename=imagefilestem + ".image", fitsimage=msfile[:-3]+'.fits', overwrite=True)
    
def fits_2_image(fitsfile):
    stokes=["I",'Q','U','V']
    hdu = fits.open(fitsfile)
    image_wcs = WCS(hdu[0].header, naxis=[1,2])
    image = hdu[0].data
    fig, axes = plt.subplots(figsize=(18,16), dpi=300, nrows=2, ncols=2, sharex=True, sharey=True, 
                             subplot_kw = {'projection' : image_wcs})
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.2)
    for i in range(4):
        ax = axes[i//2,i%2]
        im = ax.imshow(image[i,0,:,:], origin='lower', aspect='auto')
        ax.grid(color='white',alpha=0.75)
        cbar = fig.colorbar(im, ax=ax)
        cbar.ax.tick_params(labelsize=16)
        ax.text(0.1,0.9, "{}".format(stokes[i]), bbox=dict(facecolor='red', alpha=0.1), transform=ax.transAxes, fontsize=20, color='w')
        ax.tick_params(labelsize=16)
        # There's unquestionably a better way to get matplotlib to share axes
        if (i == 2 or i == 3):
            ax.set_xlabel('RA (J2000)', fontsize=16)
        else:
            ax.set_xlabel(' ')
        if (i == 0 or i == 2):
            ax.set_ylabel('Dec (J2000)', fontsize=16)
        else:
            ax.set_ylabel(' ')
        ax.tick_params(axis='both', which='both', labelsize=16)
    fig.savefig(fitsfile[:-5]+".png", dpi=300, bbox_inches='tight')
    plt.show()
    
def image_2_video(imagefiles, output_filename, output_format=".avi", fps=4.0):
    frame = cv2.imread(imagefiles[0])
    height, width, channels = frame.shape
    if output_format == ".avi":
        fourcc = cv2.VideoWriter_fourcc('M','J','P','G') # Be sure to use lower case
    if output_format == ".mp4":
        fourcc = cv2.VideoWriter_fourcc(*'DIVX')
        fps = 20
    out = cv2.VideoWriter(output_filename+output_format, fourcc, fps, (width, height))
    for imagefile in imagefiles:
        frame = cv2.imread(imagefile)
        out.write(frame) # Write out frame to video
    out.release()