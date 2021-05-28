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
def overlay_dirty_images(fitsfile_1, fitsfile_2, title=None, save_fig=False):
    stokes=["I",'Q','U','V']
    hdu_1 = fits.open(fitsfile_1)
    image_wcs = WCS(hdu_1[0].header, naxis=[1,2])
    image_1 = hdu_1[0].data
    
    hdu_2 = fits.open(fitsfile_2)
    image_2 = hdu_2[0].data
    
    fig, axes = plt.subplots(figsize=(15,10), dpi=300, nrows=2, ncols=2, sharex=True, sharey=True, 
                             subplot_kw = {'projection' : image_wcs})
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.2)
    for i in range(4):
        ax = axes[i//2,i%2]
        vmin, vmax = np.min(image_1[i,0,:,:]), np.max(image_1[i,0,:,:])
        im_1 = ax.imshow(image_1[i,0,:,:], origin='lower', interpolation='nearest', cmap='Reds', aspect='auto', vmin=vmin, vmax=vmax)
        im_2 = ax.imshow(image_2[i,0,:,:], origin='lower', interpolation='nearest', cmap='Blues', alpha=0.8, aspect='auto', vmin=vmin, vmax=vmax)            
        ax.grid(color='white',alpha=0.5)
        cbar_2 = fig.colorbar(im_2, ax=ax, pad=0.08)
        cbar_2.ax.tick_params(labelsize=10)
        cbar_1 = fig.colorbar(im_1, ax=ax, pad=0.08)
        cbar_1.ax.tick_params(labelsize=10)

        ax.text(0.1,0.9, "{}".format(stokes[i]), bbox=dict(facecolor='white', alpha=0.1), transform=ax.transAxes, fontsize=20, color='k')
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
    if title is not None:
        fig.suptitle(str(title), fontsize=16, y=0.95)
    if save_fig:
        fig.savefig(fitsfile[:-5]+".png", dpi=300, bbox_inches='tight')
    plt.show()
    
def uvh5_2_fits(filename, outputpath, ra=0, dec=0, niter=0, \
           weighting='briggs', robust=0, imsize=[512,512], pbcor=False, cell=['500 arcsec'], \
           specmode='mfs', nterms=1, stokes='IQUV', interactive=False, pblimit=-1, **kwargs):
    uvd = UVData()
    uvd.read(filename)
    if uvd.phase_type != 'phased':
        uvd.phase(ra, dec)
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
           phasecenter='J2000 %srad %srad' % (ra, dec), **kwargs)
    exportfits(imagename=imagefilestem + ".image", fitsimage=msfile[:-3]+'.fits', overwrite=True)
    
def fits_2_image(fitsfile, title=None, save_fig=False):
    stokes=["I",'Q','U','V']
    hdu = fits.open(fitsfile)
    image_wcs = WCS(hdu[0].header, naxis=[1,2])
    image = hdu[0].data
    fig, axes = plt.subplots(figsize=(14,12), dpi=300, nrows=2, ncols=2, sharex=True, sharey=True, 
                             subplot_kw = {'projection' : image_wcs})
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.2)
    for i in range(4):
        ax = axes[i//2,i%2]
        im = ax.imshow(image[i,0,:,:], origin='lower', aspect='auto', interpolation='nearest',)
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
    if title is not None:
        fig.suptitle(str(title), fontsize=16, y=0.95)
    if save_fig:
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