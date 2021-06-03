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
def plot_image_difference(fitsfile_1, fitsfile_2, reference="1", title=None, save_fig=False):
    stokes=["I",'Q','U','V']
    hdu_1 = fits.open(fitsfile_1)
    image_wcs = WCS(hdu_1[0].header, naxis=[1,2])
    image_1 = hdu_1[0].data
    
    hdu_2 = fits.open(fitsfile_2)
    image_2 = hdu_2[0].data
    
    fig, axes = plt.subplots(figsize=(15,20), dpi=300, nrows=4, ncols=3, sharex=True, sharey=True, 
                             subplot_kw = {'projection' : image_wcs})
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=0.2)
    for i in range(4):
        if reference == "1":
            vmin, vmax = np.min(image_1[i,0,:,:]), np.max(image_1[i,0,:,:])
        elif reference == "2":
            vmin, vmax = np.min(image_2[i,0,:,:]), np.max(image_2[i,0,:,:])
        else:
            raise ValueError("Reference value must be either '1' or '2'.")
        ax = axes[i,0]
        im_1 = ax.imshow(image_1[i,0,:,:], origin='lower', interpolation='nearest', cmap='Blues', aspect='auto', vmin=vmin, vmax=vmax)
        ax.grid(color='white',alpha=0.5)
        ax.set_title("image1", fontsize=16)
        ax.text(0.1,0.9, "{}".format(stokes[i]), bbox=dict(facecolor='white', alpha=0.1), transform=ax.transAxes, fontsize=20, color='k')
        ax.set_ylabel('Dec (J2000)', fontsize=16)
        if i == 3:
            ax.set_xlabel('RA (J2000)', fontsize=16)
        else:
            ax.set_xlabel(' ')
        ax.tick_params(axis='both', which='both', labelsize=14)
        
        ax = axes[i,1]
        im_2 = ax.imshow(image_2[i,0,:,:], origin='lower', interpolation='nearest', cmap='Blues', aspect='auto', vmin=vmin, vmax=vmax)
        ax.grid(color='white',alpha=0.5)
        ax.set_title("image2", fontsize=16)
        if i == 3:
            ax.set_xlabel('RA (J2000)', fontsize=16)
        else:
             ax.set_xlabel(' ')
        ax.set_ylabel(' ')
        ax.tick_params(axis='both', which='both', labelsize=14)
        
        ax = axes[i,2]
        im_3 = ax.imshow(image_1[i,0,:,:]-image_2[i,0,:,:], origin='lower', interpolation='nearest', cmap='Blues', aspect='auto', vmin=vmin, vmax=vmax)
        ax.grid(color='white',alpha=0.5)
        ax.set_title("image1-image2", fontsize=16)
        if i == 3:
            ax.set_xlabel('RA (J2000)', fontsize=16)
        else:
            ax.set_xlabel(' ')
        ax.set_ylabel(' ')
        cax = ax.inset_axes([1.04, 0.1, 0.05, 0.8], transform=ax.transAxes)
        cbar = fig.colorbar(im_3, ax=ax, cax=cax)
        cbar.ax.tick_params(labelsize=10)
        ax.tick_params(axis='both', which='both', labelsize=14)
    if title is not None:
        fig.suptitle(str(title), fontsize=16, y=0.92)
    if save_fig==True:
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
    
from cycler import cycler

def diff_uv_1d(uvd1, uvd2, check_metadata=True, bins=50, title=None):
    """
    Parameters
    ----------
    uvd1, uvd2 : pyuvdata.UVData
        Input UVData objects which contain the visibilities to be differenced
        and any relevant metadata.
    check_metadata : bool, optional
        Whether to check that the metadata for `uvd1` and `uvd2` match.
        See ``utils.check_uvd_pair_metadata`` docstring for details on
        how the metadata are compared. If `check_metadata` is set to
        False, but the metadata don't agree, then the plotter may or
        may not error out, depending on how the metadata disagree.
        Default behavior is to check the metadata.
    bins : int, optional
        Number of bins to use for regridding the u and v arrays.
    """
    # check the metadata unless instructed otherwise
    if check_metadata:
        uvtools.utils.check_uvd_pair_metadata(uvd1, uvd2)

    # convert polarization to index

    # load in relevant metadata
    bl_vecs = uvd1.uvw_array
    freqs = uvd1.freq_array[0]

    # import astropy constants to convert freq to wavelength
    from astropy.constants import c
    wavelengths = c.value / freqs

    # get uvw vectors; shape = (Nfreq, Nblts, 3)
    uvw_vecs = np.array([bl_vecs / wavelength for wavelength in wavelengths])

    # reshape uvw vectors to (Nblts, Nfreq, 3)
    uvw_vecs = np.swapaxes(uvw_vecs, 0, 1)

    # get the u and v arrays, flattened
    uvals, vvals = uvw_vecs[:,:,0].flatten(), uvw_vecs[:,:,1].flatten()
    uv_radius = np.sqrt(uvals**2+vvals**2)
    
    # get the visbilities
    vis1, vis2 = np.zeros((4, *uvals.shape)).astype(np.complex128), np.zeros((4, *uvals.shape)).astype(np.complex128) 
    vis1[0,:], vis1[1,:], vis1[2,:], vis1[3,:] = uvd1.data_array[:,0,:,0].flatten(), uvd1.data_array[:,0,:,1].flatten(), uvd1.data_array[:,0,:,2].flatten(), uvd1.data_array[:,0,:,3].flatten()
    vis2[0,:], vis2[1,:], vis2[2,:], vis2[3,:] = uvd2.data_array[:,0,:,0].flatten(), uvd2.data_array[:,0,:,1].flatten(), uvd2.data_array[:,0,:,2].flatten(), uvd2.data_array[:,0,:,3].flatten()
    
    abs_vis1, abs_vis2 = np.abs(vis1), np.abs(vis2)
    phase_vis1, phase_vis2 = np.angle(vis1), np.angle(vis2)
    # get the regridded radius
    uv_radius_bins = np.linspace(np.min(uv_radius[uv_radius>0]), uv_radius.max(), bins+1)
    uv_radius_bin_width = np.mean(np.diff(uv_radius_bins))
    hist_bins_vis1, hist_bins_vis2 = np.zeros((2,4, bins+1)), np.zeros((2,4,bins+1)) 
    
    for (i,bin_radius) in enumerate(uv_radius_bins):
        for pol in range(4):
            hist_bins_vis1[0,pol,i] += np.sum(abs_vis1[pol,:] [(uv_radius >= bin_radius) & (uv_radius < bin_radius+uv_radius_bin_width)])
            hist_bins_vis1[1,pol,i] += np.sum(phase_vis1[pol,:] [(uv_radius >= bin_radius) & (uv_radius < bin_radius+uv_radius_bin_width)])
            hist_bins_vis2[0,pol,i] += np.sum(abs_vis2[pol,:] [(uv_radius >= bin_radius) & (uv_radius < bin_radius+uv_radius_bin_width)])
            hist_bins_vis2[1,pol,i] += np.sum(phase_vis2[pol,:] [(uv_radius >= bin_radius) & (uv_radius < bin_radius+uv_radius_bin_width)])
    # setup the figure
    fig, axes = plt.subplots(figsize=(14,12), ncols=2, nrows=2)
    
    ax = axes[0,0]
    ax.set_prop_cycle(cycler('color', ['c', 'm']))
    l1 = ax.plot(uv_radius_bins,  hist_bins_vis1[0,:2,:].T, ls="--")
    l2 = ax.plot(uv_radius_bins,  hist_bins_vis2[0,:2,:].T, ls=":")
    ymax = np.max([np.max(hist_bins_vis1[0,:2,:]), np.max(hist_bins_vis2[0,:2,:])])
    ax.set_ylim(ymax/1e8, ymax/3)
    ax.set_title("Magnitudes", fontsize=16)
    ax.set_xlabel("uv radius", fontsize=16)
    ax.tick_params(axis='both', which='both', bottom=False, labelsize=12)
    
    ax = axes[0,1]
    ax.set_prop_cycle(cycler('color', ['r', 'g']))
    l3 = ax.plot(uv_radius_bins,  hist_bins_vis1[0,2:,:].T, ls="--")
    l4 = ax.plot(uv_radius_bins,  hist_bins_vis2[0,2:,:].T, ls=":")
    ymax = np.max([np.max(hist_bins_vis1[0,2:,:]), np.max(hist_bins_vis2[0,2:,:])])
    ax.set_ylim(ymax/1e8, ymax/3)
    ax.set_title("Magnitudes", fontsize=16)
    ax.set_xlabel("uv radius", fontsize=16)
    ax.tick_params(axis='both', which='both', bottom=False, labelsize=12)
    
    fig.legend([*l1,*l3,*l2,*l4], ["uvd1_xx", "uvd1_yy", "uvd1_xy", "uvd1_yx", "uvd2_xx", "uvd2_yy", "uvd2_xy", "uvd2_yx"], bbox_to_anchor=(0.53,1.0), loc='upper center', fontsize=14, ncol=4,)
   
    ax = axes[1,0]
    ax.set_prop_cycle(cycler('color', ['c', 'm']))
    ax.plot(uv_radius_bins,  hist_bins_vis1[1,:2,:].T, ls="--")
    ax.plot(uv_radius_bins,  hist_bins_vis2[1,:2,:].T, ls=":", )
    ax.set_ylim(np.min([np.min(hist_bins_vis1[1,:2,:]), np.min(hist_bins_vis2[1,:2,:])]), np.max([np.max(hist_bins_vis1[1,:2,:]), np.max(hist_bins_vis2[1,:2,:])]))
    ax.set_title("Phases", fontsize=16)
    ax.set_xlabel("uv radius", fontsize=16)
    ax.tick_params(axis='both', which='both', bottom=False, labelsize=12)
    
    ax = axes[1,1]
    ax.set_prop_cycle(cycler('color', ['r', 'g']))
    ax.plot(uv_radius_bins,  hist_bins_vis1[1,2:,:].T, ls="--")
    ax.plot(uv_radius_bins,  hist_bins_vis2[1,2:,:].T, ls=":", )
    ax.set_ylim(np.min([np.min(hist_bins_vis1[1,2:,:]), np.min(hist_bins_vis2[1,2:,:])]), np.max([np.max(hist_bins_vis1[1,2:,:]), np.max(hist_bins_vis2[1,2:,:])]))
    ax.set_title("Phases", fontsize=16)
    ax.set_xlabel("uv radius", fontsize=16)
    ax.tick_params(axis='both', which='both', bottom=False, labelsize=12)
    
    if title is not None:
        fig.suptitle(str(title), y=1.08, fontsize=16)

def plot_uvw_temporal_behavior(uvd, baseline):
    if isinstance(baseline, (int, np.int)):
        baseline_num = np.int64(baseline)
    elif (isinstance(baseline, tuple) and isinstance(baseline[0], (int, np.int))):
        baseline_num = np.int64(2048*(baseline[0]+1)+(baseline[1]+1)+2**16)
    else:
        raise TypeError("baseline must be an int or a tuple of int.")
    uvw_baseline = uvd.uvw_array[np.array(np.where(uvd.baseline_array == baseline_num)).flatten(),:]
    time_array = np.unique(uvd.time_array)
    Julian_Date = np.floor(time_array[0])
    time_array -= Julian_Date
    fig, ax = plt.subplots(figsize=(6,4))
    ax.plot(time_array, uvw_baseline[:,0], label="u", ls="-")
    ax.plot(time_array, uvw_baseline[:,1], label="v", ls="--")
    ax.plot(time_array, uvw_baseline[:,2], label="w", ls=":")
    ax.set_xlabel("time JD{}".format(Julian_Date), fontsize=16)
    ax.set_ylabel("uvw", fontsize=16)
    ax.legend(loc='best', fontsize=14)
    ax.tick_params(axis='both', which='both', labelsize=16)
    ax.set_title("baseline {}".format(uvd.baseline_to_antnums(baseline_num)),fontsize=16)