import numpy as np, matplotlib.pyplot as plt
import argparse
import glob
import aplpy

from astropy.io import fits
from astropy.wcs import WCS
from matplotlib import gridspec
from astropy.utils.data import get_pkg_data_filename

"""
    The Purpose of this script it to find the proper scaling factor S and convert unity to Jansky. This script also makes Stokes images with two type of corrd axies [galactic long. {Deg}, galactic lat. {Deg}] and [RA {Hr}, Dec {Deg}].
"""


def S_sgrA(nu): # Flux density of GC
    # Make Amplitude Scale Flux Desnity of the Non-Linear Galactic Center
    # SSgrA∗ (ν) ≈ 3709 Jy × (ν/408 MHz)−0.5 ---> Saul's Paper [Eq.11]
    
    flux_density= 3709 # [Jansky]
    nu_0 = 408 # [MHz]
    alpha = -05
    
    S_of_nu = flux_density*np.power(nu/nu_0,-0.5)
    hdu = fits.PrimaryHDU(S_of_nu)
    hdu.writeto("Scale_Spectrum_of_GC.fits",overwrite=True) #save
    
    return S_of_nu

"""
    # This makes a spectrum and plots it.
    
    nimage = 150
    freq = np.linspace(115,188,nimage)
    S = S_sgrA(freq)
    
    for ind in range(nimage):
    
    plt.plot(chan[ind],S[ind].max(),'.')
    plt.xlabel("Channel [Arb]")
    #plt.plot("Freq [MHz]")
    plt.ylabel("Flux Density [Arb]")
    
    plt.title('Flux Density Spectrum Scale Value')
    plt.savefig('{}_ScaleValue_FluxDensitySpectrum.png'.format(filename))
    plt.show()
    plt.close()
    """

files_from_fits = glob.glob('zen.2457755.89889.conj.HH.uv.TrueVis.TrueVis_NoDeconvolution_Saul_paper_Fig4.image.fits')[0]

my_sim = fits.open(files_from_fits)
freq_value = my_sim[0].header['CRVAL3']*1e-6 # In units on MHz
copy_fits = my_sim[0].copy()

S = S_sgrA(np.abs(freq_value)) #insert the freq given in GHz. Convert to MHz and take the absolute value
copy_fits.data *= S

# Write it out to a new fits file
# You probably don't need to do this step, just fine the factor array and then mult it by the data but it might be a good idea to have the numbers saved somewhere.

#hdu1 = fits.PrimaryHDU(S)
#hdu1.writeto("factor.fits",overwrite=True)

copy_fits.writeto("scaled_{}".format(files_from_fits),overwrite=True)

# New STOKES Visibility Plots
# Mult my fits data by this factor and create new product fits file

#----------------------------------------------------------------------------------

#==================================================================================
# No WCS header added
#factor = fits.open("factor.fits")[0].data
new_my_sim = S*my_sim_data
wcs = WCS(my_sim[0].header,naxis=2)

iarr = [0,0,1,1]
jarr = [0,1,0,1]
nrow,ncol = 2,2
stoke = ['I','Q','U','V']
name = "scaled_{}".format(fitsfile)


for fitsfile in files_from_fits:
    
    f,axarr = plt.subplots(nrow,ncol,sharex=True,sharey=True,figsize=(15,8))
    f.suptitle("Scaled STOKES Visibility Plots: {}".format(fitsfile), fontsize=14)
    gs = gridspec.GridSpec(nrow,ncol)
    
    for pol in range(4):
        ax=plt.subplot(gs[iarr[pol],jarr[pol]],projection=wcs)
        #ax[iarr[pol],jarr[pol]].subplot(projection=wcs)
        
        if pol == 0:
            vmax=3.e4
            vmin=0
            cmap="viridis"
            cax=ax.imshow(new_my_sim[pol,0,:,:], vmax=vmax, vmin=vmin, cmap=cmap, origin='lower')
            f.colorbar(cax,label='Stokes '+stoke[pol])
            
            ax.grid(color='white', ls='solid',which='major')
            ax.set_xlabel('Galactic Longitude [Deg]')
            ax.set_ylabel('Galactic Latitude [Deg]')
        
        if pol == 1:
            vmax=3.e2
            vmin=-3.e2
            cmap="RdYlGn"#"PRGn"
            cax=ax.imshow(new_my_sim[pol,0,:,:], vmax=vmax, vmin=vmin, cmap=cmap, origin='lower')
            f.colorbar(cax,label='Stokes '+stoke[pol])
            
            ax.grid(color='k', ls='solid',which='major')
            ax.set_xlabel('Galactic Longitude [Deg]')
            ax.set_ylabel('Galactic Latitude [Deg]')
        if pol == 2:
            vmax=3.e2
            vmin=-3.e2
            cmap="RdYlGn"#"PRGn"
            cax=ax.imshow(new_my_sim[pol,0,:,:], vmax=vmax, vmin=vmin, cmap=cmap, origin='lower')
            f.colorbar(cax,label='Stokes '+stoke[pol])
            
            ax.grid(color='k', ls='solid',which='major')
            ax.set_xlabel('Galactic Longitude [Deg]')
            ax.set_ylabel('Galactic Latitude [Deg]')
        if pol == 3:
            vmax=3.e2
            vmin=-3.e2
            cmap="RdYlGn"#"PRGn"
            cax=ax.imshow(new_my_sim[pol,0,:,:], vmax=vmax, vmin=vmin, cmap=cmap, origin='lower')
            f.colorbar(cax,label='Stokes '+stoke[pol])
            
            ax.grid(color='k', ls='solid',which='major')
            ax.set_xlabel('Galactic Longitude [Deg]')
            ax.set_ylabel('Galactic Latitude [Deg]')

    f.savefig("type1_{}_scaled_STOKES.png".format(fitsfile))
    #plt.show()
    plt.close()

"""
# I think this is the orginal command line stuff

factor = fits.open("factors.fits")[0].data
new_my_sim = factor*my_sim_data
wcs = WCS(my_sim[0].header,naxis=2)

plt.subplot(projection=wcs)
plt.imshow(new_my_sim[0,0,:,:], vmin=0, vmax=3.e4, origin='lower')
plt.colorbar()
plt.grid(color='white', ls='solid')
plt.xlabel('Galactic Longitude [Deg]')
plt.ylabel('Galactic Latitude [Deg]')
plt.show()
"""
#==================================================================================

#----------------------------------------------------------------------------------

# WCS header added to new fits file. Now you can use aply to make the plots.
files_from_fits = glob.glob(name)

for fitsfile in files_from_fits:
    
    f = plt.figure(figsize=(15,8))
    for pol in np.arange(4):
        fig = aplpy.FITSFigure(fitsfile,dimensions=[0,1],slices=[0,pol],figure=f,subplot=(2,2,pol+1))
        if pol == 0:
            vmax=3.e4
            vmin=0
            cmap="viridis"
        if pol == 1:
            vmax=3.e2
            vmin=-3.e2
            cmap="RdYlGn"#"PRGn"
        if pol == 2:
            vmax=3.e2
            vmin=-3.e2
            cmap="RdYlGn"#"PRGn"
        if pol == 3:
            vmax=3.e2
            vmin=-3.e2
            cmap="RdYlGn"#"PRGn"
        
        fig.show_colorscale(cmap=cmap,vmax=vmax,vmin=vmin)#,stretch='arcsczdxcinh')
        fig.add_grid()
        fig.grid.set_color(color='black')#, linestyle='solid')
        fig.grid.set_xspacing(15)
        fig.grid.set_yspacing(15)
        fig.grid.show()
        fig.axis_labels.set_font(size='small')
        fig.tick_labels.set_font(size='small')
        fig.tick_labels.set_xformat('hh')
        fig.tick_labels.set_yformat('dd')
        fig.add_colorbar()
        fig.colorbar.set_font(size='small')

    plt.suptitle('{} STOKE Visibilities'.format(fitsfile))
    fig.savefig('type2_{}_scaled_STOKES.png'.format(fitsfile))
    plt.close()
#----------------------------------------------------------------------------------

#==================================================================================



