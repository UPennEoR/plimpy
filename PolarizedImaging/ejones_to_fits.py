def ejones_to_fits(ejones):
    """
    Make a fits file from a complex E Jones
    
    Arguments:
    ejones = in form of (nx,ny,4,1) : according to casa, (testing for this...)
    
    Returns two fits files of real and imaginary components of ejones in the working directory
    
    """
    
   
    
    #Real
    hdu = fits.PrimaryHDU()  #make empty fits header, fill it with data
    hdu.data = ejones.real
    
    #in future when casa shape is known, better to savve it at real_ejones.fits for usability
    hdu.writeto('R_{}_{}_{}_{}_jones.fits'.format(ejones.shape[0],ejones.shape[1],ejones.shape[2],ejones.shape[3]))
    
    #Imag
    hdu = fits.PrimaryHDU()
    hdu.data = ejones.imag
    hdu.writeto('I_{}_{}_{}_{}_jones.fits'.format(ejones.shape[0],ejones.shape[1],ejones.shape[2],ejones.shape[3]))
    
    
    return "Finished! Shape of input was {}_{}_{}_{}".format(ejones.shape[0],ejones.shape[1],ejones.shape[2],ejones.shape[3])

