from glob import glob
import sys

impath = sys.argv[3]
images = [x[:-9] for x in glob(impath + "*.info.txt")]

for file in images:
    summary = open(file + ".info.txt", "r")
    fitsfile, ra0, dec0 = [x[:-1] for x in summary.readlines()]
    importuvfits(fitsfile, file+".ms")
    
    tclean(vis=file+'.ms', imagename=file+'.no_deconvolution', niter=0, weighting='briggs', robust=0, imsize = [512,512], pbcor=False, cell=['500 arcsec'], specmode='mfs', nterms=1, spw='0:100~920', stokes='IQUV', interactive=False, pblimit=-1)
    
    tclean(vis=file+'.ms', imagename=file+'.deconvolved', niter=5000, weighting='briggs', robust=0, imsize = [512,512], pbcor=False, cell=['500 arcsec'], specmode='mfs', nterms=1, spw='0:100~920', stokes='IQUV', mask='%s.masks.txt', interactive=False, cycleniter=1000, threshold='1Jy/beam', pblimit=-1)

    exportfits(file+".no_deconvolution.image", file+".no_deconvolution.fits")
    exportfits(file+".deconvolved.image", file+".deconvolved.fits")
