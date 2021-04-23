import numpy as np
import os

uvfitspaths = sys.argv[3].split(',')
impath = sys.argv[4]

summary = open(impath+"summary.txt", "r")
ra0, dec0 = [x[:-1] for x in summary.readlines()]    
summary.close()

def this_tclean(vis, newimage, niter=2000, mask='', wprojplanes=1, facets=1, side=512, nterms=2, robust=0.5,
               deconvolver='clarkstokes', scales=[0], spw='0:0~100'):
    
    # add spw='0:100~920', for real data
    tclean(vis=vis, imagename=newimage, niter=niter, weighting='briggs', 
            robust=robust, imsize = [side,side], cell=['500 arcsec'], scales=scales,
            specmode='mfs', nterms=nterms, mask=mask, stokes='IQUV', spw=spw, 
            interactive=False, pblimit=-1, phasecenter='J2000 {}deg {}deg'.format(ra0, dec0),
            gridder='widefield', wprojplanes=wprojplanes, facets=facets)
    
    for s in ['image', 'model', 'residual', 'psf']:
        exportfits(imagename='{}.{}'.format(newimage, s), 
                 fitsimage='{}.{}.fits'.format(newimage, s))
   
for path in uvfitspaths:
    bpath = '/'.join(path.split('/')[:-1]) + '/' + '.'.join(path.split('/')[-1].split('.')[:-1])
    mspath = bpath + '.ms'
    importuvfits(fitsfile=path, vis=mspath)
    this_tclean(vis=mspath, newimage=bpath+'_no_deconvolution', niter=0)
    this_tclean(vis=mspath, newimage=bpath+'_deconvolved')



'''
for f in [1,5,10,20,100]:
    for t in [10]:#[1,2,4,7,10]:
        this_tclean('f_0:{:.0f}_t_{:.0f}mins'.format(f,t), 
                    spw='0:0~'+str(f), #timerange='11:40:00~11:{:.0f}:00'.format(t+40),
                    mask = "circle[[{0}deg, {1}deg], 10deg]".format(
                                        str(ra0), str(dec0))) 
'''