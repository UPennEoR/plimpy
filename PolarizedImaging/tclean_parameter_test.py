import numpy as np
import os

uvfitspath = sys.argv[3]
impath = sys.argv[4]

summary = open(impath+"summary.txt", "r")
ra0, dec0 = [x[:-1] for x in summary.readlines()]    
summary.close()

mspath = impath + '.'.join(uvfitspath.split('/')[-1].split('.')[:-1] + ['ms'])

importuvfits(fitsfile=uvfitspath, vis=mspath)

def this_tclean(newimage, niter=2000, mask='', wprojplanes=1, facets=1, side=512, nterms=2, robust=0.5,
               deconvolver='hogbom'):
    
    tclean(vis=mspath, imagename=impath+newimage, niter=niter, weighting='briggs', 
            robust=robust, imsize = [side,side], cell=['500 arcsec'],
            specmode='mfs', nterms=nterms, mask=mask, spw='0:100~920', stokes='IQUV', 
            interactive=False, pblimit=-1, phasecenter='J2000 {}deg {}deg'.format(ra0, dec0),
            gridder='widefield', wprojplanes=wprojplanes, facets=facets)
    
    for s in ['image', 'model', 'residual', 'psf']:
        exportfits(imagename='{}.{}'.format(impath+newimage,s), 
                 fitsimage='{}.{}.fits'.format(impath+newimage,s))
      
    return

this_tclean('no_deconvolution', niter=0)

this_tclean('multiscale', deconvolver='multiscale')
this_tclean('mtmfs', deconvolver='mtmfs')

for p in [1,3]:
    for f in [1,3]:
        this_tclean('planes_{}_facets_{}'.format(p, f),
                    wprojplanes=p, facets=f)

this_tclean('planes_5', wprojplanes=5)
this_tclean('facets_5', facets=5)
this_tclean('1024x1024', side=1024)
this_tclean('nterms_4', nterms=4)

for r in [0,0.25,0.75,1,-2,-1,2]:
    this_tclean('r_' + str(r), robust=r) 
