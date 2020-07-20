from glob import glob
import numpy as np
import os

mspath = sys.argv[3]
impath = sys.argv[4]

total_niter = [1000, 3000, 5000]
mask_proportion = [0, 0.5, 0.8, 1]

summary = open(impath+"info.txt", "r")
ra0, dec0 = [x[:-1] for x in summary.readlines()]    

#tclean(vis=mspath,imagename=impath+'no_deconvolution',
#                niter=1, imsize = [512,512], cell=['500 arcsec'],
#                specmode='mfs', spw='0:100~920', stokes='IQUV', 
#                interactive=False, pblimit=-1, #phasecenter='J2000 {}deg {}deg'.format(ra0, dec0),
#               gridder='widefield')

masks = glob(impath + "*mask.txt")    

for n in total_niter:
    for p in mask_proportion:
        n1 = int(n * p)
        n2 = n - n1
        
        for mpath in masks:
            m = mpath[-14:-8]
            tclean(vis=mspath, imagename=impath+'n={}_p={}_mask={}'.format(n,p,m),
                    niter=n1, weighting='briggs', robust=0.5, imsize = [512,512], cell=['500 arcsec'],
                    specmode='mfs', nterms=2, mask=mpath, spw='0:100~920',stokes='IQUV', 
                    interactive=False, pblimit=-1, #phasecenter='J2000 {}deg {}deg'.format(ra0, dec0),
                    gridder='widefield')

            os.system("rm -rf {}n={}_p={}_mask={}.mask".format(impath,n,p,m))
            
            for s in ['residual', 'model', 'image']:
                exportfits(imagename=impath+'n={}_p={}_mask={}.{}'.format(n,p,m,s), 
                       fitsimage=impath+'n={}_p={}_mask={}_just_mask.{}.fits'.format(n,p,m,s))
            
            tclean(vis=mspath, imagename=impath+'n={}_p={}_mask={}'.format(n,p,m),
                    niter=n2, weighting='briggs', robust=0.5, imsize = [512,512], cell=['500 arcsec'],
                    specmode='mfs', nterms=2, spw='0:100~920',stokes='IQUV', 
                    interactive=False, pblimit=-1, #phasecenter='J2000 {}deg {}deg'.format(ra0, dec0),
                    gridder='widefield', usemask='pb', pbmask=0.0)
           
            
            for s in ['residual', 'model', 'image']:
                exportfits(imagename=impath+'n={}_p={}_mask={}.{}'.format(n,p,m,s), 
                       fitsimage=impath+'n={}_p={}_mask={}_final.{}.fits'.format(n,p,m,s))
                