from glob import glob
import numpy as np
import os

mspath = sys.argv[3]
impath = sys.argv[4]

total_niter = [1000, 3000, 5000]
mask_proportion = [0.0, 0.5, 0.8, 1.0]

summary = open(impath+"summary.txt", "r")
ra0, dec0 = [x[:-1] for x in summary.readlines()]    

os.system("cp -r {} {}".format(mspath, impath))

mspath = impath + mspath.split('/')[-1]

def this_tclean(newimage, niter, mpath=''):
    tclean(vis=mspath, imagename=newimage,
            niter=niter, weighting='briggs', robust=0.5, imsize = [512,512], cell=['500 arcsec'],
            specmode='mfs', nterms=2, mask=mpath, spw='0:100~920',stokes='IQUV', 
            interactive=False, pblimit=-1, phasecenter='J2000 {}deg {}deg'.format(ra0, dec0),
            gridder='widefield')

this_tclean('no_deconvolution', 0)
exportfits(imagename='no_deconvolution.image', fitsimage='no_deconvolution.fits')
    
masks = sorted(glob(impath + "*mask.txt"))   

for n in total_niter:
    for p in mask_proportion:
        n1 = int(n * p)
        n2 = n - n1
        
        for mpath in masks:
            mname = mpath.split('/')[-1].split('_')[-1]
            m = '.'.join(mname.split('.')[:-1])
            newimage = impath+'n{}_p{:.2f}_mask{}'.format(n,p,m)
            print(newimage.split('/')[-1])
            
            this_tclean(newimage, n1, mpath)

            for s in ['residual', 'model', 'image', 'mask']:
                exportfits(imagename='{}.{}'.format(newimage,s), 
                       fitsimage='{}_just_mask.{}.fits'.format(newimage,s))
            
            os.system("rm -rf {}.mask".format(newimage))
            
            this_tclean(newimage, n2)
            
            for s in ['residual', 'model', 'image']:
                exportfits(imagename='{}.{}'.format(newimage,s), 
                       fitsimage='{}_final.{}.fits'.format(newimage,s))
                
                
os.system("rm -rf {}summary.txt".format(impath))