from glob import glob
import numpy as np

file = '/lustre/aoc/projects/hera/nkern/CHAMP_Bootcamp/Lesson10_HERADataPartII/data/zen.2458116.24482.xx.HH.uvOCR.ms/'
mspath = sys.argv[3]
impath = sys.argv[4]

total_niter = [1000,3000,5000]
mask_proportion = [0,0.5,0.8,1]

summary = open(file+".info.txt", "r")
ra0, dec0 = [x[:-1] for x in summary.readlines()]    

tclean(vis=mspath,imagename=savepath+'no_deconvolution',
                niter=0, imsize = [512,512], cell=['500 arcsec'],
                specmode='mfs', spw='0:100~920', stokes='IQUV', 
                interactive=False, pblimit=-1, phasecenter='J2000 {}deg {}deg'.format(ra0, dec0),
                gridder='widefield')

masks = [x[:-4] for x in glob(impath + "*mask.txt")]

for n in total_niter:
    for p in mask_proportion:
        n1 = n * p
        n2 = n - n1
        
        for m in masks:
            tclean(vis=mspath,imagename=impath+'n={}_p={}_{}_just_mask'.format(n,p,m),
                    niter=n1, weighting='briggs', robust=-1, imsize = [512,512], cell=['500 arcsec'],
                    specmode='mfs', nterms=2, mask=m+'.txt', spw='0:100~920',stokes='IQUV', 
                    interactive=False, pblimit=-1, phasecenter='J2000 {}deg {}deg'.format(ra0, dec0),
                    gridder='widefield')

            tclean(vis=impath+'n={}_p={}_{}_just_mask'.format(n,p,m),
                    imagename=impath+'n={}_p={}_{}_final'.format(n,p,m),
                    niter=n2, weighting='briggs', robust=-1, imsize = [512,512], cell=['500 arcsec'],
                    specmode='mfs', nterms=2, spw='0:100~920',stokes='IQUV', 
                    interactive=False, pblimit=-1, phasecenter='J2000 {}deg {}deg'.format(ra0, dec0),
                    gridder='widefield')
           
            exportfits(imagename=impath+'n={}_p={}_{}_final.residual'.format(n,p,m), 
                       fitsimage=impath+'n={}_p={}_{}_final.residual.fits'.format(n,p,m))

                         
            exportfits(imagename=impath+'n={}_p={}_{}_final.image'.format(n,p,m), 
                       fitsimage=impath+'n={}_p={}_{}_final.image.fita'.format(n,p,m))
            