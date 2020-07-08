#!/usr/bin/env python
# coding: utf-8

# In[38]:


# make a bunch of different tclean options to run and test each of the residuals...
from glob import glob
import numpy as np


file = '/lustre/aoc/projects/hera/nkern/CHAMP_Bootcamp/Lesson10_HERADataPartII/data/zen.2458116.24482.xx.HH.uvOCR.ms/'
beam = '/lustre/aoc/projects/hera/gtucker/repositories/team_polarization/readable_herapb.tab/'
savepath='/lustre/aoc/projects/hera/gtucker/repositories/team_polarization/casa_fits/correctpb/'
path='/lustre/aoc/projects/hera/gtucker/repositories/team_polarization'
niter = [1000,3000,5000]
logical = [True,False]
nterms = [1,2]
g = 'widefield'

#add multiscale to this

data = np.genfromtxt('/lustre/aoc/projects/hera/gtucker/notebooks/cords.txt') #for phase center
ra0,dec0 = data[0],data[1]
robust = [-1,-0.5,0]
for i in niter:
        
    for r in robust:

        for t in nterms:
          

            tclean(vis=file,imagename=savepath+'{}_{}_{}_{}'.format(i,r,g,t),
                niter=i, weighting='briggs', robust=r,imsize = [512,512], cell=['500 arcsec'],
                specmode='mfs', nterms=t,mask='masks2.txt', spw='0:100~920',stokes='IQUV', 
                interactive=False, pblimit=-1, phasecenter='J2000 %sdeg %sdeg' % (ra0, dec0),gridder=g,vptable=beam)


            exportfits(imagename=savepath+'{}_{}_{}_{}.residual'.format(i,r,g,t),fitsimage=savepath+'{}_{}_{}_{}.fits'.format(i,r,g,t),overwrite=True)
