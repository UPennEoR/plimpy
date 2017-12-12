#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    Created on Thu Nov  2 17:18:23 2017
    
    @author: tashalee
    """

import numpy as np
from astropy.io import fits #for python only
import matplotlib.pyplot as plt
import glob
import os

from fitsconverter import * #for python only

#Convert from .uvfits to .uvfitsMS
#*******************************************************************
uvfits = glob.glob('*.uvfits')
for uvfit in uvfits:
    msfile=uvfit.strip('uvfits') + 'MS'
    importuvfits(vis=msfile,fitsfile=uvfit)
#*******************************************************************


def reorrder(msname):
    import casac
    ms=casac.casac.table()
    ms.open(msname,nomodify=False)
    a1 , a2 , data = [ms.getcol(x) for x in [ "ANTENNA1" , "ANTENNA2" , "DATA" ]]
    m = a1 > a2
    data[:, :, m]= data[:, :, m].conj()
    x = a2[m]
    a2[m] = a1[m]
    a1[m] = x
    ms.putcol("ANTENNA1", a1)
    ms.putcol("ANTENNA2", a2)
    ms.putcol("DATA", data)
    ms.flush()
    ms.close()

def flag(msname): #You might have to update this
    flagdata(msname, flagbackup=True, mode='manual',antenna="22" )
    flagdata(msname, flagbackup=True, mode='manual',antenna="43" )
    flagdata(msname, flagbackup=True, mode='manual',antenna="80" )
    flagdata(msname, flagbackup=True, mode='manual',antenna="81" )

def gen_image(msname,imagename): # Makes image
    clean(msname,imagename=imagename, niter=500, weighting='briggs', robust=0, imsize=[512,512], cell=['500 arcsec'], mode='mfs', nterms=1, spw='0:150~900', stokes='IQUV')

def mkinitmodel(msname,ext): #Model you give casa
    cl.addcomponent(flux=1.0, fluxunit='Jy', shape='point', dir='J2000 17h45m40.0409s -29d0m28.118s')
    cl.rename('GC'+ext+'.cl')
    cl.close()
    ft(msname, complist='GC'+ext+'.cl', usescratch=True)

def clear_cal(msname):
    clearcal(msname)

def phscal(msname):
    import os,sys
    kc = os.path.basename(msname) + "K.cal"
    bc = os.path.basename(msname) + "B.cal"
    gaincal(msname,caltable=kc, gaintype='K', solint='inf', refant='10')
    applycal(msname,gaintable=[kc])
    bandpass(msname,caltable=bc, solint='inf', combine='scan', refant='10')
    applycal(msname,gaintable=[bc])


#%% CREATE IMAGE
MSfilelist = glob.glob('*.MS')
name=MSfilelist
imgnam1='clean1_'
imgnam2='clean2_'

nwms=[]
for ms in MSfilelist:
    msf=ms.strip('zen.'+'HH.uvcRU.MS')
    nwms.append(msf+'.R')

i=0
for i in np.arange(len(MSfilelist)):
    reorrder(msname=name[i])# Reorder antenna 1 & 2 in each correlation so 2 is greater than 1
    flag(msname=name[i])# Flagging bad frequency channels or antennas and autoccorelations in CASA
    gen_image(msname=name[i],imagename=imgnam1+nwms[i])# Imaging
    mkinitmodel(msname=name[i],ext=nwms[i])# Calibration:Assume a calibrator to generate a model spectrum. We are going to use the Galactic Center in our case.
    clear_cal(msname=name[i])# Calibration:Solving for delays and gain solutions
    phscal(msname=name[i])
    gen_image(msname=name[i],imagename=imgnam2+nwms[i])# Imaging again
#clean(msname,imagename=imagename,niter =500,weighting = 'briggs',robust =0,imsize =[512 ,512] ,cell=['500 arcsec'] ,mode='mfs',nterms =1,spw='0:150~900',stokes='IQUV')

#%% CREATE .NPZ FILE IN CASA
MSBCALlist=glob.glob('*MSB.cal')

i=0
for i in np.arange(len(MSBCALlist)):
    tb.open(MSBCALlist[i])
    gain = tb.getcol('CPARAM') # use 'FPARAM' when using K.cal files
    np.savez(nwms[i]+'.npz',gains=gain)

d = np.load(nwms[i]+'.npz')
d.keys()

plt.imshow(np.abs(gain[0,:,:]).T, aspect='auto', interpolation='nearest');plt.colorbar();plt.show()
plt.plot(np.abs(gain[0,:,0]));plt.show()

#%%
# You can export .image file as a .fits file then read fits file into python normally
exportfits('imagename.image','fitsimagename.fits') # Done in CASA for a single file

imagename = glob.glob('clean2*.image')
for im in imagename:
    ext = '.fits'
    fitsimagename = im.strip('.image')
    exportfits(imagename,fitsimagename+ext)
fitsimage = glob.glob(clean2*.fits)

# This section is to be done in python
test=fits.open('new_name.fits')
test.info()
test[0].header

gen_image(msname=name[0],imagename=imgnam2+nwms[0])
gen_image(name[1],imgnam1+nwms[1])
gen_image(msname=name[1],imagename=imgnam2+nwms[1])
gen_image(name[2],imgnam1+nwms[2])
gen_image(msname=name[2],imagename=imgnam2+nwms[2])
gen_image(name[3],imgnam1+nwms[3])
gen_image(msname=name[3],imagename=imgnam2+nwms[3])
