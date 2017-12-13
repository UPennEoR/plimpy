# -*- coding: utf-8 -*-
"""
    Created on Thu Nov 02 10:09:49 2017
    
    @author: tashaleebillings
    """
#printf "\e[?2004l"
import numpy as np
from pyuvdata import UVData
import glob, os

#-------------------------------------------------------------------------------------------------------
#IF YOU HAVE ONE POLARIZATION USE THIS
#-------------------------------------------------------------------------------------------------------
pollist = glob.glob('zen.2457548.46619.*.HH.uvcR')
for pol in pollist:
    os.system('python add_uvws.py '+ pol+' -C hsa7458_v001')
    os.system('python miriad2uvfits.py ' + pol+'U')
#-------------------------------------------------------------------------------------------------------
#IF YOU HAVE ONE GROUP OF POLARIZATIONS (4 polarizations) USE THIS
#-------------------------------------------------------------------------------------------------------
pollist = glob.glob('zen.2457548.46619.*.HH.uvcR')
for pol in pollist:
    os.system('python add_uvws.py '+ pol+' -C hsa7458_v001')

files = ['zen.2457548.46619.xx.HH.uvcR','zen.2457548.46619.yy.HH.uvcR','zen.2457548.46619.xy.HH.uvcR','zen.2457548.46619.yx.HH.uvcR'] #['zen.2457548.*.{xx,yy,xy,yx}.HH.uvcRU'] you can't use glob.glob because it will mess up the order of the polarizations
uv = UVData()
uv.read_miriad(files)
uv.write_miriad('zen.2457548.46619.HH.uvcU')
#-------------------------------------------------------------------------------------------------------
#IF YOU HAVE MORE THAN ONE GROUP OF POLARIZATIONS USE THIS
#-------------------------------------------------------------------------------------------------------

pollist = glob.glob('zen.2457548.46619.*.HH.uvcR')
for pol in pollist:
    os.system('python add_uvws.py '+ pol+' -C hsa7458_v001')

#You have to make list and lst.
#list1= ['zen.2457548.*.{xx,yy,xy,yx}.HH.uvcRU'] you can't use glob.glob because it will mess up the order of the polarizations
# lists=[list1,list2,list3,....] these contain your list of grouped miriad files which will always have length 4 since there are only 4 polarizations.

for lst in lists:
    mf = lst[0].strip('xx.HH.uvcRU')
    files=lst
    uv = UVData()
    uv.read_miriad(files)
    uv.write_miriad(mf+'.HH.uvcRU')
    os.system('python miriad2uvfits.py ' + mf+'.HH.uvcRU')

