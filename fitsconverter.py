# -*- coding: utf-8 -*-
"""
    Created on Thu Nov 02 10:09:49 2017
    
    @author: tashaleebillings
    """

import numpy as np
from pyuvdata import UVData
import glob
import os

pollist = glob.glob('zen.2457548.46619.*.HH.uvcR')
for i in np.arange(len(pollist)):
    os.system('python add_uvws.py '+ pollist[i]+' -C hsa7458_v001')
    os.system('python miriad2uvfits.py ' + pollist[i]+'U')
i=0
#Now make a single miriad file
files = glob.glob('*U')

uv = UVData()
uv.read_miriad(files)
uv.write_miriad('zen.2457548.46619.HH.uvcU')
#OR IF YOU HAVE A LONG LIST USE THIS
#You have to make list and lst.
#list1= glob.glob('zen.--.--.{xx,yy,xy,yx}.HH.uvcRU')
# lst=[list1,list2,list3,....] these contain your list of grouped miriad files which will always have length 4 since there are only 4 polarizations.
for i in np.arange(len(lst)):
    mf = lst[i][0].strip('xx.HH.uvcRU')
    files=lst[i]
    uv = UVData()
    uv.read_miriad(files)
    uv.write_miriad(mf+'.HH.uvcRU')
    os.system('python miriad2uvfits.py ' + mf+'.HH.uvcRU')
