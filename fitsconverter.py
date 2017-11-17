# -*- coding: utf-8 -*-
"""
    Created on Thu Nov 02 10:09:49 2017
    
    @author: tashaleebillings
    """

import numpy as np
import glob
import os



"""
    The purpose of this section is to convert .uvc files to .uvfits
    This script needs to be ran in the same directory as
    """
daydir = '2457550'
datalistname = 'list.txt'
pythonscript = ['add_uvws.py','miriad2uvfits.py','uvfits2ms.py']

pathtodatadir = '/data4/paper/plaplant/HERA2015/upenn_projects/'
pathtodatatxtdir = 'visibility/'+daydir+'/'
pathtoformatscript = 'visibility/pythonscripts/'

ra = '17:45.6' # Right Ascension (h:m)
dec = '-28:56' # Declination (deg:m)

#Generate list of data names
with open(datalistname,'r') as myfile:
    datalist = myfile.read().splitlines()

daydirlist = ['2457550','2457550','2457550','2457550']
#daydirlist = ['2457548','2457549','2457550','2457551',
#              '2457551','2457552','2457552','2457553',
#              '2457554','2457555','2457548','2457549',
#              '2457550','2457551','2457551','2457552',
#              '2457552','2457553','2457554','2457555']

#dattimelist= ['.45227.','.44531.','','','','','','','','','','','','','','','','','','']

#os.system('mkdir '+pathtodatatxtdir)
os.system('cp '+pathtoformatscript+'hsa7458_v001.py '+pathtodatatxtdir )

i,j = 0,0

for i in np.arange(len(pythonscript)):
    os.system('cp '+pathtoformatscript+pythonscript[i]+' '+pathtodatatxtdir )
i=0
fitsdatalist = []

path2data=[]
for i in np.arange(len(datalist)):
    path2data.append(pathtodatadir+daydirlist[i]+'/'+datalist[i])

for j in np.arange(len(datalist)):
    os.system('cp -r '+path2data[j]+' .')
    os.system('python '+pythonscript[i]+' '+datalist[j]+' -C '+'hsa7458_v001') #.uvcRO -> .uvcROU
    os.system('python '+pythonscript[i+1]+' '+datalist[j]+'U') # .uvcROU -> .uvcROU.uvfits
    fitsdatalist.append(datalist[j]+'U.uvfits')

print('You can now switch to CASA and execute casascript.py')















"""
    2457548/zen.2457548.45227.xx.HH.uvcRO
    2457549/zen.2457549.44531.xx.HH.uvcRO
    2457550/zen.2457550.44531.xx.HH.uvcRO
    2457551/zen.2457551.43835.xx.HH.uvcRO
    2457551/zen.2457551.44531.xx.HH.uvcRO
    2457552/zen.2457552.43836.xx.HH.uvcRO
    2457552/zen.2457552.44532.xx.HH.uvcRO
    2457553/zen.2457553.43835.xx.HH.uvcRO
    2457554/zen.2457554.43139.xx.HH.uvcRO
    2457555/zen.2457555.43139.xx.HH.uvcRO
    2457548/zen.2457548.45227.yy.HH.uvcRO
    2457549/zen.2457549.44531.yy.HH.uvcRO
    2457550/zen.2457550.44531.yy.HH.uvcRO
    2457551/zen.2457551.43835.yy.HH.uvcRO
    2457551/zen.2457551.44531.yy.HH.uvcRO
    2457552/zen.2457552.43836.yy.HH.uvcRO
    2457552/zen.2457552.44532.yy.HH.uvcRO
    2457553/zen.2457553.43835.yy.HH.uvcRO
    2457554/zen.2457554.43139.yy.HH.uvcRO
    2457555/zen.2457555.43139.yy.HH.uvcRO
    """
