import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys
from pyuvdata import UVData, UVCal
#import hera_cal as hc
#from hera_cal.data import DATA_PATH
#from collections import OrderedDict as odict
from pyuvdata import utils as uvutils
import copy
#import uvtools as uvt
#from hera_cal.datacontainer import DataContainer
#import aipy
#import operator
from astropy import units as u
from astropy import constants as c
import pandas as pd
from astropy.time import Time
from astropy.coordinates import get_body, SkyCoord, AltAz, EarthLocation, Angle

# Calculate sun position
def SunAlt(times, location):
    sun = get_body('sun', times)
    sun_altaz = sun.transform_to(AltAz(location = location))
    return sun_altaz.alt.deg

#JD = '2457552'
JD = str(sys.argv[1])
print('Processing JD', JD)
rawpath = '/lustre/aoc/projects/hera/plaplant/HERA19Golden/RawData/'+JD+'/'
outpath = '/lustre/aoc/projects/hera/jaguirre/HERA19Summer2020/RawData/'+JD+'/'
pols = ['xx','yy','xy','yx']

#load in all files for that day
filenames = sorted(glob.glob(os.path.join(rawpath, 'zen.' + JD +'*xx.HH.uvcRP')))

# There's no point in doing this multiple times
uvd = UVData()
uvd.read(filenames[0])
# HERA Lat / Lon
lat, lon, alt = uvd.telescope_location_lat_lon_alt_degrees
latitude = Angle(lat*u.deg)
longitude = Angle(lon*u.deg)
# This is GEODETIC
karoo = EarthLocation(longitude, latitude, height=alt*u.m)

data = {'Filename': [], 
        'Julian Date': [], 
        'LST (rad)': [], 
        'UTC': [],
        'LST (hr)': [],
        'Sun alt (deg)': []
       }

for filename in filenames:
    
    filename_only = filename.split('/')[-1]
    tmp = filename_only.split('.')
    juldate = tmp[1]+'.'+tmp[2]
    print(juldate)
    preamble = rawpath+'zen.'+juldate
    #print(preamble)
    outfile = 'zen.'+juldate+'.uvcRP.uvh5'
    
    polfiles = [preamble+'.'+p+'.HH.uvcRP' for p in pols]
    calfits = preamble + '.HH.uvcRP.calfits'
       
    uvd = UVData()
    uvd.read(filename)
    ntimes = uvd.Ntimes
    lst_tmp = np.unique(uvd.lst_array)
    time_tmp = np.unique(np.unique(uvd.time_array))
    times = Time(time_tmp, format='jd')
    
    """ This approach makes a list out of however many time samples are in a file (a variable
    number, otherwise we would just declare the size of the array at the outset)
    """
    data['Filename'] += [outfile for i in np.arange(ntimes)]
    data['Julian Date'] += list(time_tmp)  
    data['LST (rad)'] += list(lst_tmp)
    data['UTC'] += list(Time(time_tmp, format='jd').iso)
    data['LST (hr)'] += list(Angle(lst_tmp, u.rad).to_string(unit=u.hour))
    data['Sun alt (deg)'] += list(SunAlt(times, karoo))
    
    uvdallpol = UVData()
    uvdallpol.read(polfiles)
    uvdallpol.write_uvh5(outpath+outfile, clobber=True)
    

dataframe = pd.DataFrame(data)
dataframe.to_csv(outpath+JD+'_summary.csv')