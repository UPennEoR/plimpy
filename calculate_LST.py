from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy import units as u

def calculate_LST(time_in,lon):
    loc = EarthLocation(lat=0*u.deg,lon=lon,height=0*u.m)
    time = Time(time_in,location=loc)
    lst = time.sidereal_time('mean')
    
    return lst


"""
#EXAMPLE

from calculate_LST import *

lst = calculate_LST('2017-01-02 09:06:11.5','21d25m42.3s')
print lst
#17h20m41.9048s

jd = Time(2457755.87936,format='jd')
jd
#<Time object: scale='utc' format='jd' value=2457755.87936>
lst = calculate_LST(jd,'21d25m41.0s')
lst
#<Longitude 17.34639899731436 hourangle>
print lst
#17h20m47.0364s


# Seach for Files that have GC at Transit
$ lst_select.py --ra=17.4_17.6 -C hsa7458_v001 zen.*xx.*uv
"""
