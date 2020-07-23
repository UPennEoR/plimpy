#!/usr/bin/env python
# coding: utf-8

# In[44]:


import numpy as np
import matplotlib.pyplot as plt
import glob
import os
#import sys
from pyuvdata import UVData #, UVCal
#import hera_cal as hc
#from pyuvdata import utils as uvutils
import copy
#import uvtools as uvt
#from hera_cal.datacontainer import DataContainer
#import aipy
#import operator
from astropy import units as u
#from astropy import constants as c
import pandas as pd
from astropy.coordinates import Angle


# In[65]:


JD = '2457549'
#rawpath = '/lustre/aoc/projects/hera/plaplant/HERA19Golden/RawData/'+JD+'/'
path = '/lustre/aoc/projects/hera/jaguirre/HERA19Summer2020/RawData/'+JD+'/'


# In[66]:


summary = pd.read_csv(os.path.join(path,JD+'_summary.csv'))


# In[67]:


filenames = np.array(summary['Filename'])
lsts = np.array(summary['LST (rad)'])
jds = np.array(summary['Julian Date'])


# In[68]:


# Set up LST grid
lst_bin_size = 6.*u.min
n_lst_bins = int((24.*u.hr/lst_bin_size).to(u.dimensionless_unscaled).value)
lst_edges = np.linspace(0, 2.*np.pi, n_lst_bins+1, endpoint=True)
lst_start = lst_edges[0:-1]
lst_end = lst_edges[1:]
lst_mid = (lst_start + lst_end)/2.

assert len(lst_start) == n_lst_bins
assert len(lst_end) == n_lst_bins
assert np.isclose((lst_end-lst_start)*24.*60./(2.*np.pi), lst_bin_size.value).sum() == n_lst_bins


# In[69]:


one_second = 1/(24.*3600.)


# In[70]:


get_ipython().run_cell_magic('time', '', "# Probably want to calculate how many time samples *should* appear ...\nmin_time_samp = 1\nfor i in np.arange(n_lst_bins):\n    lst_min = lst_start[i]\n    lst_max = lst_end[i]\n    lst_range = np.logical_and(lsts >= lst_min, lsts < lst_max)\n    #print(lst_min*24./(2.*np.pi), lst_range.sum())\n    \n    if lst_range.sum() > min_time_samp:\n        lst_min_str = Angle(lst_min, u.rad).to_string(unit=u.hr)\n        lst_max_str = Angle(lst_max, u.rad).to_string(unit=u.hr)\n        lst_filename = 'lst.'+JD+'.'+lst_min_str+'.'+lst_max_str+'.uvcRP.drift.uvh5'\n        print(lst_filename)\n        jd_to_select = jds[lst_range]\n        files_to_read = path + np.unique(filenames[lst_range])\n        uvd = UVData()\n        uvd.read(files_to_read)\n        lsts_in_file = np.unique(uvd.lst_array)\n        uvd.select(time_range=[jd_to_select.min()-one_second, jd_to_select.max()+one_second])\n        uvd.write_uvh5(os.path.join(path, lst_filename), clobber=True)\n        # Gonna make this just drift for now\n        #uvd.phase(lst_mid, lat, allow_rephase = True)")


# In[ ]:




