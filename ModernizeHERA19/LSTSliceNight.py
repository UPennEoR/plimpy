#!/usr/bin/env python
# coding: utf-8

# An example for finding raw data which can be interpolated onto a specific LST range.

import numpy as np
#import matplotlib.pyplot as plt
from pyuvdata import UVData
from scipy.interpolate import CubicSpline
from astropy import units as u, constants as c
import pandas as pd
import os
from astropy.coordinates import Angle
import time
from LSTSlicer import BuildNewUVData


class Timer():
    
    def start(self):
        self.t0 = time.time()
        #print(self.t0)
        
    def stop(self, message):
        self.t1 = time.time()
        print(message, self.t1 - self.t0, 'sec')


# Define the LST grid.

# Set up LST grid; in radians
# A trifle worried about numerical precision, but let's plunge ahead

# Files will be this long
lst_bin_size = 6.*u.min
n_lst_bins = int((24.*u.hr/lst_bin_size).to(u.dimensionless_unscaled).value)
lst_edges = np.linspace(0, 2.*np.pi, n_lst_bins+1, endpoint=True)
lst_start = lst_edges[0:-1]
lst_end = lst_edges[1:]
lst_mid = (lst_start + lst_end)/2.

assert len(lst_start) == n_lst_bins
assert len(lst_end) == n_lst_bins
assert np.isclose((lst_end-lst_start)*24.*60./(2.*np.pi), lst_bin_size.value).sum() == n_lst_bins

# Need to specify the time sampling within a "bin"
lst_time_sample = 10*u.s # This is actually every 10 sidereal seconds
n_samples = int((lst_bin_size/lst_time_sample).to(u.dimensionless_unscaled).value)
lst_sampling = ((np.arange(n_samples)*lst_time_sample)/(24*u.hr)).to(u.dimensionless_unscaled).value*2.*np.pi

# 40 LST seconds ... a bit arbitrary
pad = 40/3600.*2.*np.pi/24.

# This is highly specific to running at NRAO; should probably give an option to
# specify paths to the script
if JD == '2457756'
    inpath = '/lustre/aoc/projects/hera/jaguirre/HERA19Summer2020/Simulation/'
    outpath = '/lustre/aoc/projects/hera/jaguirre/HERA19Summer2020/LSTSlice/Simulation/'
else:
    inpath = '/lustre/aoc/projects/hera/jaguirre/HERA19Summer2020/CalibratedData/'+JD+'/'
    outpath = '/lustre/aoc/projects/hera/jaguirre/HERA19Summer2020/LSTSlice/'+JD+'/'

# Grab the summary file for the day in question
summary = pd.read_csv(os.path.join(inpath,JD+'_summary.csv'))
filenames_all = np.array(summary['Filename'])
lsts_all = np.array(summary['LST (rad)'])
jds_all = np.array(summary['Julian Date'])

def make_lst_str(lst_rad):
    return Angle(lst_rad, u.rad).to_string(unit=u.hr)

timer = Timer()

# Probably want to calculate how many time samples *should* appear ...
min_time_samp = 36
for i in np.arange(n_lst_bins):
    
    # We're going to pick this LST bin
    lst_min = lst_start[i]
    lst_max = lst_end[i]
    lsts_new = lst_min + lst_sampling

    lst_range = np.logical_and(lsts_all >= lst_min-pad, lsts_all < lst_max+pad)
    
    nsamples = lst_range.sum()

    if nsamples > min_time_samp:

        timer.start()
        # Output filename
        lst_min_str = Angle(lst_min, u.rad).to_string(unit=u.hr)
        lst_max_str = Angle(lst_max-lst_sampling[1], u.rad).to_string(unit=u.hr)
        lst_filename = 'lst.'+JD+'.'+lst_min_str+'.'+lst_max_str+'.uvcRP.drift.uvh5'

        print('Making file for', lst_min_str,'with',nsamples,'nsamples of raw data')
        
        files_to_read = inpath + np.unique(filenames_all[lst_range])
        uvd_raw = UVData()
        uvd_raw.read(files_to_read, time_range=[jds_all[lst_range].min(), jds_all[lst_range].max()])

        uvd_new = BuildNewUVData(uvd_raw, lsts_new)

        ants = uvd_new.get_ants()
        uvd_new.select(antenna_nums=ants)
        uvd_new.write_uvh5(os.path.join(outpath, lst_filename), clobber=True)
        timer.stop(lst_filename)





