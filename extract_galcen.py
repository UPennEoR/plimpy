#! /usr/bin/env python2.7
# -*- mode: python; coding: utf-8 -*-
"""
Code for extracting a ~10 minute file surrounding the galactic center transit.

Caveats: this script assumes a specific directory structure (it must be run from
the upenn_projects folder, and applies a filter for getting the directory
names), and the LST catalog-maker from wedgie to have been run already. It also
assumes that the data contains the LST of the galactic center is somewhere in
the files, and far enough away from the edges that it can extract enough time
samples (that is to say, no bounds- checking is done by the software). There is
also no file checking to ensure that output files do not exist already. But the
script is more comments than code at this point, so hopefully it's easy
enough to follow.

"""
from __future__ import print_function, division, absolute_import
import numpy as np
import os
import re
from pyuvdata import UVData
import catalog_LST

# GC coordinates in J2000, in fractional hours
gc_ra = 17 + 45/60 + 40.04/3600

# Length of a sidereal day
sd = 23 + 56/60 + 4.0916/3600

# GC coordinates in LST
gc_lst = gc_ra / sd * 2 * np.pi
gc_min = gc_lst - 0.001
gc_max = gc_lst + 0.001

# loop over days
jd_dirs = [d for d in sorted(os.listdir(os.getcwd())) if '24575' in d]
for jd in jd_dirs:
    # status update
    print(jd)
    abspath = os.path.abspath(jd)
    basename = os.path.basename(jd)
    os.chdir(jd)

    # use the LST catalog to find the relevant files
    # assumes the catalogs have already been built
    lst_rng = "{0}_{1}".format(gc_min, gc_max)
    fns = catalog_LST.find_LST(lst_rng)

    # loop over polarizations
    for pol in ['xx', 'yy', 'xy', 'yx']:
        fns_dir = [d for d in sorted(os.listdir(os.getcwd()))
                   if d[-10:] == '{}.HH.uvcR'.format(pol)]
        fn_pol = re.sub('\.xx\.', '.{}.'.format(pol), fns[0][1])
        idx_fn = fns_dir.index(fn_pol + 'R')

        idx_min = idx_fn - 1
        idx_max = idx_fn + 1

        # read in all files
        uvd = UVData()
        print("Reading {}, {}, {}".format(*fns_dir[idx_min:idx_max+1]))
        uvd.read_miriad(fns_dir[idx_min:idx_max+1])

        # get the right LST entry
        idx_lst = np.argmin(np.abs(np.unique(uvd.lst_array) - gc_lst))

        # get bracketing indices
        imin = idx_lst - 30
        imax = idx_lst + 30

        # select the data
        times = np.unique(uvd.time_array)[imin:imax]
        uvd.select(times=times)

        # save the data
        fn_out = 'gc.{0}.{1}.uvcR'.format(basename, pol)
        print("Saving {}".format(fn_out))
        uvd.write_miriad(fn_out)

        # clean up
        del uvd

    # go back up a level
    os.chdir('..')
