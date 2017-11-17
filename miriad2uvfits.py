#!/usr/bin/env python

from pyuvdata import UVData
import optparse
import os,sys

o=optparse.OptionParser()
o.set_usage('python [options] *.uv')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])

for uvfile in args:
   print 'Reading ', uvfile
   UV = UVData()
   UV.read_miriad(uvfile,run_check=False,run_check_acceptability=False)
   UV.write_uvfits(uvfile+'.uvfits',spoof_nonessential=True,force_phase=True)
