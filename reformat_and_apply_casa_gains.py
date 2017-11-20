#!/usr/bin/env python
"""
This script reformats npz gains tables from CASA into dictionary
objects, which can subsequently be read-in to pyuvdata and written
into calfits files, and applied to a MIRIAD file.
"""
import numpy as np
import sys, os, optparse
from hera_cal.cal_formats import AbsCal

o = optparse.OptionParser()
o.add_option('-p','--pol',dest='pol',default=None,help='polarization of intended calfits file.')
o.add_option('--miriad',dest='miriad_file',default=None,help='miriad file to calibrate')
o.add_option('--x_npz',default=None,help='gains npz derived for the x polarization')
o.add_option('--y_npz',default=None,help='gains npz derived for the y polarization')
o.add_option('--xy_npz',default=None,help='gains npz derived for the x and y polarizations together (assumes X in the first dimension, Y in the second).')
o.add_option('--cal_scheme',default='divide',help='multiply or divide by gain solutions?')
o.add_option('--ex_ants',default=None,help='antenna numbers to exclude from the calfits files.')
o.add_option('--apply',action='store_true',help='Toggle run omni_apply after file reformatting stage.')
o.add_option('--path2hera_cal',default='/home/saulkohn/githubs/hera_cal/',help='If opts.apply is toggled, a path to the hera_cal repo is required.')
opts,args = o.parse_args(sys.argv[1:])

xgo,ygo=False,False
abscal_list = []
if not opts.ex_ants:
    ex_ants = []
else:
    ex_ants = opts.ex_ants.split(',')
    ex_ants = [a.strip() for a in ex_ants]

assert opts.pol, 'polarization required'
assert opts.pol in opts.miriad_file, 'polarization not present in MIRIAD filename provided.'

if not opts.xy_npz:
    if 'x' in opts.pol:
        assert opts.x_npz, 'provide xy_npz or x_npz for polarization %s'%opts.pol
        xcal = np.load(opts.x_npz)['gains'][0,:,:]
        nants = xcal.shape[2]
        fx_out = '{0}/cal.{1}'.format(os.path.dirname(opts.x_npz),os.path.basename(opts.x_npz))
        abscal_list.append(fx_out)
        xgo = True
    if 'y' in opts.pol:
        assert opts.y_npz, 'provide xy_npz or y_npz for polarization %s'%opts.pol
        ycal = np.load(opts.y_npz)['gains'][0,:,:]
        nants = ycal.shape[2]
        fy_out = '{0}/cal.{1}'.format(os.path.dirname(opts.y_npz),os.path.basename(opts.y_npz))
        abscal_list.append(fy_out)
        ygo = True
else:
    if opts.x_npz or opts.y_npz:
        print('single and dual-pol solutions provided. Defaulting to use dual-pol')
    xycal = np.load(opts.xy_npz)['gains']
    xcal = xycal[0,:,:]
    ycal = xycal[1,:,:]
    nants = xycal.shape[2]
    fx_out = '{0}/cal.xx.{1}'.format(os.path.dirname(opts.xy_npz),os.path.basename(opts.xy_npz))
    fy_out = '{0}/cal.yy.{1}'.format(os.path.dirname(opts.xy_npz),os.path.basename(opts.xy_npz))
    if 'x' in opts.pol:
        abscal_list.append(fx_out)
    if 'y' in opts.pol:
        abscal_list.append(fy_out)
    xgo,ygo = True,True

if xgo:
    d = {}
    ng = np.isnan(xcal)
    xcal[ng] = 1.+0j
    for i in range(nants):
        d[str(i)] = xcal[:,i]
    print('   Saving reformatted gain npz %s...'%fx_out)
    np.savez(fx_out,**d)

if ygo:
    d = {}
    ng = np.isnan(xcal)
    ycal[ng] = 1.+0j
    for i in range(nants):
        d[str(i)] = ycal[:,i]
    print('   Saving reformatted gain npz %s...'%fy_out)
    np.savez(fy_out,**d)

ac = AbsCal(opts.miriad_file,abscal_list,opts.pol,opts.cal_scheme,ex_ants=ex_ants,append2hist=' '.join(sys.argv))
outname = opts.miriad_file+'.calfits'
print('    Writing %s...'%outname)
ac.write_calfits(outname) #XXX not clobbering, collisions possible (handled by pyuvdata)

if opts.apply:
    assert os.path.exists(opts.path2hera_cal)
    omnipath=outname
    outpath=os.path.dirname(opts.miriad_file)
    cmd = '{0}/scripts/omni_apply.py -p {1} --omnipath={2} --extension="K" --outpath={3} {4}'.format(opts.path2hera_cal,opts.pol,omnipath,outpath,opts.miriad_file)
    os.system(cmd)


