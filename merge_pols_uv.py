from pyuvdata import UVData
import optparse
import os,sys

o=optparse.OptionParser()
o.set_usage('python merge_pols.py prefix')
o.set_description(__doc__)
opts,args = o.parse_args(sys.argv[1:])

# Making this general is a huge pain in the ass, because one has to
# parse filenames, or deduce polarization from the file
prefix = args[0]

pols = ['xx','yy','xy','yx']

files = []
for p in pols:
    files.append(prefix+'.'+p+'.HH.uv')

uv = UVData()
print 'Reading data'
for f in files:
    print f
uv.read_miriad(files)
uvfitsfile = prefix+'.HH.uv'
print 'Writing data to '+uvfitsfile
uv.write_miriad(uvfitsfile)
