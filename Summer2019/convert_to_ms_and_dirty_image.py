from glob import glob
import os, sys, time

#sys.path.append('/lustre/aoc/projects/hera/jaguirre/gitrepos/plimpy/PolarizedMosaics/')

#JD = '2458098'

JD = os.getenv('JD')
inpath = os.getenv('INPATH')
outpath = os.getenv('OUTPATH')

inpath = inpath+JD+'_sliced/'
outpath = outpath+JD+'_sliced/'

class Timer():
    
    def start(self):
        self.t0 = time.time()
        #print(self.t0)
        
    def stop(self, message):
        self.t1 = time.time()
        print(message, self.t1 - self.t0, 'sec')

timer = Timer()

uvfitsfiles = glob(inpath+'zen*.uvfits')
uvfitsfiles.sort()

def get_filestem(fullfilepath):
    fparts = os.path.basename(fullfilepath).split('.')
    fs = ''
    for fpart in fparts[0:-1]:
        fs += fpart+'.'
    return fs

for uvfitsfile in uvfitsfiles:
    filestem = get_filestem(uvfitsfile)
    
    msfile = outpath+filestem+'ms'
    print(uvfitsfile)
    timer.start()
    importuvfits(fitsfile = uvfitsfile, vis = msfile)
    timer.stop('Converting to MS')
