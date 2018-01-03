import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

gleam = Table.read('/Users/jaguirre/Data/gleam.fits')

wh = np.where(gleam['Fp158'] > 4.)[0]


f = open('gleam.reg','w')
f.write('# Region file format: DS9 version 4.1\n')
f.write('global color=red dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
f.write('fk5\n')
for w in wh:
    #text = 'circle(%.3f,%.3f,2400\") # text={%.1f Jy} \n' % (gleam['RAJ2000'][w],gleam['DEJ2000'][w],gleam['Fp158'][w])
    text = 'point(%.3f,%.3f) # text={%.1f Jy} point=diamond\n' % (gleam['RAJ2000'][w],gleam['DEJ2000'][w],gleam['Fp158'][w])
    f.write(text)
#f.write('circle(%.3f,%.3f,2400\") # text={%.1f Jy} \n')
f.close()
