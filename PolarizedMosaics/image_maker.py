import numpy as np
from astropy.io import fits
from pyuvdata import UVData
from glob import glob
import os

def within_bounds(di, dra, dde, ra, de, max_angle, min_brightness):
    #angular distance on a sphere from ra and dec isn't just the standard distance formula
    cosa = np.sin(dde)*np.sin(de) + np.cos(dde)*np.cos(de)*np.cos(dra - ra)
    return np.logical_and(cosa > np.cos(max_angle), di > min_brightness)

#path = '/lustre/aoc/projects/hera/jaguirre/PolarizedMosaics/2458098/'
#impath = '/lustre/aoc/projects/hera/aseidel/test_script/'

path = '/home/alex/gleam clean boxes/script_test/'
impath = '/home/alex/gleam clean boxes/script_test/images/'

images = [x[len(path):-7] for x in glob(path+'zen.*.HH.calibrated.uvfits')]
images.sort()

for file in images:
    imnamepath = impath+file
    hdul = UVData().read_uvfits(path+file+".uvfits", read_data=False, read_metadata=False)
    dec0, ra0 = np.deg2rad(np.array((hdul[0].header['CRVAL2'], hdul[0].header['CRVAL1'])))

    data_hdul = fits.open('asu.fit')
    data = data_hdul[2].data
    in_bounds = data[within_bounds(data['Fintwide'], np.deg2rad(data['RAJ2000']), np.deg2rad(data['DEJ2000']), ra0, dec0, np.deg2rad(10.), 1)]

    fornax = np.array([[(3+24/60)/24*360, -(37+16/60), 260], [(3+(21+40/60)/60)/24*360, -(37+10/60), 490], [(3+(22+43/60)/60)/24*360, -(37+(12+2/60)/60), 2]])
    fornax_inb = fornax[within_bounds(fornax[:,2], np.deg2rad(fornax[:,0]), np.deg2rad(fornax[:,1]), ra0, dec0, np.deg2rad(30.), 0)]

    mask_file = open(imnamepath+".masks.txt", "w")
    mask_file.write("#CRTFv0\n")
    mask_file.writelines(["circle[[" + str(x[0]) + "deg, " + str(x[1]) + "deg], 0.5deg]\n" for x in in_bounds])
    mask_file.writelines(["circle[[" + str(x[0]) + "deg, " + str(x[1]) + "deg], 0.5deg]\n" for x in fornax_inb])
    mask_file.close()

    summary = open(imnamepath+".info.txt", "w")
    summary.write(path+file[:32]+".uvfits\n")
    summary.write(str(ra0) + "\n")
    summary.write(str(dec0) + "\n")
    summary.close()

os.system("casa -c image_maker_casa.py \"%s\"" % impath)
