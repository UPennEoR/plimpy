import numpy as np
from astropy.io import fits
from pyuvdata import UVData
from glob import glob
import os

def within_bounds(di, min_brightness, dra, ddec, ra, dec, max_angle):
    # a = angular distance on a sphere from ra and dec
    cosa = np.sin(ddec)*np.sin(dec) + np.cos(ddec)*np.cos(dec)*np.cos(dra - ra)
    return np.logical_and(cosa > np.cos(max_angle), di > min_brightness)

path = ''
impath = ''

#path = '/lustre/aoc/projects/hera/jaguirre/PolarizedMosaics/2458098/'
#impath = '/lustre/aoc/projects/hera/aseidel/test_script/'

images = [x[len(path):-7] for x in glob(path+'*.uvfits')]
images.sort()

# check this
dec0 = -30.8

data_hdul = fits.open('asu.fit')
data = data_hdul[2].data

for img in images:
    imnamepath = impath + img
    uvd = UVData()
    uvd.read_uvfits(path + img + ".uvfits", read_data=False, keep_all_metadata=True)
    ra0 = np.rad2deg(np.median(np.unique(uvd.lst_array)))

    all_angles = [np.deg2rad(x) for x in [data['RAJ2000'], data['DEJ2000'], ra0, dec0, 10]]
    in_bounds = data[within_bounds(data['Fintwide'], 1, *all_angles)]

    # additional sources not in GLEAM [ra, dec, intensity]
    fornax = np.array([[(3 + 24         /60)/24 * 360, -(37 + 16/60), 260],
                       [(3 +(21 + 40/60)/60)/24 * 360, -(37 + 10/60), 490],
                       [(3 +(22 + 43/60)/60)/24 * 360, -(37 +(12+2/60)/60), 2]])
    f_angles = [np.deg2rad(x) for x in [fornax[:,0], fornax[:,1], ra0, dec0, 30]]
    f_inb = fornax[within_bounds(fornax[:,2], 1, *f_angles)]

    make_circles = lambda d : ["circle[[{0}deg, {1}deg], 0.5deg]\n".format(
                                        str(x[0]), str(x[1])) for x in d]
    mask_file = open(imnamepath+".masks.txt", "w")
    mask_file.write("#CRTFv0\n")
    mask_file.writelines(make_circles(in_bounds) + make_circles(f_inb))
    mask_file.close()

    summary = open(imnamepath + ".info.txt", "w")
    summary.write(path + img[:32] + ".uvfits\n")
    summary.write(str(ra0) + "\n")
    summary.write(str(dec0) + "\n")
    summary.close()

os.system("casa -c image_maker_casa.py \"%s\"" % impath)
