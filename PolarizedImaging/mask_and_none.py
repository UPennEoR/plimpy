import numpy as np
from astropy.io import fits
from pyuvdata import UVData
from glob import glob
import os

def within_bounds(di, min_brightness, dra, ddec, ra, dec, max_angle):
    # a = angular distance on a sphere from ra and dec
    cosa = np.sin(ddec)*np.sin(dec) + np.cos(ddec)*np.cos(dec)*np.cos(dra - ra)
    return np.logical_and(cosa > np.cos(max_angle), di > min_brightness)

mspath = '/lustre/aoc/projects/hera/nkern/CHAMP_Bootcamp/Lesson10_HERADataPartII/data/zen.2458116.24482.xx.HH.uvOCR.ms'
impath = '/lustre/aoc/projects/hera/aseidel/mask_and_none_test/'

os.system("rm -rf {}*".format(impath))
os.system("casa -c \"tclean(vis='{}',imagename='{}no_deconvolution',".format(mspath, impath) + 
                "niter=0, imsize = [512,512], cell=['500 arcsec']," +
                "specmode='mfs', spw='0:100~920', stokes='IQUV'," + 
                "interactive=False, pblimit=-1, gridder='widefield') | " +
                "exportfits(imagename='{}no_deconvolution.image',fitsimage='{}')\"".format(
                    impath, impath+'kerndata.fits'))

header = fits.open(impath+'kerndata.fits')[0].header
ra0, dec0 = header['CRVAL1'], header['CRVAL2']

data_hdul = fits.open('/lustre/aoc/projects/hera/aseidel/summer_2019/asu.fit')
data = data_hdul[2].data

for minb in [0.5, 1.0]:
    all_angles = [np.deg2rad(x)
                    for x in [data['RAJ2000'], data['DEJ2000'], ra0, dec0, 10]]
    in_bounds = data[within_bounds(data['Fintwide'], minb, *all_angles)]

    # additional sources not in GLEAM [ra, dec, intensity]
    fornax = np.array([[(3 + 24         /60)/24 * 360, -(37 + 16/60), 260],
                       [(3 +(21 + 40/60)/60)/24 * 360, -(37 + 10/60), 490],
                       [(3 +(22 + 43/60)/60)/24 * 360, -(37 +(12+2/60)/60), 2]])
    f_angles = [np.deg2rad(x)
                    for x in [fornax[:,0], fornax[:,1], ra0, dec0, 30]]
    f_inb = fornax[within_bounds(fornax[:,2], 1, *f_angles)]

    make_circles = lambda d : ["circle[[{0}deg, {1}deg], 0.5deg]\n".format(
                                        str(x[0]), str(x[1])) for x in d]
    mask_file = open(impath+"{:.2f}Jymask.txt".format(minb), "w")
    mask_file.write("#CRTFv0\n")
    mask_file.writelines(make_circles(in_bounds) + make_circles(f_inb))
    mask_file.close()

summary = open(impath + "info.txt", "w")
summary.write(str(ra0) + "\n")
summary.write(str(dec0) + "\n")
summary.close()

os.system("casa -c mask_and_none_casa.py \"{}\" \"{}\"".format(mspath, impath))
