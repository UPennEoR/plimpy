import numpy as np
from astropy.io import fits
from glob import glob
import os

def within_bounds(di, dra, dde, ra, de, max_angle, min_brightness):
    #anglular distance on a sphere from ra and dec isn't just the standard distance formula
    cosa = np.sin(dde)*np.sin(de) + np.cos(dde)*np.cos(de)*np.cos(dra - ra)
    return np.logical_and(cosa > np.cos(max_angle), di > min_brightness)

path = '/lustre/aoc/projects/hera/H1C_IDR2/IDR2_2/2458098/'
impath = '/lustre/aoc/projects/hera/aseidel/'

images = [x[len(path):-22] for x in glob(path+'zen.*.HH.calibrated.uvh5_image/*image.image.fits')]

for im in images:
    imnamepath = impath+im[43:]
    hdul = fits.open(path+im+".uvh5.image.image.fits")

    data_hdul = fits.open('asu.fit')

    dec0, ra0 = np.deg2rad(np.array((hdul[0].header['CRVAL2'], hdul[0].header['CRVAL1'])))

    data = data_hdul[2].data
    in_bounds = data[within_bounds(data['Fintwide'], np.deg2rad(data['RAJ2000']), np.deg2rad(data['DEJ2000']), ra0, dec0, np.deg2rad(10.), 1)]

    mask_file = open(imnamepath+".masks.txt", "w")
    mask_file.write("#CRTF\n")
    mask_file.writelines(["circle[[" + str(x[0]) + "deg, " + str(x[1]) + "deg], 0.5deg]\n" for x in in_bounds])
    mask_file.close()

    os.system("casa -c \"clean(vis='%s.ms', imagename='%s.deconvolved', niter=10000, weighting='briggs', robust=0, imsize = [512,512], pbcor=False, cell=['500 arcsec'], mode='mfs', nterms=1, spw='0:150~900', stokes='IQUV', mask='%s.masks.txt', interactive=False, npercycle=5, threshold='0.1mJy/beam')\"" % (path+im[43:], imnamepath, imnamepath))
