import numpy as np
from astropy.io import fits
from astropy import wcs
from pyuvdata import UVData
from glob import glob
import os


uvh5path = '/lustre/aoc/projects/hera/jaguirre/HERA19Summer2020/Simulation/zen.2457755.74025.uvCP.uvh5'
impath = '/lustre/aoc/projects/hera/aseidel/mask_and_none_simulated/'

os.system("rm -rf {}*".format(impath))

uvd = UVData()
uvd.read_uvh5(uvh5path)

name = '.'.join(uvh5path.split('/')[-1].split('.')[:-1])
uvfitspath = impath + name + '.uvfits'
mspath = impath + name + '.ms'

uvd.write_uvfits(uvfitspath, spoof_nonessential=True, force_phase=True)

#mspath = '/lustre/aoc/projects/hera/plaplant/HERA19Golden/CalibratedData/2457548/zen.2457548.26437.HH.uvcRPCS.MS'
#mspath = '/lustre/aoc/projects/hera/nkern/CHAMP_Bootcamp/Lesson10_HERADataPartII/data/zen.2458116.24482.xx.HH.uvOCR.ms'

ra0 = np.median(np.unique(uvd.lst_array))
dec0 = uvd.phase_center_dec

#generate image without deconvolution
os.system("casa -c \"importuvfits(fitsfile='{}',vis='{}') ; ".format(uvfitspath, mspath) +
                "tclean(vis='{}',imagename='{}no_deconvolution',".format(mspath, impath) + 
                "niter=0, imsize = [512,512], cell=['500 arcsec']," +
                "specmode='mfs', spw='0:100~920', stokes='IQUV'," + 
                "interactive=False, pblimit=-1, gridder='widefield'," + 
                "phasecenter='J2000 {}deg {}deg') ; ".format(ra0, dec0) +
                "exportfits(imagename='{}no_deconvolution.image',fitsimage='{}')\"".format(
                    impath, impath+'no_deconvolution.fits'))

"""
#make wcs solution
header = fits.open(impath+'no_deconvolution.fits')[0].header
w = wcs.WCS(header, naxis=2)

#open the primary beam fits file
primary_beam = fits.open("muellerbeam.fits")[0].data[0,0,0,:,:]

data_hdul = fits.open('/lustre/aoc/projects/hera/aseidel/summer_2019/asu.fit')
gleam_data = data_hdul[2].data


def effective_intensity(ra, dec, intensity):   
    # get primary beam at each source within the picture
    pbx, pby = np.array(w.all_world2pix(ra, dec, 0)).astype(int)
    pb_inb = np.logical_and.reduce(
        (pbx >= 0, pby >= 0, pbx < 512, pby < 512))
    pb = primary_beam[pbx * pb_inb, pby * pb_inb] * pb_inb
    
    return intensity * pb 


def within_bounds(ra, dec, intensity, min_brightness):
    return effective_intensity(ra, dec, intensity) > min_brightness


make_circles = lambda d : ["circle[[{0}deg, {1}deg], 0.5deg]\n".format(
                                        str(x[0]), str(x[1])) for x in d]

minbs = [0.25, 0.5] #min effective brightness (Jy)

for minb in minbs:
    in_bounds = gleam_data[within_bounds(gleam_data['RAJ2000'], gleam_data['DEJ2000'], gleam_data['Fintwide'], minb)]

    # additional sources not in GLEAM [intensity, ra, dec]
    fornax = np.array([[(3 + 24         /60)/24 * 360, -(37 + 16/60), 260],
                       [(3 +(21 + 40/60)/60)/24 * 360, -(37 + 10/60), 490],
                       [(3 +(22 + 43/60)/60)/24 * 360, -(37 +(12+2/60)/60), 2]])
    f_inb = fornax[within_bounds(*fornax.T, minb)]

    name = impath + "new_{:.2f}Jymask.txt".format(minb)
    os.system('rm -f {}'.format(name))
    mask_file = open(name, "w")
    mask_file.write("#CRTFv0\n")
    mask_file.writelines(make_circles(in_bounds) + make_circles(f_inb))
    mask_file.close()
"""


mask_file = no_mask = open(impath+"10deg_mask.txt", "w")
mask_file.write("#CRTFv0\n")
mask_file.write("circle[[{0}deg, {1}deg], 10deg]\n".format(
                                        str(ra0), str(dec0)))
mask_file.close()

no_mask = open(impath+"no_msk.txt", "w")
no_mask.write("#CRTFv0\n")
no_mask.write("circle[[{0}deg, {1}deg], 50deg]\n".format(
                                        str(ra0), str(dec0)))
no_mask.close()
    
summary = open(impath + "summary.txt", "w")
summary.write(str(ra0) + "\n")
summary.write(str(dec0) + "\n")
summary.close()

os.system("casa -c mask_and_none_casa.py \"{}\" \"{}\"".format(mspath, impath))
