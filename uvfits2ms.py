import glob

uvfits = glob.glob('*.uvfits')
for uvfit in uvfits:
   msfile=uvfit.strip('uvfits') + 'MS'
   importuvfits(vis=msfile,fitsfile=uvfit)
