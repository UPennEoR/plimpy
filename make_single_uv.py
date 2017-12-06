from pyuvdata import UVData
import glob

files = glob.glob('*U') #[xx,yy,xy,yx]

uv = UVData()
uv.read_miriad(files)
uv.write_miriad('zen.2457548.46619.HH.uvcU')
