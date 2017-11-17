from pyuvdata import UVData

files = ['zen.2457548.45923.xx.HH.uvcU',
         'zen.2457548.45923.yy.HH.uvcU',
         'zen.2457548.45923.xy.HH.uvcU',
         'zen.2457548.45923.yx.HH.uvcU']

uv = UVData()
uv.read_miriad(files)
uv.write_miriad('zen.2457548.45923.HH.uvcU')
