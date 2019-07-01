# plimpy
PoLarized IMaging in PYthon

Store data used for analysis of polarized imaging.

# Install PIXELL
## Linux Users

(hera3) modv-ve503-1250:~ USER$ conda install astropy numpy matplotlib pyyaml h5py pillow
(hera3) modv-ve503-1250:~ USER$ conda install -c conda-forge libsharp healpy
(hera3) modv-ve503-1250:~ USER$ conda install -c anaconda automake
(hera3) modv-ve503-1250:~ USER$ git clone https://github.com/simonsobs/pixell
(hera3) modv-ve503-1250:~ USER$ cd pixell
(hera3) modv-ve503-1250:~ USER$ python setup.py install --user

Now, lets see if this works. Go back to you home directory and then try to import the modules from the command line.

(hera3) modv-ve503-1250:~ USER$ cd
(hera3) modv-ve503-1250:~ USER$ python -c "from pixell import enmap, utils"

## MAC-OS Users

The problem here is that you don't have the correct default C compiler. You need to install the correct C compiler, gcc8, and point to the correct C and Fortran compiler.

(hera3) modv-ve503-1250:~ USER$ conda install astropy numpy matplotlib pyyaml h5py pillow
(hera3) modv-ve503-1250:~ USER$ conda install -c conda-forge libsharp healpy
(hera3) modv-ve503-1250:~ USER$ conda install -c anaconda automake
(hera3) modv-ve503-1250:~ USER$ git clone https://github.com/simonsobs/pixell

#### Install C compiler
(hera3) modv-ve503-1250:~ USER$ sudo port install gcc8

#### Point to correct compilers
(hera3) modv-ve503-1250:~ USER$ cd pixell
(hera3) modv-ve503-1250:~ USER$ FC=/opt/local/bin/gfortran-mp-8 CC=/opt/local/bin/gcc-mp-8  python setup.py build_ext
python setup.py install

#### Install as usual
(hera3) modv-ve503-1250:~ USER$ python setup.py install --user

Now, lets see if this works. Go back to you home directory and then try to import the modules from the command line.

(hera3) modv-ve503-1250:~ USER$ cd
(hera3) modv-ve503-1250:~ USER$ python -c "from pixell import enmap, utils"
