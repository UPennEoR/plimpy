#all imports
import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import sys
from pyuvdata import UVData, UVCal
import hera_cal as hc
from hera_cal.data import DATA_PATH
from collections import OrderedDict as odict
from pyuvdata import utils as uvutils
import copy
import uvtools as uvt
from hera_cal.datacontainer import DataContainer
import aipy
import operator
import astropy

path = sys.argv[0]
pathToCalibrated = sys.argv[1]
J = sys.argv[2]
JD = str(J)

'''def contain_final_flags(path):
    flagged = False
    
    flags = glob.glob(path + "zen.*.final_flags.h5")
    if len(flags) == 0: flagged = False
    else: flagged = True
        
    return flagged'''

def imaged_files (path):
    imagefits = glob.glob(str(path) + 'zen.*.HH.calibrated.uvh5_image/*.image.image.fits')
    uvh5files = []
    for imf in imagefits:
        tmp = os.path.split(imf)[-1].split('.')
        uvh5files.append(str(path) + tmp[0]+'.'+ tmp[1]+'.'+tmp[2]+'.HH.uvh5')
    return uvh5files

def calibration(pathToData, pathToCalibrated):
    cmd = "apply_cal.py --new_cal={} --filetype_in uvh5 --filetype_out uvh5 {} {}"
    
    pathToData = str(path)
    JDs = []
    im_files = imaged_files(str(path))
    for f in im_files:
        JD = os.path.split(f)[0].split('zen.')[1].split('.HH')[0]
        JDs.append(JD)

    for f in JDs:
        uvh5 = str(pathToData) + "zen." + str(f) + ".HH.uvh5"
        cal =  pathToCalibrated + "zen." + str(f) + ".HH.calibrated.uvh5"
        refca = str(pathToData) + "zen." + str(f) + ".HH.smooth_abs.calfits'"
    
        os.system(cmd.format(refca, uvh5, cal))

def make_dictionary (JD, path):
    
    #load in all files for that day
    filenames = sorted(glob.glob(os.path.join(pathToCalibrated, JD + '/zen.' + JD +'*.HH.calibrated.uvh5')))
    
    #create dictionary
    master_dict = {}
    
    for file in filenames:
        uvd = UVData()
        uvd.read(file, read_data=False, read_metadata=True)
        LST = np.unique(uvd.lst_array)
        JD = np.unique(uvd.time_array)
    
        for (J, L) in zip(JD, LST): 
            master_dict.update({(J, L): file})
        
    return master_dict

def find_files(JD, LSTstart, LSTend, master_Dict=None):
    files = []
    
    if master_Dict is None:
        master_Dict = make_dictionary (JD, str(path))
    else:
        if type(master_Dict) is not dict:
            raise TypeError
    
    dict_list = list(master_Dict)
    for i,(J,LST) in enumerate(dict_list):
        if ((min(LSTstart, LSTend) <= dict_list[i][1] and dict_list[i][1] <= max(LSTstart, LSTend))):
            files.append(list(master_Dict.values())[i])
    
    return np.unique(files)

def checkImageFiles(foundFiles, imagedFiles):
    check = False
    for f in foundFiles:
        for fi in imagedFiles:
            if ((f == fi)): 
                check = True
    return check

def max_min_LSTs(dict_list):
    Max = 0
    Min = 1e99
    for i, (J, LST) in enumerate(dict_list):
        if (Max < dict_list[i][1]): 
            Max = dict_list[i][1]
        if (Min > dict_list[i][1]): 
            Min = dict_list[i][1]
    return Max, Min

def slice_day(JD, interval):
    master_Dict = make_dictionary (JD, str(path))
    dictionary_list = list(master_Dict)
    
    #max, min LST values
    Max, Min = max_min_LSTs(dictionary_list)
    
    sliced = {}
    
    #this for-loop isn't correct it doesn't include all the ranges. what about last value?
    for i in np.arange(Min, Max, interval):
        sliced.update({(JD, i, i+interval): find_files(JD, i, i+interval, master_Dict)})
        
    sliced.update({(JD, i, i+interval): find_files(JD, i, Max, master_Dict)})
    
    return sliced, master_Dict

def splitDay (JD, interval):
    pathname = str(path) + JD + '/'
    sliced, master_dict = slice_day (JD, interval)
    imagefiles = imaged_files(pathname)
    
    dict_list = list(master_dict)
    
    radDay = 2*np.pi
    for i in np.arange(0, radDay, interval):
        files = find_files(JD, i, i+interval, master_dict)
        fi_li = list(files)
        
        #check if files have imaged files
        check = checkImageFiles(fi_li, imagefiles)
        
        if len(fi_li) == 0 or (not check):
            continue
        else:
            uvd = UVData()
            uvd.read(list(files), read_data=True, read_metadata=True)
            LST = np.unique(uvd.lst_array)
            times = np.unique(uvd.time_array)

            indexs = np.logical_and(LST >= i, LST <=i+interval)

            oneNfiftyScs = (110./60.)/(24*60.)*2*np.pi

            if (LST[indexs].max() - LST[indexs].min() >= oneNfiftyScs):
                uvd.select(times=times[indexs])
                uvd.write_uvfits('zen.' + JD + '.' + str(i) + '-' + (str(i+interval)) + '.HH.calibrated.uvfits', 
                                 force_phase=True, spoof_nonessential=True)

#calls
#flagbool = contain_final_flags(path)

#constant(s)
interval = 2./(24*60.)*2*np.pi

calibration(path, pathToCalibrated)
splitDay(JD, interval)
