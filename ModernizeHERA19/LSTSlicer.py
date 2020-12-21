import numpy as np
from pyuvdata import UVData

def jd_at_lst(uvdata, lst0):
    """ Given a UVData object, find the JD corresponding to the requested LST.  
    Uses simple linear interpolation. """
    jd = np.unique(uvdata.time_array)
    lst = np.unique(uvdata.lst_array)
    
    return np.interp(lst0, lst, jd)

def BuildNewUVData(uvd_raw, lsts_new):
    
    Ntimes = len(lsts_new)

    # I hope linear interpolation is sufficiently accurate.
    jds_new = jd_at_lst(uvd_raw, lsts_new) 
    
    lsts = np.unique(uvd_raw.lst_array) # For some God-forsaken reason, LSTs are in radians
    # We're going to need the nearest index to the original JDs 
    jds = np.unique(uvd_raw.time_array)
    jd_nearest_indx = np.zeros_like(jds_new, dtype='int64')
    for i,jd in enumerate(jds_new):
        jd_nearest_indx[i] = np.argmin(np.abs(jds-jd))
    
    # Build up the list of required attributes
    attributes = list(uvd_raw.__dict__.keys())
    for i,attribute in enumerate(attributes):
        if attribute[0] == '_':
            attributes[i] = attribute[1:]
            
     # Define new array dimensions
    Nbls = uvd_raw.Nbls
    Nblts = Ntimes * Nbls
    Nspws = uvd_raw.Nspws
    Nfreqs = uvd_raw.Nfreqs
    Npols = uvd_raw.Npols   
    
    # Define new arrays
    data_array_new = np.zeros_like(uvd_raw.data_array, shape = [Nblts, Nspws, Nfreqs, Npols])
    flag_array_new = np.zeros_like(uvd_raw.flag_array, shape = [Nblts, Nspws, Nfreqs, Npols])
    nsample_array_new = np.ones_like(uvd_raw.nsample_array, shape = [Nblts, Nspws, Nfreqs, Npols])
    uvw_array_new = np.zeros_like(uvd_raw.uvw_array, shape=[Nblts, 3])
    time_array_new = np.zeros_like(uvd_raw.time_array, shape=[Nblts])
    lst_array_new = np.zeros_like(uvd_raw.lst_array, shape=[Nblts])
    ant_1_array_new = np.zeros_like(uvd_raw.ant_1_array, shape=[Nblts])
    ant_2_array_new = np.zeros_like(uvd_raw.ant_2_array, shape=[Nblts])
    baseline_array_new = np.zeros_like(uvd_raw.baseline_array, shape=[Nblts])
    integration_time_new = np.zeros_like(uvd_raw.integration_time, shape=[Nblts])
    
    # Copy over
    for itime in np.arange(Ntimes):
        for ibl in np.arange(Nbls):
            indx = itime * Nbls + ibl
            # If in drift mode, they're all the same
            uvw_array_new[indx, :] = uvd_raw.uvw_array[ibl, :]
            time_array_new[indx] = jds_new[itime]
            lst_array_new[indx] = lsts_new[itime]
            # we're really assuming the pattern repeats every time sample ...
            ant_1_array_new[indx] = uvd_raw.ant_1_array[ibl]
            ant_2_array_new[indx] = uvd_raw.ant_2_array[ibl]  
            
    # What do we need to change?
    to_change = {'Ntimes' : Ntimes,
                 'Nblts' : Nblts,
                 'data_array' : data_array_new,
                 'flag_array' : flag_array_new, # ... hoo boy
                 'nsample_array' : nsample_array_new, # what?
                 'uvw_array' : uvw_array_new, #... need to confirm we're drift scanning
                 'time_array' : time_array_new, 
                 'lst_array' : lst_array_new, 
                 'integration_time' : integration_time_new,
                 'ant_1_array' : ant_1_array_new, # why God, why?
                 'ant_2_array' : ant_2_array_new,
                 'baseline_array' : baseline_array_new # Oh, fuck my life
                }
    
    # Define the new object
    uvd_new = UVData()

    for attribute in attributes:
        if attribute in to_change.keys():
            #print('Changing '+attribute)
            setattr(uvd_new, attribute, to_change[attribute])
        else:
            #print('Keeping '+attribute)
            setattr(uvd_new, attribute, getattr(uvd_raw, attribute))
        try:
            blt_ind1, blt_ind2, pol_ind = uvd_new._key2inds(bl)
            #print('blt_ind1', blt_ind1.shape)
        except:
            pass
            #print('Cannot get key2ind')
    
    # Hopefully robust
    uvd_new.integration_time[:] = uvd_new._calc_single_integration_time()
    
    # Don't know how to generate these magic nubmers otherwise
    uvd_new.baseline_array = uvd_new.antnums_to_baseline(uvd_new.ant_1_array, uvd_new.ant_2_array)
    
    # Some checking
    assert uvd_new.Nants_data == uvd_raw.Nants_data
    
    # OK, now go through baseline by baseline to fix up the changed data ..
    for bl in uvd_new.get_antpairpols():
        #print(bl)
        data = uvd_raw.get_data(bl)
        # There might not be a good way to handle flagging that differs in time between frequencies ...
        # i.e., I think when the size of the unflagged data differs, you might just have to loop.  Ugh.
        # And then, if there is no good data ... double Ugh.
        flags = uvd_raw.get_flags(bl)
        data_interp = np.zeros_like(data, shape=[uvd_new.Ntimes, uvd_new.Nfreqs])
        flags_interp = np.zeros_like(flags, shape=[uvd_new.Ntimes, uvd_new.Nfreqs])
    
        for ifreq in np.arange(uvd_new.Nfreqs):
            fflags = ~flags[:,ifreq].squeeze()
            flags_interp[:, ifreq] = flags[jd_nearest_indx, ifreq]
            # If there is not enough data to spline interpolate, skip
            if fflags.sum() >= 2:
                cs = CubicSpline(lsts[fflags], data[fflags, ifreq])
                data_interp[:, ifreq] = cs(lsts_new)
              
                # This was super clever for vectorizing the cubic spline across frequencies, but I don't think
                # it generalizes with different flags for every frequency
                #cs = CubicSpline(lsts, data, axis=0)
                #interp = cs(lsts_new)
    
        # Need that magic incantation for shoving a visibility's data back into the UVData object
        inds = uvd_new._key2inds(bl)
        blt_inds = inds[0]
        pol_ind = inds[2][0]
        uvd_new.data_array[blt_inds,0,:,pol_ind] = data_interp
        uvd_new.flag_array[blt_inds,0,:,pol_ind] = flags_interp
    
    return uvd_new
