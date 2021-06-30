import numpy as np
from pyuvdata import UVData
import copy
from scipy import interpolate
    
class uvdata_pol_calibrator():
    
    def __init__(self, model_data=None, real_data=None, mode="model_based"):
        """
        Parameters
        ----------
        model_data, real_data : UVData
        
        mode : str
            Specify how gains are defined. Choices are "model_based" and "data_based".
            For "model_based", V_ij_data = G_i V_ij_model G_j^H;
            For "data_based", V_ij_model = G_i V_ij_data G_j^H.
            
        """
        if mode == "model_based":
            self.base_data, self.plex_data = model_data, real_data
        if mode == "data_based":
            self.base_data, self.plex_data = real_data, model_data
        self.gain_array = np.zeros((self.base_data.Nants_data, self.base_data.Nfreqs, 2, 2)).astype(np.complex128)
        
    def data_select(self, use_all_times=False, time_range=[], use_all_frequencies=False, freq_range=[], use_all_ants=False, ants=[], bl_range=[]):
        """
        Select the subset of the data for running the calibration.
        
        Parameters
        ----------
        use_all_times : bool
        
        use_all_frequencies : bool
        
        use_all_ants : bool
        
        ants : list
            list of antenna nums
            
        time_range : list
            list of length-2, specifying the start and the end of the time range
            
        freq_range : list
            list of length-1 or length-2;
            If the length is 1, it just contains a single frequency channel.
            If the length is 2, it contains the first and the last frequency channel. 
            
        bl_range : list
            Only select baselines falling into the range of length.  
        
        """
        if use_all_ants:
            assert np.all(np.sort(np.unique(self.plex_data.ant_1_array)) == np.sort(np.unique(self.base_data.ant_1_array))), "ants must be the same."
            ants = np.sort(np.unique(self.plex_data.ant_1_array))
        else:
            assert isinstance(ants, list), "ants must be a list." 
            ants = list(np.sort(ants))
        self.ants2cal = ants
        
        if use_all_times:
            assert np.all(np.isclose(np.unique(self.base_data.time_array), np.unique(self.plex_data.time_array))), "time_arrays must be the same."
            time_range = [np.unique(self.plex_data.time_array)[0]-1e-2, np.unique(self.plex_data.time_array)[-1]+1e-2]
        else:
            assert isinstance(time_range, list), "time_range must be a list."
            assert len(time_range) == 2, "Length of time_range must be 2."
        
        if use_all_frequencies:
            assert self.base_data.Nfreqs == self.plex_data.Nfreqs, "Nfreqs must be the same."
            freqs = list(range(self.base_data.Nfreqs))
        else:
            assert isinstance(freq_range, list), "freq_range must be a list."
            if len(freq_range) == 1:
                freqs = freq_range
            elif len(freq_range) == 2:
                freqs = list(range(freq_range[0], freq_range[1]+1))
            else:
                raise ValueError("Length of freq_range must be either 1 or 2.")
        self.freqs2cal = freqs
        
        base_data_select = self.base_data.select(antenna_nums=ants, freq_chans=freqs, time_range=time_range, inplace=False) 
        plex_data_select = self.plex_data.select(antenna_nums=ants, freq_chans=freqs, time_range=time_range, inplace=False)
        
        assert isinstance(bl_range, list), "bl_range must be a list."
        if len(bl_range) != 0:
            assert len(bl_range) == 2, "Length of bl_range must be 2."
            # Form baseline length splits on the data array
            bl_length = ( np.linalg.norm(base_data_select.uvw_array[:,:2], axis=1) +np.linalg.norm(plex_data_select.uvw_array[:,:2], axis=1) )/2.
            bl_splits = (bl_length > np.min(bl_range)) * (bl_length < np.max(bl_range)) 
            base_data_select.data_array *= bl_splits[:,None, None, None]
            plex_data_select.data_array *= bl_splits[:,None, None, None]
        
        # We are now reordering the data_array
        # original shape of data_array: (Nblts, Nspws, Nfreqs)
        # intermidiate shape: (Nants, Nants, Nspws, Ntimes, Nfreqs, Npols)
        # final shape: (Nants, Nants, :, 2, 2)
        assert base_data_select.data_array.shape == plex_data_select.data_array.shape
        base_data_array_select = np.zeros((len(ants), len(ants), base_data_select.Nspws, base_data_select.Ntimes, base_data_select.Nfreqs, base_data_select.Npols)).astype(np.complex128)
        plex_data_array_select = np.zeros((len(ants), len(ants), plex_data_select.Nspws, plex_data_select.Ntimes, plex_data_select.Nfreqs, plex_data_select.Npols)).astype(np.complex128)
        
        for (i, ant1) in enumerate(ants):
            for (j,ant2) in enumerate(ants):
                for spw in range(base_data_select.Nspws): 
                    if ant1 <= ant2:
                        baseline_number = 2048*(ant1+1)+(ant2+1)+2**16 
                        # baseline index = 2048 * (ant1+1) + (ant2+1) + 2**16
                        base_data_array_copy = copy.deepcopy(base_data_select.get_data(baseline_number))
                        base_flag_array_copy = copy.deepcopy(base_data_select.get_flags(baseline_number))
                        # uvdata.get_data() or uvdata.get_flags() returns an array with a shape (Ntimes, Nfreqs, Npols)
                        plex_data_array_copy = copy.deepcopy(plex_data_select.get_data(baseline_number))
                        plex_flag_array_copy = copy.deepcopy(plex_data_select.get_flags(baseline_number))
                        
                        flag_array_copy = np.logical_or(base_flag_array_copy, plex_flag_array_copy)
                        # zero the flagged data 
                        # the flags are base.flags | plex.flags
                        base_data_array_copy *= ~flag_array_copy
                        plex_flag_array_copy *= ~flag_array_copy
                        
                        base_data_array_select[i,j,spw, :,:,0], base_data_array_select[i,j,spw, :,:,1],  base_data_array_select[i,j,spw, :,:,2],  base_data_array_select[i,j,spw, :,:,3] = base_data_array_copy[:,:,0], base_data_array_copy[:,:,2],  base_data_array_copy[:,:,3],  base_data_array_copy[:,:,1]
                        plex_data_array_select[i,j,spw, :,:,0], plex_data_array_select[i,j,spw, :,:,1], plex_data_array_select[i,j,spw, :,:,2], plex_data_array_select[i,j,spw, :,:,3] = plex_data_array_copy[:,:,0], plex_data_array_copy[:,:,2], plex_data_array_copy[:,:,3], plex_data_array_copy[:,:,1]
                        """
                        data_array orginally stores the pseudo Stokes parameters in the order of a 1d array [-5, -6, -7, -8], or [XX, YY, XY, YX].
                        Here we first reorder it to [-5,-7, -8,-6], 
                        and then reshape it into a 2d array like [[-5,-7],[-8,-6]], or [[XX, XY], [YX, YY]].
                        """
                        
                    else:
                        baseline_number = 2048*(ant2+1)+(ant1+1)+2**16
                        base_data_array_copy = copy.deepcopy(np.conj(base_data_select.get_data(baseline_number)))
                        base_flag_array_copy = copy.deepcopy(base_data_select.get_flags(baseline_number))
                        plex_data_array_copy = copy.deepcopy(np.conj(plex_data_select.get_data(baseline_number)))
                        plex_flag_array_copy = copy.deepcopy(plex_data_select.get_flags(baseline_number))
                        
                        flag_array_copy = np.logical_or(base_flag_array_copy, plex_flag_array_copy)
                        # the flags are base.flags | plex.flags
                        base_data_array_copy *= ~flag_array_copy
                        plex_flag_array_copy *= ~flag_array_copy
                        
                        base_data_array_select[i,j,spw, :,:,0], base_data_array_select[i,j,spw, :,:,1],  base_data_array_select[i,j,spw, :,:,2],  base_data_array_select[i,j,spw, :,:,3] = base_data_array_copy[:,:,0], base_data_array_copy[:,:,3], base_data_array_copy[:,:,2], base_data_array_copy[:,:,1]
                        plex_data_array_select[i,j,spw, :,:,0], plex_data_array_select[i,j,spw, :,:,1], plex_data_array_select[i,j,spw, :,:,2], plex_data_array_select[i,j,spw, :,:,3] = plex_data_array_copy[:,:,0], plex_data_array_copy[:,:,3], plex_data_array_copy[:,:,2], plex_data_array_copy[:,:,1]
                        """
                        Since V_{ji} = V_{ij}^H, we should take conjugate values here and then reorder the order of pols as
                        [-5,-8, -7,-6], which become [[-5,-8],[-7,-6]] after converted into a 2d array.
                        """
                        
        data_shape = base_data_array_select.shape
        base_data_array_select = base_data_array_select.reshape((data_shape[0], data_shape[1], data_shape[2], data_shape[3], data_shape[4], 2, 2))
        plex_data_array_select = plex_data_array_select.reshape((data_shape[0], data_shape[1], data_shape[2], data_shape[3], data_shape[4], 2, 2))
        # reshape the pols
        base_data_array_select = base_data_array_select.reshape((data_shape[0], data_shape[1], data_shape[2]*data_shape[3]*data_shape[4], 2, 2))
        plex_data_array_select = plex_data_array_select.reshape((data_shape[0], data_shape[1], data_shape[2]*data_shape[3]*data_shape[4], 2, 2))
        # concatenate axis to a final shape (Nants, Nants, :, 2, 2)

        self.base_data_array_select,  self.plex_data_array_select = base_data_array_select, plex_data_array_select
        
    def Wirtinger_lm_cal(self, damping_para=0., update_damping_para_per_loop=False, initial_gain_amp=1.0, diagonalize=False, Niteration=50, verbose=False, epsilon=1e-12):
        """
        Using Newton-Gauss method to iteratively obtain gains G which minimizing \sum{D_[ij]-G_i M_{ij} G_j^H}, where D, G and M are all 2*2 matrices. 
        Update each step: G_{k+1} = [J(G_k)^H J(G_k)]^{-1} * J(G_k)^H * D, where J is the Jacobian matrix. 

        Parameters
        ----------
        damping_para : int
            
        initial_gain_amp : int
            Matrix to start with take the form of [[initial_guess_amp, 0],[0, initial_guess_amp]]
            
        Returns
        -------
        gain_array_iter : ndarray
            Iterative gains. 
            
        sum_of_residuals : ndarray
            Iterative residuals.
      
        """
        ants, plex, base = self.ants2cal, self.plex_data_array_select, self.base_data_array_select
        gain_prev = np.array([[initial_gain_amp,0],[0,initial_gain_amp]]).astype(np.complex128)
        gain_prev = np.repeat(gain_prev[np.newaxis,:,:], len(ants), axis=0)
        gain_H_prev = np.copy(gain_prev)
        gain_next, gain_H_next = np.zeros_like(gain_prev), np.zeros_like(gain_prev)
        
        gain_array_iter = np.zeros((Niteration, *gain_prev.shape), dtype=np.complex128)
        sum_of_residuals = np.zeros((Niteration,2,2))

        for iteration in range(Niteration):
            # update gains
            for (i,ant) in enumerate(ants):
                JH_J = np.zeros((2,2)).astype(np.complex128)
                JH_D = np.zeros((2,2)).astype(np.complex128)
                for (j,ant_q) in enumerate(ants):
                    if j!=i:
                        # sum over baselines, frequencies and times
                        JH_J += np.sum(np.matmul(base[i,j], np.matmul(gain_H_prev[j][None], np.matmul(gain_prev[j][None], base[j,i]))), axis=0)
                        JH_D += np.sum(np.matmul(plex[i,j], np.matmul(gain_prev[j][None], base[j, i])), axis=0)
                        
                # For possible singular matrix, we add a very small diagonal matrix onto it to make it invertible.
                if np.isclose(np.linalg.det(JH_J),0):
                    JH_J += np.array([[epsilon,0],[0,epsilon]])
                JH_D_inv_JH_J = np.matmul(JH_D, np.linalg.inv(JH_J))
                if np.isclose(np.linalg.det(JH_D_inv_JH_J),0):
                    gain_next[i] = np.zeros((2,2)).astype(np.complex128)
                else:
                    gain_next[i] = damping_para/(1.+damping_para)*gain_prev[i] + 1./(1.+damping_para)*JH_D_inv_JH_J
                
                if diagonalize==True:
                    gain_next[i] = np.diag(np.diag(gain_next[i]))
                gain_H_next[i] = np.transpose(np.conj(gain_next[i]))
            
            # update residuals                         
            for (i,ant) in enumerate(ants):
                for (j,ant_q) in enumerate(ants):
                    if i != j:
                        sum_of_residuals[iteration] += np.sum(np.abs(plex[i,j] - np.matmul(gain_next[i][None], np.matmul(base[i,j], gain_H_next[j][None])))**2, axis=0)
            
            # update damping parameter if needed
            if update_damping_para_per_loop and iteration >= 1:
                if np.sum(sum_of_residuals[iteration]) >= np.sum(sum_of_residuals[iteration-1]):
                    damping_para *= 1.5
                    gain_prev = (gain_next + gain_prev)/2. 
                    gain_H_prev = (gain_H_next + gain_H_prev)/2.
                else:
                    damping_para /= 1.5
                    gain_prev = gain_next
                    gain_H_prev = gain_H_next
            else: 
                gain_prev = gain_next
                gain_H_prev = gain_H_next
            gain_array_iter[iteration] = gain_prev
        
        # Update results to self.gain_array
        # Rephase to one antenna feed to reduce the degeneracies
        phase_constant = np.angle(gain_prev[0,0,0])   
        gain_prev = gain_prev*np.exp(-1.j*phase_constant)
        gain_H_prev = gain_H_prev*np.exp(1.j*phase_constant)
        
        ants_all = np.sort(np.unique(self.plex_data.ant_1_array))
        #???
        for (i, ant) in enumerate(ants):
            ant_index = np.where(ants_all==ant)
            for freq_index in self.freqs2cal:
                self.gain_array[ant_index, freq_index,:,:] = gain_prev[i,:,:]

        if verbose:
            return gain_array_iter,sum_of_residuals  

def gain_array_interpolate(gain_array, interpolate_freq_index, approach="real_imag", need_smoothing=True):
    """
    Parameters
    ----------
    gain_array : ndarray
        Shape (Nants, Nfreqs, 2, 2)
        
    interpolate_freq_index : ndarray
        1d array containing the frequency indices for interpolation along the axis of frequency
        
    approach : str
        "amp_phase" and "real_imag"
        
    need_smoothing : bool
    
    Returns
    -------
    gain_array_interp : ndarray
    
    interpolate_freq_index : ndarray
    """
    
    assert len(np.shape(gain_array)) == 4, "gain_array must be a 4d array."
    assert np.shape(gain_array)[3] == 2 and np.shape(gain_array)[2] == 2, "gain_array must have a shape as (Nants, Nfreqs, 2, 2)."
    if need_smoothing == True:
        ## get rid of freq_indices where gains are close to zero 
        interpolate_freq_index = interpolate_freq_index[ np.all(abs(np.linalg.det(gain_array[:,interpolate_freq_index,:,:]))>0, axis=0) ]
        ## get rid of large values
        gain_array_freq_avg = np.mean(gain_array[:,interpolate_freq_index,:,:], axis=1)
        interpolate_freq_index = interpolate_freq_index[ np.all(abs(gain_array[:, interpolate_freq_index,:,:]) < 5*abs(gain_array_freq_avg[:,None,:,:]), axis=(0,2,3)) ]
    
    freqs = np.arange(gain_array.shape[1])
    
    # Get amplitudes and phases of the gain_array 
    gain_amp = abs(gain_array)
    gain_phs = np.angle(gain_array)
    gain_phs[np.where(gain_phs<0)] += 2*np.pi
    gain_array_interp = np.zeros_like(gain_array).astype(np.complex128)
    
    for ant in range(gain_array.shape[0]):
        for i in range(2):
            for j in range(2):
                ## get rid of bad channels too big or too small
                if approach == "amp_phase":
                    interpolate_amp = interpolate.interp1d(interpolate_freq_index, gain_amp[ant,interpolate_freq_index,i,j], bounds_error=False, fill_value=[0.])
                    interpolate_phs = interpolate.interp1d(interpolate_freq_index, gain_phs[ant,interpolate_freq_index,i,j], bounds_error=False, fill_value=[0.])
                    gain_array_interp[ant,:,i,j] = interpolate_amp(freqs)*np.exp(1.j*interpolate_phs(freqs))
                if approach == "real_imag":
                    interpolate_real = interpolate.interp1d(interpolate_freq_index, np.real(gain_array[ant,interpolate_freq_index,i,j]), bounds_error=False, fill_value=[0.])
                    interpolate_imag = interpolate.interp1d(interpolate_freq_index, np.imag(gain_array[ant,interpolate_freq_index,i,j]),bounds_error=False, fill_value=[0.])
                    gain_array_interp[ant,:,i,j] = interpolate_real(freqs) + 1.j * interpolate_imag(freqs)
    return gain_array_interp, interpolate_freq_index

def apply_gains(uvd, gain_array_, vis_units="Jy", need_interpolate=False, interpolate_freq_index=None, need_inverse=True, approach="real_imag", need_smoothing=True, epsilon=1e-12):
    """
    uvd_cal = gain_array@uvd@gain_array
    Parameters
    ----------
    uvd : UVData
    
    gain_array_ : ndarray
        Shape (Nants, Nfreqs, 2, 2)
    
    need_interpolate : bool
    
    interpolate_freq_index : ndarray
        1d array containing the frequency indices for interpolation along the axis of frequency
        
    need_inverse : bool
    
    approach : str
        "amp_phase" and "real_imag"
    
    need_smoothing : bool
    
    Returns
    -------
    uvd_cal : UVData
        calibrated result after applying gain_array_ on uvd
    """
    gain_array = copy.deepcopy(gain_array_)
    if need_interpolate:
        if interpolate_freq_index is not None:
            gain_array = gain_array_interpolate(gain_array, interpolate_freq_index, approach=approach, need_smoothing=need_smoothing)
        else:
            raise ValueError("interpolate_freq_index must not be None.")
    if need_inverse:
        for ant_index in range(gain_array.shape[0]):
            for freq_index in range(gain_array.shape[1]):
                # For possible singular matrix, we add a very small diagonal matrix onto it to make it invertible.
                if np.isclose(np.linalg.det(gain_array[ant_index,freq_index]),0):
                    gain_array[ant_index,freq_index] += np.array([[epsilon,0],[0,epsilon]]) 
                gain_array[ant_index,freq_index] = np.linalg.inv(gain_array[ant_index,freq_index])
        
    
    ant_dicts = dict(enumerate(np.sort(np.unique(uvd.ant_1_array))))
    ant_dicts = dict((v,k) for k,v in ant_dicts.items())
    
    ant_1_gain_array = []
    for ant in uvd.ant_1_array[0:uvd.Nbls]:
        ant_index = ant_dicts[ant]
        ant_1_gain_array.append(gain_array[ant_index,:,:,:])
    ant_1_gain_array = np.array(ant_1_gain_array)
    ant_1_gain_array = ant_1_gain_array[:,None,:,:,:]
    ant_1_gain_array = np.tile(ant_1_gain_array, (uvd.Ntimes,1,1,1,1))
    
    ant_2_gain_array = []
    for ant in uvd.ant_2_array[0:uvd.Nbls]:
        ant_index = ant_dicts[ant]
        ant_2_gain_array.append(gain_array[ant_index,:,:,:])
    ant_2_gain_array = np.array(ant_2_gain_array)
    ant_2_gain_array = ant_2_gain_array[:,None,:,:,:]
    ant_2_gain_array = np.tile(ant_2_gain_array,(uvd.Ntimes,1,1,1,1))
    ant_2_gain_array = np.conj(ant_2_gain_array)
    ant_2_gain_array = np.transpose(ant_2_gain_array, (0,1,2,4,3))
    
    raw_data_array = np.zeros_like(uvd.data_array).astype(np.complex128)
    raw_data_array[:,:,:,0], raw_data_array[:,:,:,1],  raw_data_array[:,:,:,2],  raw_data_array[:,:,:,3] = uvd.data_array[:,:,:,0], uvd.data_array[:,:,:,2], uvd.data_array[:,:,:,3], uvd.data_array[:,:,:,1]
    raw_data_array = raw_data_array.reshape((raw_data_array.shape[0],raw_data_array.shape[1], raw_data_array.shape[2],2,2))
    
    cal_data_array = np.matmul(ant_1_gain_array, np.matmul(raw_data_array, ant_2_gain_array))
    uvd_cal = copy.deepcopy(uvd)
    uvd_cal.data_array[:,:,:,0], uvd_cal.data_array[:,:,:,1], uvd_cal.data_array[:,:,:,2], uvd_cal.data_array[:,:,:,3] = cal_data_array[:,:,:,0,0], cal_data_array[:,:,:,1,1], cal_data_array[:,:,:,0,1], cal_data_array[:,:,:,1,0]
    uvd_cal.vis_units = vis_units
    return uvd_cal
        
















