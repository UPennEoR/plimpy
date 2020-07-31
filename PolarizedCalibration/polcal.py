import numpy as np
import matplotlib.pyplot as plt
from pyuvdata import UVData
import scipy
import copy
import time
import glob
from scipy import interpolate

class uvdata_pol_calibrator():
    
    def __init__(self, model_data=None, real_data=None, mode="model_based"):
        """
        Parameters
        ----------
        model_data, real_data : UVData
        
        mode : str
            Specify the way gains are aligned. Choices are "model_based" and "data_based".
            "model_based" means V_ij_data = G_i V_ij_model G_j^H,
            and "data_based" means V_ij_model = G_i V_ij_data G_j^H.
            
        """
        if mode == "model_based":
            self.base_data, self.prod_data = model_data, real_data
        if mode == "data_based":
            self.base_data, self.prod_data = real_data, model_data
        self.gain_array = np.zeros((self.base_data.Nants_data, self.base_data.Nfreqs, 2, 2)).astype(np.complex128)
        self.gain_H_array = np.zeros((self.base_data.Nants_data,self.base_data.Nfreqs, 2, 2)).astype(np.complex128)
        
    def data_slice(self, use_all_times=False, base_time_range=[], prod_time_range=[], use_all_frequencies=False, base_freq_range=[], prod_freq_range=[], 
                   use_all_ants=False, base_ants=[], prod_ants=[]):
        """
        Get a slice of the data.
        
        Parameters
        ----------

  
        """
        if use_all_ants:
            ants = np.unique(self.prod_data.ant_1_array)
        else:
            assert isinstance(base_ants, list) and isinstance(prod_ants, list), "ants must be a list." 
            base_ants, prod_ants= list(np.sort(base_ants)), list(np.sort(prod_ants))
            assert base_ants == prod_ants
            ants = prod_ants
        self.ants2cal = ants
        
        if use_all_times:
            base_time_range = [np.unique(self.base_data.time_array)[0]-1e-2, np.unique(self.base_data.time_array)[-1]+1e-2]
            prod_time_range = [np.unique(self.prod_data.time_array)[0]-1e-2, np.unique(self.prod_data.time_array)[-1]+1e-2]
        else:
            assert isinstance(base_time_range, list) and isinstance(prod_time_range, list), "time_range must be a list"
            assert len(base_time_range) == len(prod_time_range), "Length of time_ranges must be the same."
            assert len(base_time_range) == 2, "Length of time_range must be 2."
        
        if use_all_frequencies:
            base_freqs = list(range(self.base_data.Nfreqs))
            prod_freqs = list(range(self.prod_data.Nfreqs))
        else:
            assert isinstance(base_freq_range, list) and isinstance(prod_freq_range, list), "freq_range must be a list"
            assert base_freq_range == prod_freq_range, "freq_ranges must be the same."
            if len(base_freq_range) == 1:
                base_freqs = base_freq_range
                prod_freqs = prod_freq_range
            elif len(base_freq_range) == 2:
                base_freqs = list(range(base_freq_range[0],base_freq_range[1]+1))
                prod_freqs = list(range(prod_freq_range[0],prod_freq_range[1]+1))
            else:
                raise ValueError("Length of freq_range must be 1 or 2.")
        self.freqs2cal = prod_freqs
        
        base_data_slice = self.base_data.select(antenna_nums=ants, freq_chans=base_freqs, time_range=base_time_range, inplace=False) 
        prod_data_slice = self.prod_data.select(antenna_nums=ants, freq_chans=prod_freqs, time_range=prod_time_range, inplace=False)             
        
        # original data_array shape: (Nblts, Nspws, Nfreqs, )
        # intermidiate data_array shape: (Nants, Nants, Nspws, Ntimes, Nfreqs, Npols)
        # final data_array shape: (Nants, Nants, :, 2, 2)
        assert base_data_slice.data_array.shape == prod_data_slice.data_array.shape
        base_data_array = np.zeros((len(ants), len(ants), base_data_slice.Nspws, base_data_slice.Ntimes, base_data_slice.Nfreqs, base_data_slice.Npols)).astype(np.complex128)
        prod_data_array = np.zeros((len(ants), len(ants), prod_data_slice.Nspws, prod_data_slice.Ntimes, prod_data_slice.Nfreqs, prod_data_slice.Npols)).astype(np.complex128)
        
        for (i, ant1) in enumerate(ants):
            for (j,ant2) in enumerate(ants):
                for spw in range(base_data_slice.Nspws): 
                    if ant1 <= ant2:
                        baseline_number = 2048*(ant1+1)+(ant2+1)+2**16
                        # baseline index = 2048 * (ant1+1) + (ant2+1) + 2**16
                        base_data_array_copy = base_data_slice.get_data(baseline_number)
                        # uvdata.get_data() returns an array with a shape (Ntimes, Nfreqs, Npols) 
                        base_data_array[i,j,spw, :,:,0], base_data_array[i,j,spw, :,:,1],  base_data_array[i,j,spw, :,:,2],  base_data_array[i,j,spw, :,:,3] = base_data_array_copy[:,:,0], base_data_array_copy[:,:,2],  base_data_array_copy[:,:,3],  base_data_array_copy[:,:,1]
                        """
                        data_array orginally stores the polarization information as a 1d array [-5, -6, -7, -8], corresponding to [XX, YY, XY, YX].
                        Here we first modify it into [-5,-7, -8,-6], 
                        then we will modify it into a 2d array[[-5,-7],[-8,-6]], corresponding to [[XX, XY], [YX, YY]].
                        """
                        prod_data_array_copy = prod_data_slice.get_data(baseline_number)
                        prod_data_array[i,j,spw, :,:,0], prod_data_array[i,j,spw, :,:,1], prod_data_array[i,j,spw, :,:,2], prod_data_array[i,j,spw, :,:,3] = prod_data_array_copy[:,:,0], prod_data_array_copy[:,:,2], prod_data_array_copy[:,:,3], prod_data_array_copy[:,:,1]

                    if ant1 > ant2:
                        baseline_number = 2048*(ant2+1)+(ant1+1)+2**16
                        base_data_array_copy = np.conj(base_data_slice.get_data(baseline_number))
                        base_data_array[i,j,spw, :,:,0], base_data_array[i,j,spw, :,:,1],  base_data_array[i,j,spw, :,:,2],  base_data_array[i,j,spw, :,:,3] = base_data_array_copy[:,:,0], base_data_array_copy[:,:,3], base_data_array_copy[:,:,2], base_data_array_copy[:,:,1]
                        """
                        Since V_{ji} = V_{ij}^H, we should take conjugate values here and then reorder the pols as
                        [-5,-8, -7,-6], which becomes [[-5,-8],[-7,-6]] after converted into a 2d array.
                        """
                        prod_data_array_copy = np.conj(prod_data_slice.get_data(baseline_number))
                        prod_data_array[i,j,spw, :,:,0], prod_data_array[i,j,spw, :,:,1], prod_data_array[i,j,spw, :,:,2], prod_data_array[i,j,spw, :,:,3] = prod_data_array_copy[:,:,0], prod_data_array_copy[:,:,3], prod_data_array_copy[:,:,2], prod_data_array_copy[:,:,1]


        data_shape = base_data_array.shape
        base_data_array = base_data_array.reshape((data_shape[0], data_shape[1], data_shape[2], data_shape[3], data_shape[4], 2, 2))
        prod_data_array = prod_data_array.reshape((data_shape[0], data_shape[1], data_shape[2], data_shape[3], data_shape[4], 2, 2))
        # reshape the pols
        base_data_array = base_data_array.reshape((data_shape[0], data_shape[1], data_shape[2]*data_shape[3]*data_shape[4], 2, 2))
        prod_data_array = prod_data_array.reshape((data_shape[0], data_shape[1], data_shape[2]*data_shape[3]*data_shape[4], 2, 2))
        # concatenate axis to a shape (Nants, Nants, :, 2, 2)

        self.base_data_array,  self.prod_data_array = base_data_array, prod_data_array
        
    def Wirtinger_lm_cal(self, diagonalize=False, Niteration=50, including_autobaseline=False, verbose=False, epsilon=1e-12):
        """
        Using Newton-Gauss method to obtain calibration gains G which minimizing \sum{D_[ij]-G_i M_{ij} G_j^H}, where D, G and M are all 2*2 matrices. 
        Update each step: G_{k+1} = [J(G_k)^H J(G_k)]^{-1} * J(G_k)^H * D, where J is the Jacobian matrix. 

        Parameters
        ----------

       

        Returns
        -------

        gain : dict
            calibration gains

        residual : 

        """
        ants, prod, base = self.ants2cal, self.prod_data_array, self.base_data_array
        gain_prev = np.array([[1,0],[0,1]]).astype(np.complex128)
        gain_prev = np.repeat(gain_prev[np.newaxis,:,:], len(ants), axis=0)
        gain_H_prev = np.copy(gain_prev)
        gain_next, gain_H_next = np.zeros_like(gain_prev), np.zeros_like(gain_prev)
        
        print('base', base.shape)
        print('prod', prod.shape)
        print('freqs2cal', self.freqs2cal)
        
        tmp = list(self.gain_array.shape)
        print('gain_array shape', tmp)
        tmp.append(Niteration)
        print('wtf', tmp)
        gain_array_iter = np.zeros(tmp, dtype=np.complex128)
        print('gain_array_iter shape', gain_array_iter.shape)
        
        residual = np.zeros(Niteration)

        for iteration in range(Niteration):
            for (i,ant) in enumerate(ants):
                JH_J = np.zeros((2,2)).astype(np.complex128)
                JH_D = np.zeros((2,2)).astype(np.complex128)
                for (j,ant_q) in enumerate(ants):
                    # sum over baselines, frequencies and times
                    JH_J += np.sum(np.matmul(base[i,j], np.matmul(gain_H_prev[j][None], np.matmul(gain_prev[j][None], base[j,i]))), axis=0)
                    JH_D += np.sum(np.matmul(prod[i,j], np.matmul(gain_prev[j][None], base[j, i])), axis=0)
                if not including_autobaseline:
                    # if not including auto-baseline
                    JH_J -= np.sum(np.matmul(base[i,i], np.matmul(gain_H_prev[i][None], np.matmul(gain_prev[i][None], base[i,i]))), axis=0)
                    JH_D -= np.sum(np.matmul(prod[i,i], np.matmul(gain_prev[i][None], base[i, i])), axis=0)         
                
                # For possible singular matrix, we add a very small diagonal matrix onto it to make it invertible.
                if np.isclose(np.linalg.det(JH_J),0):
                    JH_J += np.array([[epsilon,0],[0,epsilon]])   
                if diagonalize==True:
                    gain_next[i] = np.diag(np.diag(np.matmul(JH_D, np.linalg.inv(JH_J))))
                else:
                    gain_next[i] = np.matmul(JH_D, np.linalg.inv(JH_J))
                
                JH_J = np.zeros((2,2)).astype(np.complex128)
                JH_D = np.zeros((2,2)).astype(np.complex128)
                for (k,ant_p) in enumerate(ants):
                    JH_J += np.sum(np.matmul(base[i,k], np.matmul(gain_H_prev[k][None], np.matmul(gain_prev[k][None], base[k,i]))), axis=0)
                    JH_D += np.sum(np.matmul(base[i,k], np.matmul(gain_H_prev[k][None], prod[k, i])), axis=0)     
                if not including_autobaseline:
                    JH_J -= np.sum(np.matmul(base[i,i], np.matmul(gain_H_prev[i][None], np.matmul(gain_prev[i][None], base[i,i]))), axis=0)
                    JH_D -= np.sum(np.matmul(base[i,i], np.matmul(gain_H_prev[i][None], prod[i, i])), axis=0)  
                
                if np.isclose(np.linalg.det(JH_J),0):
                    JH_J += np.array([[epsilon,0],[0,epsilon]])
                if diagonalize==True:
                    gain_H_next[i] = np.diag(np.diag(np.matmul(np.linalg.inv(JH_J), JH_D)))
                else:
                    gain_H_next[i] = np.matmul(np.linalg.inv(JH_J), JH_D)
            
            # Rephase to one antenna feed to reduce the degeneracies
            gain_prev = gain_next*np.exp(-1.j*np.angle(gain_next[0,0,0]))
            gain_H_prev = gain_H_next*np.exp(1.j*np.angle(gain_next[0,0,0]))

            for ant in range(len(ants)):
                for ant_r in range(len(ants)):
                    if including_autobaseline:
                         residual[iteration] += np.linalg.norm(prod[ant,ant_r] - np.matmul(gain_prev[ant][None], np.matmul(base[ant,ant_r], gain_H_prev[ant_r][None])))
                    else:
                        if ant_r != ant:
                            residual[iteration] += np.linalg.norm(prod[ant,ant_r] - np.matmul(gain_prev[ant][None], np.matmul(base[ant,ant_r], gain_H_prev[ant_r][None])))
            
            print('gain_prev shape', gain_prev.shape)
            for freq_index in self.freqs2cal:
                gain_array_iter[:,freq_index,:,:,iteration] = gain_prev
        
        # Add results to self.gain_array
        ants_all = np.unique(self.prod_data.ant_1_array)
        for i, ant in enumerate(ants):
            ant_index = np.where(ants_all==ant)
            for freq_index in self.freqs2cal:
                self.gain_array[ant_index, freq_index,:,:] = gain_prev[i,:,:]
                self.gain_H_array[ant_index, freq_index,:,:] = gain_H_prev[i,:,:]

        if verbose:
            print(residual)                                         
                              
        return gain_array_iter


def gain_array_interpolate(gain_array, interpolate_freq_index):
    """
    Parameters
    ----------
    
    """
    # Get amplitudes and phases of the gain_array 
    gain_amp = abs(gain_array)
    gain_phs = np.angle(gain_array)
    gain_phs[np.where(gain_phs<0)] += 2*np.pi
    gain_array_interp = np.zeros_like(gain_array).astype(np.complex128)
    freqs = np.arange(gain_array.shape[1])
    for ant in range(gain_array.shape[0]):
        for i in range(2):
            for j in range(2):
                interpolate_amp = interpolate.interp1d(interpolate_freq_index, gain_amp[ant,interpolate_freq_index,i,j],fill_value="extrapolate")
                interpolate_phs = interpolate.interp1d(interpolate_freq_index, gain_phs[ant,interpolate_freq_index,i,j],fill_value="extrapolate")
                gain_array_interp[ant,:,i,j] = interpolate_amp(freqs)*np.exp(1.j*interpolate_phs(freqs))
    return gain_array_interp


def apply_gains(uvd, gain_array_, need_interpolate=False, interpolate_freq_index=None, need_inverse=True, epsilon=1e-12):
    """
    Parameters
    ----------
    
    """
    gain_array = copy.deepcopy(gain_array_)
    if need_interpolate:
        if interpolate_freq_index is not None:
            gain_array = gain_array_interpolate(gain_array, interpolate_freq_index)
    if need_inverse:
        for ant_index in range(gain_array.shape[0]):
            for freq_index in range(gain_array.shape[0]):
                # For possible singular matrix, we add a very small diagonal matrix onto it to make it invertible.
                if np.isclose(np.linalg.det(gain_array[ant_index,freq_index]),0):
                    gain_array[ant_index,freq_index] += np.array([[epsilon,0],[0,epsilon]]) 
                gain_array[ant_index,freq_index] = np.linalg.inv(gain_array[ant_index,freq_index])
        
    
    ant_dicts = dict(enumerate(np.unique(uvd.ant_1_array)))
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
    
    return uvd_cal


