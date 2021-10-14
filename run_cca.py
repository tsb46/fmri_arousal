import argparse
import numpy as np
import fbpca
import pickle

from run_pca import pca
from scipy.signal import hilbert
from scipy.stats import zscore
from utils.butterworth_filters import filter_functional_data
from utils.load_write import load_data, write_nifti
from utils.kernel_cca import KernelCCA


def calculate_components_time_delay(ts, n_ts, window, weights):
    n_comp = weights.shape[1]
    indx=0
    comps = []
    if window > 0:
        for t in range(window*2+1):
            x_t = ts[:, indx:(indx+n_ts)]
            x_c = np.dot(x_t.T, weights)
            comps.append(x_c)
            indx+= n_ts
        return np.stack(comps, axis=2)
    else:
        return np.dot(ts.T, weights)


def compute_kernel(ts, window, time_delay=False):
    if time_delay:
        ts = time_delay_matrix(ts, window)
    ts = zscore(ts)
    k = np.dot(ts, ts.T)
    return k, ts


def time_delay_matrix(ts, window):
    N, P = ts.shape
    if window > 0:
        window_arr = np.arange(-window, window+1)
        time_delay_mat = [np.pad(ts, ((p,0), (0,0)), mode='constant')[:-p,:] if p > 0 
                          else np.pad(ts, ((0,np.abs(p)), (0,0)), mode='constant')[np.abs(p):,:]
                          for p in window_arr]
        return np.concatenate(time_delay_mat, axis=1)
    else:
        return ts


def write_results(level, mask, kcca_obj, cca_type, comp_weights, 
                  subj_n, scan, zero_mask, n_vert):
    if level == 'group':
        analysis_str = 'cca_group'
    elif level == 'subject':
        analysis_str = f'cca_s{subj_n}'
    if scan is not None:
        analysis_str += f'_{scan}'

    # Write nifti
    if cca_type == 'tkcca': 
        analysis_str += '_tkcca'
        write_nifti(comp_weights, mask, f'{analysis_str}', zero_mask, n_vert)        
    elif pca_type == 'mca':
        analysis_str += '_mca'
        write_nifti(comp_weights, mask, f'{analysis_str}', zero_mask, n_vert)
    pickle.dump(kcca_obj, open(f'{analysis_str}_results.pkl', 'wb'))


def run_main(n_comps, level, subj_n, scan_n, low_cut, high_cut,
             cca_type, delay_choice, dim_reduce, dim_n, reg, window, 
             physio_params_fp, mask):
    if level == 'subject' and subj_n is None:
        raise Exception('Subject number must be provided for subject-level analysis')

    func_data, eeg, rv, hr, zero_mask, n_vert = load_data(level, mask, physio_params_fp, subj_n, scan_n)
    # Filter functional data if specified
    func_data = filter_functional_data(func_data, low_cut, high_cut, physio_params_fp)
    if dim_reduce:
        func_data_orig = zscore(func_data)
        pca_res = pca(func_data_orig, dim_n)
        func_data = pca_res['pc_scores']

    physio_data = np.stack([eeg, rv, hr], axis=1)
    # Normalize data
    if cca_type == 'tkcca':
        kcca = KernelCCA(n_components=n_comps, kernel='precomputed', kapa=reg)
        if delay_choice == 'physio':
            k_func, func_data = compute_kernel(func_data, window)
            k_physio, physio_data_delay = compute_kernel(physio_data, window, time_delay=True)
            kcca.fit(k_physio, k_func)
        elif delay_choice == 'func':
            k_physio, physio_data = compute_kernel(physio_data, window)
            k_func, func_data_delay = compute_kernel(func_data, window, time_delay=True)
            kcca.fit(k_func,k_physio)
        if delay_choice == 'physio':
            physio_weights = calculate_components_time_delay(physio_data_delay, 
                                                             physio_data.shape[1], 
                                                             window, kcca.alphas_)
            if dim_reduce:
                func_weights = np.dot(func_data_orig.T, kcca.betas_).T
            else:
                func_weights = np.dot(func_data.T, kcca.betas_).T
        elif delay_choice == 'func':
            func_weights = calculate_components_time_delay(func_data_delay, 
                                                           func_data.shape[1], 
                                                           window, kcca.alphas_)
            if dim_reduce:
                func_weights = np.stack([np.squeeze(func_weights[:,i,:]), for i in range(n_comps)]


            physio_weights = np.dot(physio_data.T, kcca.betas_)

        kcca.physio_weights = physio_weights

    write_results(level, mask, kcca, cca_type, 
                  func_weights, subj_n, scan_n,
                  zero_mask, n_vert)


if __name__ == '__main__':
    """Run main analysis"""
    parser = argparse.ArgumentParser(description='Run CCA analysis')
    parser.add_argument('-n', '--n_comps',
                        help='<Required> Number of components from CCA',
                        required=True,
                        type=int)
    parser.add_argument('-l', '--level',
                        help='subject or group level analysis',
                        default='group',
                        choices=['subject', 'group'],
                        type=str)
    parser.add_argument('-s', '--subject_n',
                        help='subject number for subject level analysis',
                        default=None,
                        type=int)
    parser.add_argument('-scan', '--scan_n',
                        help='scan number for subject level analysis (if multiple runs from same subject',
                        default=None,
                        type=int)
    parser.add_argument('-fl', '--low_cut', 
                        help='lowcut for butterworth filter for functional data',
                        default=None,
                        type=float) 
    parser.add_argument('-fh', '--high_cut', 
                        help='highcut for butterworth filter for functional data',
                        default=None,
                        type=float) 
    parser.add_argument('-t', '--cca_type',
                        help='Calculate complex or real PCA',
                        default='tkcca',
                        choices=['tkcca', 'complex'],
                        type=str)
    parser.add_argument('-d', '--dim_reduce', 
                        help='Whether to run dimension reduction on functional data (PCA) before CCA',
                        default=1,
                        type=int)
    parser.add_argument('-dn', '--dim_n', 
                        help='Number of components for dimension reduction',
                        default=100,
                        type=int)
    parser.add_argument('-r', '--regularization',
                        help='regularization parameter for tkCCA',
                        default=0,
                        type=float)
    parser.add_argument('-x', '--delay_choice', 
                        help='whether to time-delay physio signals or functional data in tkCCA',
                        default='physio',
                        choices=['physio', 'func'],
                        type=str)
    parser.add_argument('-w', '--window',
                        help='Size of left and right window for tkCCA - Must be even!!!',
                        default=10,
                        required=False,
                        type=int)
    parser.add_argument('-p', '--physio_params', 
                        help='file path to preprocessing params for physio signals',
                        default='physio_proc_params.json',
                        type=str)
    parser.add_argument('-m', '--mask', 
                        help='file path to brain mask',
                        default='masks/MNI152_T1_3mm_brain_mask.nii.gz',
                        type=str)
    args_dict = vars(parser.parse_args())
    run_main(args_dict['n_comps'], args_dict['level'], 
             args_dict['subject_n'], args_dict['scan_n'],
             args_dict['low_cut'], args_dict['high_cut'],
             args_dict['cca_type'], args_dict['delay_choice'], 
             args_dict['dim_reduce'], args_dict['dim_n'], 
             args_dict['regularization'], args_dict['window'], 
             args_dict['physio_params'], args_dict['mask'])
