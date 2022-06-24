import json
import nibabel as nb 
import numpy as np
import pandas as pd
import os

from scipy.io import loadmat 
from scipy.stats import zscore
from sklearn.linear_model import LinearRegression
from utils.signal.butterworth_filters import filter_functional_data
from utils.load_utils import find_fps, print_filter_info
from utils.load_physio import preprocess_physio

# path to the analysis parameter .json file. Never move this out of the base directory!
params_fp='analysis_params.json'
# path to 3mm dilated brain mask file. Never move this out of the masks directory!
mask_fp = 'masks/MNI152_T1_3mm_brain_mask_dilated.nii.gz'


def convert_2d(mask, nifti_data):
    nonzero_indx = np.nonzero(mask)
    nifti_2d = nifti_data[nonzero_indx]
    return nifti_2d.T


def filter_zero_voxels(nifti_data, group_method, use_first=True):
    # use_first = whether to use first subject as estimate of nan voxels, greatly speeds concatenation
    # Construct mask of vertices with std greater than 0
    if group_method == 'stack':
        orig_n_vert = nifti_data.shape[1]
        zero_mask = ~np.any(np.isnan(nifti_data), axis=0)
        zero_mask_indx = np.where(zero_mask)[0]
        return nifti_data[:, zero_mask], zero_mask_indx, orig_n_vert
    elif group_method == 'list':
        orig_n_vert = nifti_data[0].shape[1]
        if use_first:
            zero_mask = np.std(nifti_data[0], axis=0) > 0
        else: 
            std_vec = np.zeros(orig_n_vert)     
            for n in nifti_data:
                std_vec += (np.std(n, axis=0) > 0)
            zero_mask = (std_vec == len(nifti_data))
        zero_mask_indx = np.where(zero_mask)[0]
        for i, n in enumerate(nifti_data):
            nifti_data[i] = n[:, zero_mask]
        return nifti_data, zero_mask, orig_n_vert


def impute_zero_voxels(nifti_data, zero_mask, orig_n_vert):
    nifti_full = np.zeros([nifti_data.shape[0], orig_n_vert])
    nifti_full[:, zero_mask] = nifti_data
    return nifti_full


def initialize_group_func_array(fps, nscans, mask_n):
    n_t = 0
    for fp in fps:
        nifti = nb.load(fp)
        n_t += nifti.header['dim'][4]
    group_func_array = np.zeros((n_t, mask_n))
    return group_func_array


def load_chang_bh_event_file():
    # We are ASSUMING that the event timings are the same across all subjects 
    events = np.loadtxt(f'data/dataset_chang_bh/adb_onsets.txt')
    return events


def load_data(data, level, physio, load_physio, subj_n=None, scan_n=None, 
              events=False, group_method='stack', physio_group_method='stack', 
              verbose=True, filter_nan_voxels=True, regress_global=False):
    params = load_params()

    # Pull physio labels (if not already selected)
    if data == 'hcp_fix':
        data_str = 'hcp'
    else:
        data_str = data
    
    params_data = params[data_str]
    if physio is None:
        physio = params_data['physio']

    # Load mask
    mask = nb.load(mask_fp).get_fdata() > 0


    # Pull file paths
    fps = find_fps(data, level, physio, params_data, events, subj_n, scan_n)
    # Print filter parameters for functional and physio signals
    if verbose:
        print_filter_info(params_data, load_physio)

    # Pull events, if specified in options
    if events:
        event_df = load_events(fps['events'], level, data, group_method)
        params_data['events'] = event_df

    # Pull data for subject level analysis
    if level == 'subject':
        if verbose:
            print(f'Loading subject: {subj_n}')
        func_data = load_subject_func(fps['func'][0], mask, params_data, regress_global)
        # If load subject and group method = list, place in list for downstream compatibility
        if group_method == 'list':
            func_data = [func_data]
        # If load physio option is specified, load and preprocess physio
        if load_physio:
            physio_proc = []
            for p in physio:
                p_proc = load_subject_physio(fps[p][0], params_data, data, p)
                if physio_group_method == 'list':
                    p_proc = [p_proc]
                physio_proc.append(p_proc)
        else:
            physio_proc = None
    # Pull data for group level analysis
    elif level == 'group':
        if verbose:
            print('Loading all subjects')
            if regress_global:
                print('regressing out global signal')
        func_data = load_group_func(fps['func'], mask, params_data, group_method, regress_global, verbose)
        if load_physio:
            physio_proc = []
            for p in physio:
                p_proc = load_group_physio(fps[p], params_data, data, p, physio_group_method)
                physio_proc.append(p_proc)
        else:
            physio_proc = None

    # Filter voxels with no recorded BOLD signal (i.e. time series of all 0s)
    if filter_nan_voxels:
        func_data, zero_mask, n_vert_orig = filter_zero_voxels(func_data, group_method)
    else:
        n_vert_orig = np.nonzero(mask)[0].shape[0]
        zero_mask = np.tile(1, n_vert_orig).astype(bool)

    return func_data, physio_proc, physio, zero_mask, n_vert_orig, params_data


def load_events_fp(fp, data):
    return 


def load_events(fps, level, data, group_method):
    if level == 'subject':
        event_df = load_events_fp(fp, data)
        if group_method == 'list':
            event_df = [event_df]
    elif level == 'group':
        group_events = []
        for fp in fps:
            events = load_events_fp(fp, data)
            group_events.append(events)
        if group_method == 'stack':
            return pd.concat(group_events, axis=0)
        else:
            return group_events


def load_group_func(fps, mask, params, group_method, regress_global, verbose):
    if group_method == 'stack':
        mask_n = len(np.nonzero(mask)[0])
        indx=0
        group_data = initialize_group_func_array(fps, params['nscans'], mask_n) 
    elif group_method == 'list':
        group_data = []
    # Loop through files and concatenate/append
    for fp in fps:
        if verbose:
            print(fp)
        subj_data = load_subject_func(fp, mask, params, regress_global)
        subj_t = subj_data.shape[0]
        # Normalize data before concatenation
        subj_data = zscore(subj_data)
        if group_method == 'stack':
            group_data[indx:(indx+subj_t), :] = subj_data
            indx += subj_t
        elif group_method == 'list':
            group_data.append(subj_data)

    return group_data



def load_group_physio(fps, params, data, physio_label, group_method):
    physio_all = []
    for fp in fps:
        proc_sig = load_subject_physio(fp, params, data, physio_label)

        if params['data']['physio']['concat_method'] == 'zscore':
            proc_sig = zscore(proc_sig)

        physio_all.append(proc_sig)

    if group_method == 'stack':
        return np.concatenate(physio_all)
    elif group_method == 'list':
        return physio_all


def load_nki_event_file():
    # We are ASSUMING that the event timings are the same across all subjects (e.g. no counterbalancing)
    events = pd.read_csv('data/dataset_nki/events/A00057406_task_breathhold_events.tsv', sep='\t')
    return events


def load_params():
    # Load analysis parameters - filepath is defined as global variable at top of script
    with open(params_fp) as params_file:
        params = json.load(params_file)
    return params


def load_subject_func(fp, mask, params, regress_global):
    # Load scan
    nifti = nb.load(fp, keep_file_open = True)
    nifti_data = nifti.get_fdata()
    nifti.uncache()
    nifti_data = convert_2d(mask, nifti_data)
    nifti_data = filter_functional_data(nifti_data, params)
    if regress_global:
        nifti_data = regress_global_signal(nifti_data)
    return nifti_data


def load_subject_physio(fp, params, data, physio_label):    
    # Load AND preprocess physio
    physio_signal = np.loadtxt(fp)
    if physio_signal.ndim == 1:
        physio_signal = physio_signal[:, np.newaxis]
    physio_signal_proc = preprocess_physio(physio_signal, params, physio_label)
    return physio_signal_proc


def regress_global_signal(func_data, mask_nan =True):
    if mask_nan:
        nan_mask = ~np.any(np.isnan(func_data), axis=0)
        func_data_m = func_data[:, nan_mask]
    else:
        nan_mask = np.arange(func_data.shape[1])
        func_data_m = func_data[:, nan_mask]
    global_signal = func_data_m.mean(axis=1)
    lin_reg = LinearRegression()
    lin_reg.fit(global_signal.reshape(-1, 1), func_data_m)
    func_pred = lin_reg.predict(global_signal.reshape(-1,1))
    func_residual = func_data_m - func_pred
    func_data[:, nan_mask] = func_residual
    return func_data


def write_nifti(data, output_file, zero_mask, orig_n_vert):
    data_imp = impute_zero_voxels(data, zero_mask, orig_n_vert)
    mask = nb.load(mask_fp)
    mask_bin = mask.get_fdata() > 0
    nifti_4d = np.zeros(mask.shape + (data_imp.shape[0],), 
                        dtype=data_imp.dtype)
    nifti_4d[mask_bin, :] = data_imp.T

    nifti_out = nb.Nifti2Image(nifti_4d, mask.affine)
    nb.save(nifti_out, output_file)




    

