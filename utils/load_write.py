import json
import nibabel as nb 
import numpy as np
import pandas as pd
import os

from scipy.io import loadmat 
from scipy.stats import zscore
from sklearn.linear_model import LinearRegression
from utils.signal_utils import butterworth_filter

# path to the analysis parameter .json file. Never move this out of the base directory!
params_fp='analysis_params.json'

# Dilated mask that includes sinuses slightly outside gray matter tissue
mask_3mm="masks/MNI152_T1_3mm_brain_mask_dilated.nii.gz"


def convert_2d(mask_bin, nifti_data):
    # convert 4d nifti to 2d matrix using mask
    nonzero_indx = np.nonzero(mask_bin)
    nifti_2d = nifti_data[nonzero_indx]
    return nifti_2d.T


def convert_4d(mask_bin, nifti_data):
    # convert 2d matrix to 4d nifti using mask
    nifti_4d = np.zeros(mask_bin.shape + (nifti_data.shape[0],), 
                        dtype=nifti_data.dtype)
    nifti_4d[mask_bin, :] = nifti_data.T
    return nifti_4d


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


def find_fps(dataset, physio, params):
    # find file paths to functional and physio files
    # load subject list
    subj_list, scan_list = load_subject_list(dataset, params['subject_list'])
    # get functional file paths
    func_fps = [params['func'].format(subj, scan) 
                for subj, scan in zip(subj_list, scan_list)]
    fps = {
        'func': func_fps
    }
    # if physio is specified, pull physio
    if physio is not None:
        physio_fps = [params['physio'].format(subj, scan, physio) if scan is not None 
                      else params['physio'].format(subj, physio)
                      for subj, scan in zip(subj_list, scan_list)]
        fps['physio'] = physio_fps
    return fps


def get_fp_base(fp):
    # get nifti file path without extenstion
    fp_split = os.path.splitext(fp)
    if fp_split[1] == '.gz':
        fp_base = os.path.splitext(fp_split[0])[0]
    else:
        fp_base = fp_split[0]
    return fp_base


def impute_zero_voxels(nifti_data, zero_mask, orig_n_vert):
    nifti_full = np.zeros([nifti_data.shape[0], orig_n_vert])
    nifti_full[:, zero_mask] = nifti_data
    return nifti_full


def initialize_group_func(fps, nscans, mask_n):
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


def load_chang_cue_event_file(subj, scan):
    # load cue onset and reaction time for a scan in reaction time task
    task_mat = loadmat(
        f'data/dataset_chang_cue/events/sub_{subj}-mr_{scan}-ect_echo1_taskOUT.mat',
        squeeze_me=True
    )
    cue = task_mat['OUT']['stimTime_msec'].item()
    rt = task_mat['OUT']['RT_msec'].item()
    return cue, rt


def load_data(dataset, physio, group_method='stack', 
              physio_group_method='stack', 
              regress_global=False):
    # master function for loading and concatenating functional/physio files
    params = load_params()
    params_d = params[dataset]

    # Load mask
    mask_bin = nb.load(params_d['mask']).get_fdata() > 0

    # Pull file paths
    fps = find_fps(dataset, physio, params_d)

    # Pull data for group level analysis
    func_data = load_group_func(fps['func'], mask_bin, group_method, 
                                params_d, regress_global)
    if physio is not None:
        physio_sig = load_group_physio(fps['physio'], physio_group_method)
    else:
        physio_sig = None

    # Filter voxels with no recorded BOLD signal (i.e. time series of all 0s)
    func_data, zero_mask, n_vert_orig = filter_zero_voxels(func_data, group_method)
    # collect masks for future writing out of results
    params_d['zero_mask'] = zero_mask # non NaN voxel indices
    params_d['n_vert_orig'] = n_vert_orig # number of voxels before masking

    return func_data, physio_sig, params_d


def load_group_func(fps, mask_bin, group_method, params, regress_global):
    if group_method == 'stack':
        mask_n = len(np.nonzero(mask_bin)[0])
        indx=0
        group_data = initialize_group_func(fps, params['nscans'], mask_n) 
    elif group_method == 'list':
        group_data = []
    # Loop through files and concatenate/append
    for fp in fps:
        subj_data = load_subject_func(fp, mask_bin, regress_global)
        subj_t = subj_data.shape[0]
        # Normalize data before concatenation
        subj_data = zscore(subj_data, nan_policy='omit')
        # fill nans w/ 0 in regions of poor functional scan coverage
        subj_data = np.nan_to_num(subj_data)
        if group_method == 'stack':
            group_data[indx:(indx+subj_t), :] = subj_data
            indx += subj_t
        elif group_method == 'list':
            group_data.append(subj_data)

    return group_data


def load_group_physio(fps, group_method):
    physio_all = []
    for fp in fps:
        # load subject physio signal
        physio = load_subject_physio(fp)
        # z-score before concetenation
        physio = zscore(physio)
        physio_all.append(physio)
    # stack or append to list
    if group_method == 'stack':
        return np.concatenate(physio_all)
    elif group_method == 'list':
        return physio_all


def load_nki_event_file():
    # We are ASSUMING that the event timings are the same across all subjects (e.g. no counterbalancing)
    events = pd.read_csv('data/dataset_nki/events/sub-A00031219_ses-BAS2_task-BREATHHOLD_acq-1400_events.tsv', 
                         sep='\t')
    return events


def load_params():
    # Load analysis parameters - filepath is defined as global variable at top of script
    with open(params_fp) as params_file:
        params = json.load(params_file)
    return params


def load_subject_func(fp, mask_bin, regress_global):
    # Load subject functional
    nifti = nb.load(fp)
    nifti_data = nifti.get_fdata()
    nifti_data = convert_2d(mask_bin, nifti_data)
    # regress out global signal, if specified
    if regress_global:
        nifti_data = regress_global_signal(nifti_data)
    return nifti_data


def load_subject_list(dataset, subject_list_fp):
    # given dataset and filepath, get list of subjects and scan/runs
    subj_df = pd.read_csv(subject_list_fp)
    # load subject list for a dataset
    if dataset == 'chang':
        subj = subj_df.subject.tolist()
        scan = [f'000{s}' if s <10 else f'00{s}' for s in subj_df.scan] 
    elif dataset == 'chang_bh':
        subj = subj_df.subject.tolist()
        scan = [f'000{s}' if s <10 else f'00{s}' for s in subj_df.scan]
    elif dataset == 'chang_cue':
        subj = [f'000{s}' if s <10 else f'00{s}' for s in subj_df.subject]
        scan = [f'000{s}' if s <10 else f'00{s}' for s in subj_df.scan]
    elif dataset == 'hcp':
        subj = subj_df.subject.tolist()
        scan = subj_df.lr.tolist()
    elif dataset == 'natview':
        subj = [f'0{s}' if s <10 else f'{s}' for s in subj_df.subject]
        scan = subj_df.scan.tolist()
    elif dataset == 'nki': 
        subj = subj_df.subject.tolist()
        scan = subj_df.session.tolist()
    elif dataset == 'nki_rest': 
        subj = subj_df.subject.tolist()
        scan = subj_df.session.tolist()
    elif dataset == 'spreng':
        subj = subj_df.subject.tolist()
        scan = subj_df.scan.tolist()
    elif dataset == 'toronto':
        subj = subj_df.subject.tolist()
        scan = subj_df.scan.tolist()
    elif dataset == 'yale':
        subj = subj_df.subject.tolist()
        scan = subj_df.scan.tolist()
    return subj, scan


def load_subject_physio(fp):    
    # Load physio
    physio_signal = np.loadtxt(fp)
    return physio_signal[:, np.newaxis]


def regress_global_signal(func_data, mask_nan=True):
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


def write_nifti(data, output_file, params):
    data_imp = impute_zero_voxels(data, params['zero_mask'], params['n_vert_orig'])
    mask_nii = nb.load(params['mask'])
    mask_bin = mask_nii.get_fdata() > 0
    nifti_4d = np.zeros(mask_nii.shape + (data_imp.shape[0],), 
                        dtype=data_imp.dtype)
    nifti_4d[mask_bin, :] = data_imp.T

    nifti_out = nb.Nifti2Image(nifti_4d, mask_nii.affine)
    nb.save(nifti_out, output_file)

