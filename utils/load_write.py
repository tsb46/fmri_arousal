import json
import nibabel as nb 
import numpy as np
import os

from glob import glob
from utils.physio_preprocess import preprocess_physio
from scipy.io import loadmat 
from scipy.stats import zscore

# Global variable - subject number list
subj_n_list = [10, 12, 14, 16, 21, 27, 28, 29, 30, 31, 32]


def convert_2d(mask, nifti_data):
	nonzero_indx = np.nonzero(mask)
	nifti_2d = nifti_data[nonzero_indx]
	return nifti_2d.T


def filter_zero_voxels(nifti_data):
	orig_n_vert = nifti_data.shape[1]
	# Construct mask of vertices with std greater than 0
	zero_mask = np.std(nifti_data, axis=0) > 0
	zero_mask_indx = np.where(zero_mask)[0]
	return nifti_data[:, zero_mask], zero_mask_indx, orig_n_vert


def find_subj_fp(subj_n, folder):
	if subj_n < 10:
		nifti_file = glob(f'{folder}/sub_000{subj_n}*')
	else:
		nifti_file = glob(f'{folder}/sub_00{subj_n}*')
	# Ensure there is only one file associated with the supplied subject number 
	if len(nifti_file) > 1:
		raise Exception(f'More than one file in {folder} associated with input subject number {subj_n}')
	elif len(nifti_file) < 1:
		raise Exception(f'No file in {folder} associated with input subject number {subj_n}')
	return nifti_file[0]


def impute_zero_voxels(nifti_data, zero_mask, orig_n_vert):
	nifti_full = np.zeros([nifti_data.shape[0], orig_n_vert])
	nifti_full[:, zero_mask] = nifti_data
	return nifti_full


def load_data(level, mask_file, physio_params_fp, subj_n=None, func=True):
	if level == 'subject':
		print(f'Loading subject: {subj_n}')
		if func:
			func_data = load_subject_bold(subj_n, mask_file)
		eeg, rv, hr = load_subject_physio(subj_n, physio_params_fp)
	elif level == 'group':
		print('Loading all subjects')
		if func:
			func_data = load_group_bold(mask_file)
		eeg, rv, hr = load_group_physio(physio_params_fp)
	if func:
		func_data, zero_mask, n_vert_orig = filter_zero_voxels(func_data)
		return func_data, eeg, rv, hr, zero_mask, n_vert_orig
	else:
		return eeg, rv, hr


def load_group_bold(mask_file):
	group_data = []
	for subj in subj_n_list:
		subj_data = load_subject_bold(subj, mask_file)
		# Normalize data before concatenation
		subj_data = zscore(subj_data)
		group_data.append(subj_data)
	group_data = np.concatenate(group_data)
	return group_data


def load_group_physio(physio_params_fp):
	eeg_all = []; rv_all = []; hr_all = []
	for subj in subj_n_list:
		eeg, rv, hr = load_subject_physio(subj, physio_params_fp)
		eeg_all.append(zscore(eeg)) 
		rv_all.append(zscore(rv)) 
		hr_all.append(zscore(hr))
	eeg_all = np.concatenate(eeg_all)
	rv_all = np.concatenate(rv_all)
	hr_all = np.concatenate(hr_all)
	return eeg_all, rv_all, hr_all


def load_subject_bold(subj_n, mask_file, bold_dir='data/bold/proc3_filter_norm', 
					  sub_pref='sub_'):
	# Pull subject file path
	nifti_file = find_subj_fp(subj_n, bold_dir)
	# Load scan
	nifti = nb.load(nifti_file)
	nifti_data = nifti.get_fdata()
	nifti.uncache()
	# Load mask
	mask = nb.load(mask_file).get_fdata() > 0
	subj_data = convert_2d(mask, nifti_data)
	return subj_data


def load_subject_physio(subj_n, physio_params_fp, eeg_dir='data/eeg_at_raw', 
						rv_dir='data/physio_rv', hr_dir='data/physio_hr'):
	# Load physio preprocessing parameter json
	with open(physio_params_fp) as physio_file:
		physio_params = json.load(physio_file)
	# load physio signals
	physio_signals = []
	physio_keys = ['eeg_at', 'rv', 'hr']
	for signal_dir, p_key in zip([eeg_dir, rv_dir, hr_dir], physio_keys):
		physio_fp = find_subj_fp(subj_n, signal_dir)
		physio_signals.append(loadmat(physio_fp)[p_key])
	physio_signals_proc = preprocess_physio(*physio_signals, physio_params)
	return physio_signals_proc


def write_nifti(data, mask_file, output_file, zero_mask, orig_n_vert):
	data_imp = impute_zero_voxels(data, zero_mask, orig_n_vert)
	mask = nb.load(mask_file)
	mask_bin = mask.get_fdata() > 0
	nifti_4d = np.zeros(mask.shape + (data_imp.shape[0],), 
						dtype=data_imp.dtype)
	nifti_4d[mask_bin, :] = data_imp.T
	nifti_out = nb.Nifti2Image(nifti_4d, mask.affine)
	nb.save(nifti_out, output_file)




	

