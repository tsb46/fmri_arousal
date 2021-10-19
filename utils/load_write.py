import json
import nibabel as nb 
import numpy as np
import os

from utils.butterworth_filters import filter_functional_data
from utils.load_utils import find_fps
from utils.physio_preprocess import preprocess_physio
from scipy.io import loadmat 
from scipy.stats import zscore

# path to the analysis parameter .json file. Never move this out of the base directory!
params_fp='analysis_params.json'
# path to 3mm brain mask file. Never move this out of the masks directory!
mask_fp = 'masks/MNI152_T1_3mm_brain_mask.nii.gz'


def convert_2d(mask, nifti_data):
	nonzero_indx = np.nonzero(mask)
	nifti_2d = nifti_data[nonzero_indx]
	return nifti_2d.T


def filter_zero_voxels(nifti_data, group_method):
	# Construct mask of vertices with std greater than 0
	if group_method == 'stack':
		orig_n_vert = nifti_data.shape[1]
		zero_mask = np.std(nifti_data, axis=0) > 0
		zero_mask_indx = np.where(zero_mask)[0]
		return nifti_data[:, zero_mask], zero_mask_indx, orig_n_vert
	elif group_method == 'list':
		orig_n_vert = nifti_data[0].shape[1]
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


def load_data(data, level, physio, load_physio, subj_n=None, scan_n=None, group_method='stack'):
	# Load analysis parameters
	with open(params_fp) as params_file:
		params = json.load(params_file)

	# Pull physio labels (if not already selected)
	params_data = params[data]
	if physio is None:
		physio = params_data['physio']

	# Load mask
	mask = nb.load(mask_fp).get_fdata() > 0

	# Pull file paths
	fps = find_fps(data, level, physio, params_data, subj_n, scan_n)
	# Pull data for subject level analysis
	if level == 'subject':
		print(f'Loading subject: {subj_n}')
		func_data = load_subject_func(fps['func'][0], mask, params_data)
		if load_physio:
			physio_proc = []
			for p in physio:
				p_proc = load_subject_physio(fps[p][0], params_data, data, p)
				physio_proc.append(p_proc)
		else:
			physio_proc = None
	# Pull data for group level analysis
	elif level == 'group':
		print('Loading all subjects')
		func_data = load_group_func(fps['func'], mask, params_data, group_method)
		if load_physio:
			physio_proc = []
			for p in physio:
				p_proc = load_group_physio(fps[p], params_data, data, p, group_method)
				physio_proc.append(p_proc)
		else:
			physio_proc = None

	func_data, zero_mask, n_vert_orig = filter_zero_voxels(func_data, group_method)
	return func_data, physio_proc, physio, zero_mask, n_vert_orig
	

def load_group_func(fps, mask, params, group_method):
	group_data = []
	for fp in fps:
		subj_data = load_subject_func(fp, mask, params)
		# Normalize data before concatenation
		subj_data = zscore(subj_data)
		group_data.append(subj_data)

	if group_method == 'stack':
		return np.concatenate(group_data)
	elif group_method == 'list':
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


def load_subject_func(fp, mask, params):
	# Load scan
	nifti = nb.load(fp)
	nifti_data = nifti.get_fdata()
	nifti.uncache()
	nifti_data = convert_2d(mask, nifti_data)
	nifti_data = filter_functional_data(nifti_data, params)
	return nifti_data


def load_subject_physio(fp, params, data, physio_label):	
	if data == 'chang':
		physio_signal = loadmat(fp)
		if physio_label == 'eeg':
			physio_signal = physio_signal['eeg_at'][:,0]
		elif physio_label == 'rv':
			physio_signal = physio_signal['rv'][:,0]
		elif physio_label == 'hr':
			physio_signal = physio_signal['hr'][:,0]
	else:
		physio_signal = np.loadtxt(fp)
	physio_signal_proc = preprocess_physio(physio_signal, params)
	return physio_signal_proc


def write_nifti(data, output_file, zero_mask, orig_n_vert):
	data_imp = impute_zero_voxels(data, zero_mask, orig_n_vert)
	mask = nb.load(mask_fp)
	mask_bin = mask.get_fdata() > 0
	nifti_4d = np.zeros(mask.shape + (data_imp.shape[0],), 
						dtype=data_imp.dtype)
	nifti_4d[mask_bin, :] = data_imp.T

	nifti_out = nb.Nifti2Image(nifti_4d, mask.affine)
	nb.save(nifti_out, output_file)




	

