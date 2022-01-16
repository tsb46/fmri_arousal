import json
import nibabel as nb 
import numpy as np
import pandas as pd
import os

from utils.signal.butterworth_filters import filter_functional_data
from utils.load_utils import find_fps, print_filter_info
from utils.load_physio import preprocess_physio
from scipy.io import loadmat 
from scipy.stats import zscore

# path to the analysis parameter .json file. Never move this out of the base directory!
params_fp='analysis_params.json'
# path to 3mm dilated brain mask file. Never move this out of the masks directory!
mask_fp = 'masks/MNI152_T1_3mm_brain_mask_dilated2.nii.gz'


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


def load_data(data, level, physio, load_physio, subj_n=None, scan_n=None, 
              events=False, group_method='stack', verbose=True, filter_nan_voxels=True):
	params = load_params()

	# Pull physio labels (if not already selected)
	params_data = params[data]
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
		func_data = load_subject_func(fps['func'][0], mask, params_data)
		# If load subject and group method = list, place in list for downstream compatibility
		if group_method == 'list':
			func_data = [func_data]
		# If load physio option is specified, load and preprocess physio
		if load_physio:
			physio_proc = []
			for p in physio:
				p_proc = load_subject_physio(fps[p][0], params_data, data, p)
				if group_method == 'list':
					p_proc = [p_proc]
				physio_proc.append(p_proc)
		else:
			physio_proc = None
	# Pull data for group level analysis
	elif level == 'group':
		if verbose:
			print('Loading all subjects')
		func_data = load_group_func(fps['func'], mask, params_data, group_method)
		if load_physio:
			physio_proc = []
			for p in physio:
				p_proc = load_group_physio(fps[p], params_data, data, p, group_method)
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


def load_nki_event_file():
    # We are ASSUMING that the event timings are the same across all subjects (e.g. no counterbalancing)
    events = pd.read_csv('data/dataset_nki/events/A00057406_task_breathhold_events.tsv', sep='\t')
    return events


def load_params():
	# Load analysis parameters - filepath is defined as global variable at top of script
	with open(params_fp) as params_file:
		params = json.load(params_file)
	return params


def load_subject_func(fp, mask, params):
	# Load scan
	nifti = nb.load(fp)
	nifti_data = nifti.get_fdata()
	nifti.uncache()
	nifti_data = convert_2d(mask, nifti_data)
	nifti_data = filter_functional_data(nifti_data, params)
	return nifti_data


def load_subject_physio(fp, params, data, physio_label):	
	# Load AND preprocess physio
	if (data == 'chang') and (physio_label != 'csf'):
		physio_signal = loadmat(fp)
		if physio_label == 'eeg':
			physio_signal = physio_signal['eeg_at'][:,0]
		elif physio_label == 'rv':
			physio_signal = physio_signal['rv'][:,0]
		elif physio_label == 'hr':
			physio_signal = physio_signal['hr'][:,0]
	else:
		physio_signal = np.loadtxt(fp)
	physio_signal_proc = preprocess_physio(physio_signal, params, physio_label)
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




	

