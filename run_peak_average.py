import argparse
import json
import nibabel as nb 
import numpy as np
import pickle

from utils.butterworth_filters import filter_functional_data
from utils.load_write import load_data, write_nifti
from scipy.signal import find_peaks
from scipy.stats import zscore

def average_peak_window(peak_indx, group_data, l_window, 
                        r_window, return_peak_ts=False):
	windows = []
	for peak in peak_indx:
		l_edge = peak - l_window
		r_edge = peak + r_window
		windows.append(group_data[l_edge:r_edge, :])
	windows_array = np.dstack(windows)
	if return_peak_ts:
		return np.mean(windows_array, axis=2), np.vstack(windows)
	else:
		return np.mean(windows_array, axis=2)


def find_physio_peaks(physio_ts, height, pos_neg, sample_dist=10):
	norm_ts = zscore(physio_ts)
	if pos_neg == 'neg':
		norm_ts = norm_ts*-1
	peaks = find_peaks(norm_ts, height=height, distance=sample_dist)
	return peaks[0]


def select_peaks(peaks, l_window, r_window, max_sample, n_samples):
	filtered_peaks = np.array([peak for peak in peaks 
	                          if peak >= l_window
	                          if (peak+r_window) <= max_sample])
	rand_peak_select = np.random.permutation(len(filtered_peaks))[:n_samples]
	return filtered_peaks[rand_peak_select]


def write_results(peak_avg, pos_neg, physio_select, 
                  mask, level, subj_n, scan, 
                  zero_mask, n_vert):
	if level == 'group':
		analysis_str = 'peak_avg_group'
	elif level == 'subject':
		analysis_str = f'peak_avg_s{subj_n}'
	if scan is not None:
		analysis_str += f'_{scan}'
	analysis_str += f'_{physio_select}'
	analysis_str += f'_{pos_neg}'

	write_nifti(peak_avg, mask, analysis_str, zero_mask, n_vert)


def run_main(level, subj_n, scan_n, pos_neg, l_window, r_window, 
             low_cut, high_cut, peak_thres, physio_params_fp, mask, n_samples=50):
	# Load functional data
	if level == 'subject' and subj_n is None:
		raise Exception('Subject number must be provided for subject-level analysis')
	func_data, eeg, rv, hr, zero_mask, n_vert = load_data(level, mask, physio_params_fp, 
	                                                      subj_n, scan_n)
	# Filter functional data if specified
	func_data = filter_functional_data(func_data, low_cut, high_cut, physio_params_fp)
	# Normalize func data
	func_data = zscore(func_data)
	for physio_ts, label in zip([eeg, rv, hr], ['eeg', 'rv', 'hr']):
		# Find peaks of physio ts
		ts_peaks = find_physio_peaks(physio_ts, peak_thres, pos_neg)
		selected_peaks = select_peaks(ts_peaks, l_window, r_window, 
		                                  len(physio_ts), n_samples)
		# Average functional data around peaks
		peak_avg = average_peak_window(selected_peaks, func_data, l_window, r_window)
		# Write out results
		write_results(peak_avg, pos_neg, label, mask, level, 
		              subj_n, scan_n, zero_mask, n_vert)


if __name__ == '__main__':
	"""Run main analysis"""
	parser = argparse.ArgumentParser(description='Run positive and negative peak averaging on '
	                                 'physio time series')
	parser.add_argument('-p', '--pos_v_neg',
						help='whether to average around positive or negative peaks',
						default='pos',
						choices=['pos', 'neg'],
						type=str)
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
	parser.add_argument('-lw', '--left_window_size',
	                    help='Length of left window from selected peak', 
	                    required=False,
	                    default=10,
	                    type=int)
	parser.add_argument('-rw', '--right_window_size',
	                    help='Length of right window from selected peak', 
	                    required=False,
	                    default=10,
	                    type=int)
	parser.add_argument('-fl', '--low_cut', 
						help='lowcut for butterworth filter for functional data',
						default=None,
						type=float)	
	parser.add_argument('-fh', '--high_cut', 
						help='highcut for butterworth filter for functional data',
						default=None,
						type=float)	
	parser.add_argument('-t', '--peak_thres',
	                    help='height (absolute value) threshold for peak detection - set in zscore '
	                    'normalized units, i.e. std. deviations from the mean', 
	                    required=False,
	                    default=1,
	                    type=float)
	parser.add_argument('-pp', '--physio_params', 
						help='file path to preprocessing params for physio signals',
						default='physio_proc_params.json',
						type=str)
	parser.add_argument('-m', '--mask', 
						help='file path to brain mask',
						default='masks/MNI152_T1_3mm_brain_mask.nii.gz',
						type=str)
	args_dict = vars(parser.parse_args())
	run_main(args_dict['level'], args_dict['subject_n'], args_dict['scan_n'],
	         args_dict['pos_v_neg'], args_dict['left_window_size'], 
	         args_dict['right_window_size'], args_dict['low_cut'],
	         args_dict['high_cut'], args_dict['peak_thres'], 
	         args_dict['physio_params'], args_dict['mask'])


