import argparse
import json
import nibabel as nb 
import numpy as np
import pickle

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


def write_results(dataset, peak_avg, pos_neg, physio_select, 
                  level, subj_n, scan, 
                  zero_mask, n_vert):
	if level == 'group':
		analysis_str = f'{dataset}_peak_avg_group'
	elif level == 'subject':
		analysis_str = f'{dataset}_peak_avg_s{subj_n}'
	if scan is not None:
		analysis_str += f'_{scan}'
	analysis_str += f'_{physio_select}'
	analysis_str += f'_{pos_neg}'

	write_nifti(peak_avg, analysis_str, zero_mask, n_vert)


def run_main(dataset, level, physio, subj_n, scan_n, pos_neg, l_window, r_window, 
             peak_thres, n_samples=100):
	# Load data
	func_data, physio_sig, physio_labels, zero_mask, n_vert = load_data(dataset, level, physio=physio, load_physio=True, 
	                                                                    subj_n=subj_n, scan_n=scan_n) 
	# Normalize func data
	func_data = zscore(func_data)
	for physio_ts, label in zip(physio_sig, physio_labels):
		# Find peaks of physio ts
		ts_peaks = find_physio_peaks(physio_ts, peak_thres, pos_neg)
		selected_peaks = select_peaks(ts_peaks, l_window, r_window, 
		                                  len(physio_ts), n_samples)
		# Average functional data around peaks
		peak_avg = average_peak_window(selected_peaks, func_data, l_window, r_window)
		# Write out results
		write_results(dataset, peak_avg, pos_neg, label, level, 
		              subj_n, scan_n, zero_mask, n_vert)


if __name__ == '__main__':
	"""Run main analysis"""
	parser = argparse.ArgumentParser(description='Run positive and negative peak averaging on '
	                                 'physio time series')
	parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to run analysis on',
                        choices=['chang', 'choe', 'gu', 'nki', 'yale', 'hcp'], 
                        required=True,
                        type=str)
	parser.add_argument('-v', '--pos_v_neg',
						help='whether to average around positive or negative peaks',
						default='pos',
						choices=['pos', 'neg'],
						type=str)
	parser.add_argument('-l', '--level',
						help='subject or group level analysis',
						required=False,
						default='group',
						choices=['subject', 'group'],
						type=str)
	parser.add_argument('-p', '--physio',
						help='select physio - can provide multiple (separated by space)',
						required=False,
						default=None,
						action='append',
						type=str)
	parser.add_argument('-s', '--subject_n',
						help='subject number for subject level analysis',
						required=False,
						default=None,
						type=str)
	parser.add_argument('-scan', '--scan_n',
						help='scan number for subject level analysis (if multiple runs from same subject',
						required=False,
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
	parser.add_argument('-t', '--peak_thres',
	                    help='height (absolute value) threshold for peak detection - set in zscore '
	                    'normalized units, i.e. std. deviations from the mean', 
	                    required=False,
	                    default=1,
	                    type=float)
	args_dict = vars(parser.parse_args())
	run_main(args_dict['dataset'], args_dict['level'], args_dict['physio'], args_dict['subject_n'], 
	         args_dict['scan_n'], args_dict['pos_v_neg'], args_dict['left_window_size'], 
	         args_dict['right_window_size'], args_dict['peak_thres'])


