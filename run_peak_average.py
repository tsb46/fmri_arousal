import argparse
import json
import nibabel as nb 
import numpy as np
import pickle

from utils.load_write import load_data, write_nifti
from scipy.signal import find_peaks
from scipy.stats import zscore


def average_peak_window(peak_indx, group_data, l_window, r_window):
	windows = []
	for peak in peak_indx:
		l_edge = peak - l_window
		r_edge = peak + r_window
		windows.append(group_data[l_edge:r_edge, :])
	windows_array = np.dstack(windows)
	return np.mean(windows_array, axis=2)


def find_physio_peaks(physio_ts, min_thres, max_thres, peak_distance, pos_neg):
	norm_ts = zscore(physio_ts)
	if pos_neg == 'neg':
		norm_ts = norm_ts*-1
	peaks = find_peaks(norm_ts, height=min_thres, distance=peak_distance)

	if max_thres is not None:
		peaks_max_indx = np.where(peaks[1]['peak_heights']<max_thres)[0]
		peaks = peaks[0][peaks_max_indx]
	else:
		peaks = peaks[0]
	return peaks


def select_peaks(peaks, l_window, r_window, max_sample, n_samples):
	filtered_peaks = np.array([peak for peak in peaks 
	                          if peak >= l_window
	                          if (peak+r_window) <= max_sample])
	rand_peak_select = np.random.permutation(len(filtered_peaks))[:n_samples]
	return filtered_peaks[rand_peak_select]


def write_results(dataset, peak_avg, pos_neg, physio_select, 
                  zero_mask, n_vert, params):
	analysis_str = f'{dataset}_peak_avg_group_{physio_select}_{pos_neg}'
	write_nifti(peak_avg, analysis_str, zero_mask, n_vert, params['mask'])


def run_main(dataset, physio, pos_neg, l_window, r_window, 
             min_peak_thres, max_peak_thres, peak_distance, n_samples=100):
	# Load data
	func_data, physio_sig, physio_labels, zero_mask, n_vert, params = load_data(dataset, 
	                                                                            physio=[physio], 
	                                                                            load_physio=True) 
	# Normalize func data
	func_data = zscore(func_data)
	for physio_ts, label in zip(physio_sig, physio_labels):
		# squeeze out extra dimension
		physio_ts = np.squeeze(physio_ts)
		# Find peaks of physio ts
		ts_peaks = find_physio_peaks(physio_ts, min_peak_thres, max_peak_thres, 
		                             peak_distance, pos_neg)
		selected_peaks = select_peaks(ts_peaks, l_window, r_window, len(physio_ts), n_samples)
		# Average functional data around peaks
		peak_avg = average_peak_window(selected_peaks, func_data, l_window, r_window)
		# Write out results
		write_results(dataset, peak_avg, pos_neg, label, zero_mask, n_vert, params)


if __name__ == '__main__':
	"""Run main analysis"""
	parser = argparse.ArgumentParser(description='Run positive and negative peak averaging on '
	                                 'physio time series')
	parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to run analysis on',
                        choices=['chang', 'chang_bh', 'nki', 'yale', 'hcp', 'spreng'], 
                        required=True,
                        type=str)
	parser.add_argument('-v', '--pos_v_neg',
						help='whether to average around positive or negative peaks',
						default='pos',
						choices=['pos', 'neg'],
						type=str)
	parser.add_argument('-p', '--physio',
						help='select physio - can provide multiple (separated by space)',
						required=True,
						default=None,
						type=str)
	parser.add_argument('-lw', '--left_window_size',
	                    help='Length of left window from selected peak', 
	                    required=False,
	                    default=0,
	                    type=int)
	parser.add_argument('-rw', '--right_window_size',
	                    help='Length of right window from selected peak', 
	                    required=False,
	                    default=15,
	                    type=int)
	parser.add_argument('-t_min', '--min_peak_thres',
	                    help='min height (absolute value) threshold for peak detection - set in zscore '
	                    'normalized units, i.e. std. deviations from the mean', 
	                    required=False,
	                    default=2,
	                    type=float)
	parser.add_argument('-t_max', '--max_peak_thres',
	                    help='max height (absolute value) threshold for peak detection - set in zscore '
	                    'normalized units, i.e. std. deviations from the mean', 
	                    required=False,
	                    default=None,
	                    type=float)
	parser.add_argument('-dist', '--peak_distance',
	                    help='minimum distance between peaks', 
	                    required=False,
	                    default=10,
	                    type=float)
	args_dict = vars(parser.parse_args())
	run_main(args_dict['dataset'], args_dict['physio'], args_dict['pos_v_neg'], 
	         args_dict['left_window_size'], args_dict['right_window_size'], 
	         args_dict['min_peak_thres'], args_dict['max_peak_thres'], 
	         args_dict['peak_distance'])


