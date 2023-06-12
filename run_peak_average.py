import argparse
import json
import nibabel as nb 
import numpy as np
import pickle

from utils.load_write import load_data, write_nifti
from scipy.signal import find_peaks
from scipy.stats import zscore


def average_peak_window(peak_indx, group_data, l_window, r_window):
	# average BOLD signals to the left and/or right of peaks
	windows = []
	for peak in peak_indx:
		l_edge = peak - l_window
		r_edge = peak + r_window
		windows.append(group_data[l_edge:r_edge, :])
	windows_array = np.dstack(windows)
	return np.mean(windows_array, axis=2)


def find_physio_peaks(physio_ts, min_thres, peak_distance):
	# find peaks in physio signal based on peak parameters
	norm_ts = zscore(physio_ts)
	peaks = find_peaks(norm_ts, height=min_thres, distance=peak_distance)
	return peaks[0]


def select_peaks(peaks, l_window, r_window, max_sample, n_samples):
	# select peaks found from 'find_physio_peaks'
	filtered_peaks = np.array([peak for peak in peaks 
	                          if peak >= l_window
	                          if (peak+r_window) <= max_sample])
	rand_peak_select = np.random.permutation(len(filtered_peaks))[:n_samples]
	return filtered_peaks[rand_peak_select]


def write_results(dataset, peak_avg, physio_select, 
                  zero_mask, n_vert, out_dir):
	# write out results
	if out_dir is not None:
		analysis_str = f'{out_dir}/{dataset}_peak_avg_{physio_select}'
	else:
		analysis_str = f'{dataset}_peak_avg_{physio_select}'
	write_nifti(peak_avg, analysis_str, zero_mask, n_vert)


def run_peak_average(dataset, physio, l_window, r_window, min_peak_thres, 
                     peak_distance, out_dir=None, n_samples=100):
	# Load data
	func_data, physio_sig, zero_mask, n_vert = load_data(dataset, physio=physio) 
	# squeeze out extra dimension
	physio_sig = np.squeeze(physio_sig)
	# Normalize func data
	func_data = zscore(func_data)
	# Get physio peaks
	ts_peaks = find_physio_peaks(physio_sig, min_peak_thres, peak_distance)
	# select windows around peaks
	selected_peaks = select_peaks(ts_peaks, l_window, r_window, len(physio_sig), n_samples)
	# Average functional data around peaks
	peak_avg = average_peak_window(selected_peaks, func_data, l_window, r_window)
	# Write out results
	write_results(dataset, peak_avg, physio, zero_mask, n_vert, out_dir)


if __name__ == '__main__':
	"""Run main analysis"""
	parser = argparse.ArgumentParser(description='Run positive and negative peak averaging on '
	                                 'physio time series')
	parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to run analysis on',
                        choices=['chang', 'chang_bh', 'change_cue', 
                        		 'nki', 'yale', 'hcp', 'spreng', 
                        		 'natview'], 
                        required=True,
                        type=str)
	parser.add_argument('-p', '--physio',
						help='select physio',
						required=True,
						choices=['PPG_HR', 'ECG_HR', 'PPG_PEAK_AMP', 'PPG_LOW', 
                                 'RSP_RVT', 'GSR', 'ALPHA', 'THETA', 'DELTA', 
                                 'PUPIL'],
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
	                    default=10,
	                    type=int)
	parser.add_argument('-t_min', '--min_peak_thres',
	                    help='min height (absolute value) threshold for peak detection - set in zscore '
	                    'normalized units, i.e. std. deviations from the mean', 
	                    required=False,
	                    default=2,
	                    type=float)
	parser.add_argument('-dist', '--peak_distance',
	                    help='minimum distance between peaks', 
	                    required=False,
	                    default=10,
	                    type=float)
	args_dict = vars(parser.parse_args())
	run_peak_average(args_dict['dataset'], args_dict['physio'],  
	                 args_dict['left_window_size'], args_dict['right_window_size'], 
	                 args_dict['min_peak_thres'], args_dict['peak_distance'])


