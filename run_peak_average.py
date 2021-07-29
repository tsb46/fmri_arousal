import argparse
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


def find_physio_peaks(physio_ts, height, sample_dist=10):
	norm_ts = zscore(physio_ts[:,0])
	pos_peaks = find_peaks(norm_ts, height=height, distance=sample_dist)
	neg_peaks = find_peaks(-norm_ts, height=height, distance=sample_dist)
	return pos_peaks[0], neg_peaks[0]


def select_peaks(peaks, l_window, r_window, max_sample, n_samples):
	filtered_peaks = np.array([peak for peak in peaks 
	                          if peak >= l_window
	                          if (peak+r_window) <= max_sample])
	rand_peak_select = np.random.permutation(len(filtered_peaks))[:n_samples]
	return filtered_peaks[rand_peak_select]


def select_physio_ts(physio_select, eeg, rv, hr):
	physio_dict = {
		'eeg': eeg,
		'rv': rv,
		'hr': hr
	}
	return physio_dict[physio_select]


def write_results(peak_avg_pos, peak_avg_neg, physio_select, mask, level, subj_n, zero_mask, n_vert):
	if level == 'group':
		analysis_str = 'peak_avg_group'
	elif level == 'subject':
		analysis_str = f'peak_avg_s{subj_n}'
	analysis_str += f'_{physio_select}'

	pickle.dump(peak_avg_pos, open(f'{analysis_str}_pos_results.pkl', 'wb'))
	pickle.dump(peak_avg_neg, open(f'{analysis_str}_neg_results.pkl', 'wb'))
	write_nifti(peak_avg_pos, mask, f'{analysis_str}_pos', zero_mask, n_vert)
	write_nifti(peak_avg_neg, mask, f'{analysis_str}_neg', zero_mask, n_vert)



def run_main(physio_select, level, subj_n, l_window, r_window, 
             peak_thres, physio_params, mask, n_samples=100):
	# Load functional data
	if level == 'subject' and subj_n is None:
		raise Exception('Subject number must be provided for subject-level analysis')

	func_data, eeg, rv, hr, zero_mask, n_vert = load_data(level, mask, physio_params, subj_n)
	# Normalize func data
	func_data = zscore(func_data)
	# Select physio time series
	physio_ts = select_physio_ts(physio_select, eeg, rv, hr)
	# Find positive and negative peaks of physio ts
	ts_peaks_pos, ts_peaks_neg = find_physio_peaks(physio_ts, peak_thres)
	selected_peaks_pos = select_peaks(ts_peaks_pos, l_window, r_window, 
	                                  len(physio_ts), n_samples)
	selected_peaks_neg = select_peaks(ts_peaks_neg, l_window, r_window, 
	                                  len(physio_ts), n_samples)
	# Average functional data around peaks
	peak_avg_pos = average_peak_window(selected_peaks_pos, func_data, l_window, r_window)
	peak_avg_neg = average_peak_window(selected_peaks_neg, func_data, l_window, r_window)
	# Write out results
	write_results(peak_avg_pos, peak_avg_neg, physio_select, mask, level, 
	              subj_n, zero_mask, n_vert)


if __name__ == '__main__':
	"""Run main analysis"""
	parser = argparse.ArgumentParser(description='Run positive and negative peak averaging on '
	                                 'physio time series')
	parser.add_argument('-p', '--physio_ts', 
	                    help="""
	                    physiological (physio) time series label used as "trigger" for peak average - i.e. 
	                    functional data is averaged around fixed windows of its peaks. 
	                    """,
	                    choices=['eeg', 'rv', 'hr'],
	                    required=True)
	parser.add_argument('-l', '--level',
						help='subject or group level analysis',
						default='group',
						choices=['subject', 'group'],
						type=str)
	parser.add_argument('-s', '--subject_n',
						help='subject number for subject level analysis',
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
	parser.add_argument('-pp', '--physio_params', 
						help='file path to preprocessing params for physio signals',
						default='physio_proc_params.json',
						type=str)
	parser.add_argument('-m', '--mask', 
						help='file path to brain mask',
						default='masks/MNI152_T1_3mm_brain_mask.nii.gz',
						type=str)
	args_dict = vars(parser.parse_args())
	run_main(args_dict['physio_ts'], args_dict['level'], args_dict['subject_n'],
	         args_dict['left_window_size'], args_dict['right_window_size'], 
	         args_dict['peak_thres'], args_dict['physio_params'], args_dict['mask'])


