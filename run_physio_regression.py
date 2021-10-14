import argparse
import json
import numpy as np
import pickle

from scipy.stats import gamma, zscore
from sklearn.linear_model import LinearRegression
from utils.load_write import load_data, write_nifti
from utils.butterworth_filters import filter_functional_data


def convolve_hrf(hrf, ts):
	n_drop = len(hrf) - 1
	convolved_events = np.convolve(ts, hrf)
	return convolved_events[:-n_drop]


def double_gamma_hrf(t, tr, dip=0.35):
	# http://www.jarrodmillman.com/rcsds/lectures/convolution_background.html
	n_steps = np.arange(0, t, tr)
	gamma_peak = gamma.pdf(n_steps, 6)
	gamma_under = gamma.pdf(n_steps, 12)
	gamma_double = gamma_peak - dip * gamma_under
	return gamma_double / np.max(gamma_double) * 0.6


def lag(arr, num, fill_value=0):
	# https://stackoverflow.com/questions/30399534/shift-elements-in-a-numpy-array
    result = np.empty_like(arr)
    if num > 0:
        result[:num] = fill_value
        result[num:] = arr[:-num]
    elif num < 0:
        result[num:] = fill_value
        result[:num] = arr[-num:]
    else:
        result[:] = arr
    return result


def linear_regression(eeg, rv, hr, func_data, tr, lag_n, model_type, convolve):
	# Simple OLS - no mixed/multilevel model
	signal_lags_all = []
	for signal in [eeg, rv, hr]:
		signal_lags = []
		for l in range(lag_n+1):
			if convolve:
				hrf = double_gamma_hrf(30, tr)
				signal = convolve_hrf(hrf, signal)
			signal_lags.append(lag(signal, l))
		if len(signal_lags) > 1:
			signal_lags_all.append(np.stack(signal_lags, axis=1))
		else:
			signal_lags_all.append(signal_lags[0][:, np.newaxis])
	beta_maps = []
	for l in range(lag_n+1):
		if model_type == 'group':
			x = np.concatenate([signal_lags_all[0][:,l], 
			                   signal_lags_all[1][:,l],
			                   signal_lags_all[2][:,l]], axis=1)
			lin_reg = LinearRegression()
			lin_reg.fit(zscore(x), zscore(func_data))
			beta_maps.append(lin_reg.coef_)
		else:
			for signal_lags in signal_lags_all:
				x = signal_lags[:,l]
				lin_reg = LinearRegression()
				lin_reg.fit(zscore(x).reshape(-1,1), zscore(func_data))
				beta_maps.append(lin_reg.coef_)
	return beta_maps


def organize_beta_maps(beta_maps, lag_n, model_type):
	eeg_maps = []
	rv_maps = []
	hr_maps = []
	indx = 0
	for l in range(lag_n + 1):
		if model_type == 'group':
			eeg_maps.append(beta_maps[l][0])
			rv_maps.append(beta_maps[l][1])
			hr_maps.append(beta_maps[l][2])
		else:
			eeg_maps.append(beta_maps[indx])
			rv_maps.append(beta_maps[indx+1])
			hr_maps.append(beta_maps[indx+2])
			indx += 3
	if lag_n > 0:
		eeg_maps = np.concatenate(eeg_maps, axis=1)
		rv_maps = np.concatenate(rv_maps, axis=1)
		hr_maps = np.concatenate(hr_maps, axis=1)
	else:
		eeg_maps = eeg_maps[0]
		rv_maps = rv_maps[0]
		hr_maps = hr_maps[0]
	return eeg_maps, rv_maps, hr_maps


def write_results(level, subj_n, scan, mask, zero_mask, n_vert,
                  eeg_maps, rv_maps, hr_maps):
	if level == 'group':
		analysis_str = 'physio_reg_group'
	elif level == 'subject':
		analysis_str = f'physio_reg_s{subj_n}'

	if scan is not None:
		analysis_str += f'_{scan}'

	write_nifti(eeg_maps.T, mask, f'{analysis_str}_eeg', zero_mask, n_vert)
	write_nifti(rv_maps.T, mask, f'{analysis_str}_rv', zero_mask, n_vert)
	write_nifti(hr_maps.T, mask, f'{analysis_str}_hr', zero_mask, n_vert)


def run_main(level, subj_n, scan_n, model_type, convolve, lowcut, 
             highcut, lag_n, physio_params_fp, mask):
	# Load data
	if level == 'subject' and subj_n is None:
		raise Exception('Subject number must be provided for subject-level analysis')
	func_data, eeg, rv, hr, zero_mask, n_vert = load_data(level, mask, 
	                                                      physio_params_fp, subj_n, 
	                                                      scan_n)
	func_data = filter_functional_data(func_data, lowcut, highcut, physio_params_fp)
	# Get TR from physio params
	with open(physio_params_fp) as physio_file:
		physio_params = json.load(physio_file)
	tr = physio_params['tr']
	# run linear regression
	beta_maps = linear_regression(eeg, rv, hr, func_data, tr, lag_n, 
	                              model_type, convolve)
	eeg_maps, rv_maps, hr_maps = organize_beta_maps(beta_maps, lag_n, model_type) 
	write_results(level, subj_n, scan_n, mask, zero_mask, n_vert,
	              eeg_maps, rv_maps, hr_maps)


if __name__ == '__main__':
	"""Run main analysis"""
	parser = argparse.ArgumentParser(description='Regress functional data on '
	                                 'physio time series')
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
	parser.add_argument('-t', '--model_type',
						help='group regression of all physio signals, or'
						' individual regressions of each signal',
						default='individual',
						choices=['individual', 'group'],
						type=str)
	parser.add_argument('-c', '--convolve',
						help='whether to convolve physio time series with canonical '
						'hemodynamic response function',
						default=0,
						choices=[0,1],
						type=int)
	parser.add_argument('-fl', '--low_cut', 
						help='lowcut for butterworth filter for functional data',
						default=None,
						type=float)	
	parser.add_argument('-fh', '--high_cut', 
						help='highcut for butterworth filter for functional data',
						default=None,
						type=float)	
	parser.add_argument('-nlag', '--lag_number', 
						help='number of lags of physio regressor to include in model',
						default=0,
						type=int)
	parser.add_argument('-p', '--physio_params', 
						help='file path to preprocessing params for physio signals',
						default='physio_proc_params.json',
						type=str)
	parser.add_argument('-m', '--mask', 
						help='file path to brain mask',
						default='masks/MNI152_T1_3mm_brain_mask.nii.gz',
						type=str)

	args_dict = vars(parser.parse_args())
	run_main(args_dict['level'], args_dict['subject_n'], args_dict['scan_n'],
	         args_dict['model_type'], args_dict['convolve'], 
	         args_dict['low_cut'], args_dict['high_cut'],
	         args_dict['lag_number'], args_dict['physio_params'],
	         args_dict['mask'])


