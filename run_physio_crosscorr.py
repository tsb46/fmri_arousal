import argparse
import nibabel as nb 
import numpy as np
import pickle

from utils.load_write import load_data
from scipy.signal import find_peaks
from scipy.stats import zscore


def write_results(level, subj_n, eeg, rv, hr):
	physio_dict = {
		'eeg': eeg,
		'rv': rv,
		'hr': hr
	}
	if level == 'group':
		analysis_str = 'physio_group'
	elif level == 'subject':
		analysis_str = f'physio_s{subj_n}'
	pickle.dump(physio_dict, open(f'{analysis_str}_results.pkl', 'wb'))



def run_main(level, subj_n, physio_params):
	# Load data
	if level == 'subject' and subj_n is None:
		raise Exception('Subject number must be provided for subject-level analysis')

	eeg, rv, hr = load_data(level, None, physio_params, subj_n, func=False)
	
	write_results(level, subj_n, eeg, rv, hr)


if __name__ == '__main__':
	"""Run main analysis"""
	parser = argparse.ArgumentParser(description='Run positive and negative peak averaging on '
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
	parser.add_argument('-p', '--physio_params', 
						help='file path to preprocessing params for physio signals',
						default='physio_proc_params.json',
						type=str)

	args_dict = vars(parser.parse_args())
	run_main(args_dict['level'], args_dict['subject_n'], 
	         args_dict['physio_params'])


