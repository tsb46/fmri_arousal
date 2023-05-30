import numpy as np

from utils.signal_utils import butterworth_filter
from scipy.interpolate import interp1d
from scipy.signal import find_peaks_cwt
from scipy.stats import zscore


def filter_physio(signal, filter_type, lowcut, highcut, fs):
	if filter_type =='raw':
		return signal
	else:
		signal_filt = butterworth_filter(signal, lowcut, highcut, fs, filter_type)
		return signal_filt


def preprocess_physio(physio_sig, params, physio_label):
	# Define sampling rate
	fs = 1/params['tr']
	# Select filter params
	lowcut, highcut = select_filter_params(params['data']['physio']['filter_params'], physio_label)
	physio_sig_proc = filter_physio(physio_sig, params['data']['physio']['filter_params']['filter_choice'][physio_label],
			                       lowcut, highcut, fs)
	return physio_sig_proc


def select_filter_params(physio_params, physio_label):
	if physio_params['filter_choice'][physio_label] == 'raw':
		lowcut = None
		highcut = None
	else:
		if physio_params['filter_choice'][physio_label] == 'bandpass':
			lowcut = physio_params['bandpass']['low']
			highcut = physio_params['bandpass']['high']
		elif physio_params['filter_choice'][physio_label] == 'lowpass':
			highcut = physio_params['lowpass']['high']
			lowcut = None
		elif physio_params['filter_choice'][physio_label] == 'highpass':
			highcut = None
			lowcut = physio_params['highpass']['low']
	return lowcut, highcut




