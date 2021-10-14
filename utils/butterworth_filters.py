import json

from scipy.signal import butter, sosfiltfilt, sosfreqz


def butter_bandpass(lowcut, highcut, fs, order=5):
	nyq = 0.5 * fs
	low = lowcut / nyq
	high = highcut / nyq
	sos = butter(order, [low, high], analog=False, btype='band', output='sos')
	return sos


def butter_highpass(cutoff, fs, order=5):
	nyq = 0.5 * fs
	low = cutoff / nyq
	sos = butter(order, low, analog=False, btype='high', output='sos')
	return sos


def butter_lowpass(cutoff, fs, order=5):
	nyq = 0.5 * fs
	high = cutoff / nyq
	sos = butter(order, high, analog=False, btype='low', output='sos')
	return sos


def butterworth_filter(signals, lowcut, highcut, fs, filter_type):
	if filter_type == 'bandpass':
		sos = butter_bandpass(lowcut, highcut, fs)
	elif filter_type == 'lowpass':
		sos = butter_lowpass(highcut, fs)
	elif filter_type == 'highpass':
		sos = butter_highpass(lowcut, fs)
	elif filter_type == 'raw':
		return signals
	signals = sosfiltfilt(sos, signals, axis=0)
	return signals


def filter_functional_data(func_data, lowcut, highcut, params):
	# Get TR from analysis params
	tr = physio_params['tr']
	fs = 1/tr 
	if params['filter_type'] == 'bandpass':
		lowcut = params['bandpass']['low']
		highcut = params['bandpass']['high']
	elif params['filter_type'] == 'lowpass':
		lowcut = None
		highcut = params['lowpass']['high']
	elif params['filter_type'] == 'highpass':
		lowcut = params['highpass']['low']
		highcut = None
	elif params['filter_type'] == 'raw':
		lowcut = None
		highcut = None
	func_data = butterworth_filter(func_data, lowcut, highcut, 
	                               fs, params['filter_type'])
	return func_data