from PyEMD import CEEMDAN
from scipy.signal import butter, sosfiltfilt, sosfreqz
from scipy.stats import zscore


def butter_bandpass(lowcut, highcut, fs, order=5):
	nyq = 0.5 * fs
	low = lowcut / nyq
	high = highcut / nyq
	sos = butter(order, [low, high], analog=False, btype='band', output='sos')
	return sos


def butter_highpass(cutoff, fs, order=5):
	nyq = 0.5 * fs
	high = cutoff / nyq
	sos = signal.butter(order, high, analog=False, btype='high', output='sos')
	return sos


def butter_lowpass(cutoff, fs, order=5):
	nyq = 0.5 * fs
	low = cutoff / nyq
	sos = butter(order, low, analog=False, btype='low', output='sos')
	return sos


def ensemble_emd(signal, max_imf, trials, parallel, processes):
	if ~parallel:
		processes = None
	ceemd = CEEMDAN(trials, parallel=parallel, processes=processes)
	imfs = ceemd(signal[:,0], max_imf=max_imf)
	return imfs


def filter_physio(filter_type, signal, lowcut, highcut, fs, 
                  order=5, ceemd_params=None):
	non_emd = ['bandpass', 'highpass', 'lowpass']
	if filter_type in non_emd:
		if filter_type == 'bandpass':
			sos = butter_bandpass(lowcut, highcut, fs, order=order)
		elif filter_type == 'lowpass':
			sos = butter_lowpass(lowcut, fs, order=order)
		elif filter_type == 'highpass':
			sos = butter_highpass(highcut, fs, order=order)
		# Use filtfilt to avoid phase delay
		data_filt = sosfiltfilt(sos, signal, axis=0)
	elif filter_type == 'ceemd':
		data_filt = ensemble_emd(signal, ceemd_params['max_imf'], 
		                         ceemd_params['trials'], ceemd_params['parallel'],
		                         ceemd_params['processes'])
	return data_filt


def preprocess_physio(eeg, rv, hr, physio_params):
	# Define sampling rate
	fs = 1/physio_params['tr']
	# Select filter params
	if physio_params['filter_choice'] == 'ceemd':
		ceemd_params = select_params(physio_params)
		highcut = None 
		lowcut = None
	elif physio_params['filter_choice'] in ['bandpass', 'lowpass', 'highpass']:
		lowcut, highcut = select_params(physio_params)
		ceemd_params = None
	elif physio_params['filter_choice'] == 'raw':
		return [eeg, rv, hr]
		
	# Filter signals
	proc_signal = []
	for signal in [eeg, rv, hr]:
		tmp_signal = filter_physio(physio_params['filter_choice'], signal, 
								   lowcut, highcut, fs, ceemd_params=ceemd_params)
		proc_signal.append(zscore(tmp_signal))
	return proc_signal


def select_params(physio_params):
	if physio_params['filter_choice'] != 'ceemd':
		if physio_params['filter_choice'] == 'bandpass':
			lowcut = physio_params['filter_params']['bandpass']['low']
			highcut = physio_params['filter_params']['bandpass']['high']
		elif physio_params['filter_choice'] == 'lowpass':
			lowcut = physio_params['filter_params']['lowpass']['low']
			highcut = None
		elif physio_params['filter_choice'] == 'highpass':
			lowcut = None
			highcut = physio_params['filter_params']['highpass']['high']
		return lowcut, highcut
	elif physio_params['filter_choice'] == 'ceemd':
		ceemd_params = {
			'max_imf': physio_params['filter_params']['ceemd']['max_imf'],
			'trials': physio_params['filter_params']['ceemd']['trials'],
			'parallel': physio_params['filter_params']['ceemd']['parallel'],
			"processes": physio_params['filter_params']['ceemd']['processes']
		}
		return ceemd_params


