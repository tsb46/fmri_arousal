import numpy as np

from utils.butterworth_filters import butterworth_filter
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
	if params['data']['physio']['filter_params']['despike']:
		physio_sig = wavelet_despike(physio_sig, 
		                             params['data']['physio']['filter_params']['despike_params']['widths'],
		                             params['data']['physio']['filter_params']['despike_params']['min_snr'],
		                             params['data']['physio']['filter_params']['despike_params']['noise_perc'], 
		                             params['data']['physio']['filter_params']['despike_params']['window_size'], 
		                             params['data']['physio']['filter_params']['despike_params']['interpolation_window'])
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


def wavelet_despike(signal, widths, min_snr, noise_perc, window_size, 
                    interp_window, end_pad=20):
    l_w, r_w = interp_window
    signal_z = signal.copy()
    peaks = find_peaks_cwt(np.abs(signal_z), widths=widths, min_snr=min_snr, noise_perc=noise_perc, window_size=window_size)
    if len(peaks) > 0:
        for peak in peaks:
            if peak in range(end_pad,len(signal)-end_pad):
                peak_idx = (peak-(l_w+1), peak+(r_w))
                signal_z[peak_idx[0]:peak_idx[1]] = np.NaN
        not_nan = np.logical_not(np.isnan(signal_z))
        indices = np.arange(len(signal_z))
        interp = interp1d(indices[not_nan], signal_z[not_nan], fill_value="extrapolate")
        signal_despiked = interp(indices)
        return signal_despiked
    return signal_z




