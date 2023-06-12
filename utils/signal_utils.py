import json
import neurokit2 as nk
import numpy as np
import pandas as pd

from scipy.interpolate import interp1d
from scipy.signal import find_peaks, hilbert, resample, butter, sosfiltfilt
from scipy.stats import iqr, zscore


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


def butterworth_filter(signals, lowcut, highcut, fs, filter_type, npad=1000):
    # handle pandas objects
    if (isinstance(signals, pd.Series)):
        input_type = 'series'
        signal_name = signals.name
        signals = signals.values # convert to array
    elif (isinstance(signals, pd.DataFrame)):
        input_type='dataframe'
        signal_cols = signals.columns
        signals = signals.values
    else:
        input_type='array'
    # Set filter params
    if filter_type == 'bandpass':
        sos = butter_bandpass(lowcut, highcut, fs)
    elif filter_type == 'lowpass':
        sos = butter_lowpass(highcut, fs)
    elif filter_type == 'highpass':
        sos = butter_highpass(lowcut, fs)
    elif filter_type == 'raw':
        return signals
    if np.ndim(signals) == 1:
        signals = signals[:, np.newaxis]
    # Median padding to reduce edge effects
    signals = np.pad(signals,[(npad, npad), (0, 0)], 'median')
    # backward and forward filtering
    signals = sosfiltfilt(sos, signals, axis=0)
    # Cut padding to original signal
    signals = signals[npad:-npad, :]
    # return to pandas dataframe or series if input
    if input_type == 'series':
        signals = pd.Series(np.squeeze(signals), name=signal_name)
    elif input_type == 'dataframe':
        signals = pd.DataFrame(signals, columns=signal_cols)
    return signals


def clip_spikes(ts, spike_thres=5):
    # Clip time series within bounds
    # Spike thres is set in z-score units, must convert to original units
    if ~isinstance(ts, pd.Series):
        ts = pd.Series(ts)
    ts_mean = ts.mean()
    ts_std = ts.std()
    ts_spike_thres = (spike_thres*ts_std) + ts_mean
    ts_clip = ts.clip(-ts_spike_thres, ts_spike_thres)
    return ts_clip

