import neurokit2 as nk
import numpy as np
import pandas as pd

from utils.rsp_rvt import rsp_rvt
from utils.signal_utils import butterworth_filter

# EEG bands of interest
band_freqs = [
    ('DELTA', 1, 3),
    ('THETA', 4, 7),
    ('ALPHA', 8, 12),
]


def compute_vigilance_rms(eeg, sf, window_s=2, alpha_indx=2, theta_indx=1):
    # calculate vigilance from ratio of alpha to theta using rolling window root mean square
    # 2 second window
    win = window_s*sf
    # Freq ranges 
    alpha_l, alpha_h = band_freqs[alpha_indx][1], band_freqs[alpha_indx][2]
    theta_l, theta_h = band_freqs[theta_indx][1], band_freqs[theta_indx][2]
    # Filter to alpha range
    eeg_freq = eeg.copy()
    eeg_freq.filter(alpha_l, alpha_h, l_trans_bandwidth=1, h_trans_bandwidth=1)
    alpha = eeg_freq.get_data()
    # Filter to theta range
    eeg_freq = eeg.copy()
    eeg_freq.filter(theta_l, theta_h, l_trans_bandwidth=1, h_trans_bandwidth=1)
    theta = eeg_freq.get_data()
    # Create temporary dataframe to use Pandas rolling function
    df = pd.DataFrame({
                       'alpha': np.mean(alpha, axis=0), 
                       'theta': np.mean(theta, axis=0)
                       })
    # Calculate vigilance by alpha/theta
    alpha_rolling_rms = df.rolling(win)['alpha'].apply(lambda x: np.sqrt(np.sum(x**2)), raw=True)
    theta_rolling_rms = df.rolling(win)['theta'].apply(lambda x: np.sqrt(np.sum(x**2)), raw=True)
    vigilance_at = alpha_rolling_rms/theta_rolling_rms
    return vigilance_at.fillna(0)


def eeg_band_amplitudes(eeg, eeg_bands, l_trans=1, h_trans=1):
    # Extract eeg band amplitudes from MNE object
    eeg_power = []
    for i, (band, fmin, fmax) in enumerate(eeg_bands):
        # bandpass filter
        eeg_freq = eeg.copy()
        eeg_freq.filter(fmin, fmax, n_jobs=1,
                        l_trans_bandwidth=l_trans,
                        h_trans_bandwidth=h_trans,
                        verbose=False)  
        # get analytic signal (envelope)
        eeg_freq.apply_hilbert(envelope=True)
        eeg_power.append(eeg_freq)
    return eeg_power


def extract_ecg_signals(ts, sf):
    # Extract R-peaks and heart rate from ECG
    # Clean EKG signal
    ecg_cleaned = nk.ecg_clean(ts, sampling_rate=sf)
    # Detect peaks (correct artifacts)
    rpeaks, info = nk.ecg_peaks(ecg_cleaned, sampling_rate=sf, correct_artifacts=True)
    ecg_rate = nk.signal_rate(rpeaks, sampling_rate=sf, desired_length=len(rpeaks))
    rpeaks['ECG_HR'] = ecg_rate
    return rpeaks[['ECG_HR']]


def extract_eeg_signals(eeg, sf):
    eeg_freq = eeg.copy()
    # Extract EEG Bands
    eeg_bands = eeg_band_amplitudes(eeg_freq, band_freqs, l_trans='auto', h_trans='auto')
    # Compute vigilance index
    vigilance_at = compute_vigilance_rms(eeg_freq, sf)
    # Run through bands and resample
    bands_ts = []
    band_labels = [b[0] for b in band_freqs]
    for i, (band, band_label) in enumerate(zip(eeg_bands, band_labels)):
        eeg_df = band.to_data_frame()
        eeg_df.drop(columns=['time'], inplace=True)
        eeg_avg = eeg_df.mean(axis=1)
        eeg_avg = eeg_avg.rename(band_label)
        bands_ts.append(eeg_avg)

    bands_df = pd.concat(bands_ts, axis=1)
    bands_df['VIGILANCE'] = vigilance_at
    return bands_df


def extract_gsr_signals(ts, sf, bp_filt=[0.01, 0.1]):
    # extract tonic skin conductance from electrodermal recordings
    ts_detrend = nk.signal_detrend(ts, method='polynomial', order=3)
    ts_filt = butterworth_filter(ts_detrend, bp_filt[0], bp_filt[1], sf, 'bandpass')
    return pd.DataFrame({'GSR': ts_filt[:,0]})


def extract_ppg_signals(ts, sf, w=6):
    # Process PPG signal
    # Extract R-peaks and heart rate from PPG
    # Clean PPG signal
    ppg_signals, ppg_info = nk.ppg_process(ts, sampling_rate=sf)
    # Band-pass filter to 0.01 - 0.1Hz
    ppg_signals['PPG_LOW'] = butterworth_filter(ts, 0.01, 0.1, sf, 'bandpass')
    ## PPG Peak Amplitude
    ppg_peaks_loc = np.where(ppg_signals['PPG_Peaks'])[0]
    ppg_peaks_amp = np.abs(ppg_signals['PPG_Clean'].iloc[ppg_peaks_loc])
    ppg_signals['PPG_PEAK_AMP'] = nk.signal_interpolate(ppg_peaks_loc, ppg_peaks_amp.values, 
                                                        np.arange(ppg_signals.shape[0]), 
                                                        method='cubic')
    ppg_signals.rename(columns={'PPG_Rate': 'PPG_HR'}, inplace=True)
    return ppg_signals[['PPG_HR', 'PPG_LOW', 'PPG_PEAK_AMP']]


def extract_resp_signals(ts, sf):
    # extract respiration amplitude and frequency
    r_signals, r_info = nk.rsp_process(ts, sampling_rate=sf)
    # compute rvt = breathing frequency * breathing amplitude
    rvt, rvt_amp, rvt_if = rsp_rvt(r_signals['RSP_Clean'], sampling_rate=sf)
    resp_df = pd.DataFrame({
        'RSP_RVT': rvt,
        'RSP_RVT_AMP': rvt_amp,
        'RSP_RVT_IF': rvt_if
    })
    return resp_df
    