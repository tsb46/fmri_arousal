import neurokit2 as nk
import numpy as np
import pandas as pd

from scipy.interpolate import interp1d
from scipy.signal import find_peaks, hilbert, resample
from utils.signal.butterworth_filters import butterworth_filter 
from scipy.stats import iqr


def eeg_band_amplitudes(eeg, eeg_bands, l_trans=1, h_trans=1):
    # Extract eeg band amplitudes from MNE object
    eeg_freqs_all = []
    for i, (band, fmin, fmax) in enumerate(eeg_bands):
        # bandpass filter
        eeg_freq = eeg.copy()
        eeg_freq.filter(fmin, fmax, n_jobs=1,  # use more jobs to speed up.
                        l_trans_bandwidth=l_trans,  # make sure filter params are the same
                        h_trans_bandwidth=h_trans)  # in each band and skip "auto" option.
        # get analytic signal (envelope)
        eeg_freq.apply_hilbert(envelope=True)
        eeg_freqs_all.append(eeg_freq)
    return eeg_freqs_all    


def filt_resample_physio_to_func(ts, high_cut, func_len, sf_physio):
    # Lowpass filter to avoid aliasing
    ts_filt = butterworth_filter(ts, None, high_cut, sf_physio, 'lowpass')
    ts_func = nk.signal_resample(ts_filt, desired_length=func_len,
                                 method='FFT')
    return ts_func


def find_interpolate_spikes(ts, spike_thres=6):
    peaks_info = nk.signal_findpeaks(ts, relative_height_min=spike_thres, relative_median=True)
    for onset, offset in zip(peaks_info['Onsets'], peaks_info['Offsets']):
        if (~np.isnan(onset)) and (~np.isnan(offset)):
            onset = int(onset)
            offset = int(offset)
            ts[onset:offset] = np.nan

    not_nan = np.logical_not(np.isnan(ts))
    indices = np.arange(len(ts))
    interp = interp1d(indices[not_nan], ts[not_nan], kind='cubic', 
                      fill_value="extrapolate")
    ts_cleaned = interp(indices)
    return ts_cleaned


def get_peaks_ppg(ppg_sig, ppg_rng, sf_physio):
    minHeight = 0.05*ppg_rng
    minDist = (sf_physio*.1)+(sf_physio/2)
    locs, pks_d = find_peaks(ppg_sig, height = minHeight, distance = minDist)
    return locs, pks_d 


def hilbert_resp_amplitude(ts, sf, bp_filt=[0.05, 3]):
    ts_filt = butterworth_filter(ts, bp_filt[0], bp_filt[1], sf, 'bandpass')    
    ts_amp = np.abs(hilbert(ts_filt))
    return ts_amp


def interbeat_interval_ppg(ts, sf_physio, bp_filt=[0.5,2]):
    # Get TR of physio
    tr_physio = 1 / sf_physio

    # Bandpass filter PPG signal
    ppg_bpf = butterworth_filter(ts, bp_filt[0], bp_filt[1], sf_physio, 'bandpass')

    # Get IQR OF PPG (25%, 75%)
    ppg_rng = iqr(ppg_bpf)

    # Get peaks of PPG
    locs, pk_d = get_peaks_ppg(ppg_bpf, ppg_rng, sf_physio)
    ppg_samples = locs
    ppg_trig_times = ppg_samples*tr_physio

    # Calculate interbeat interval
    IBI = (np.diff(ppg_trig_times))
    t_IBI = 0.5*(ppg_trig_times[1:] + ppg_trig_times[0:-1])

    return IBI, t_IBI


def nk_extract_ecg_signals(ts, sf):
    # Extract R-peaks and heart rate from ECG
    # Clean EKG signal
    ecg_cleaned = nk.ecg_clean(ts, sampling_rate=sf)
    # Detect peaks (correct artifacts)
    rpeaks, info = nk.ecg_peaks(ecg_cleaned, sampling_rate=sf, correct_artifacts=True)
    ecg_rate = nk.signal_rate(rpeaks, sampling_rate=sf, desired_length=len(rpeaks))
    rpeaks['ECG_Rate'] = ecg_rate
    return rpeaks


def nk_extract_ppg_signals(ts, sf):
    # Process PPG signal
    # Extract R-peaks and heart rate from PPG
    # Clean PPG signal
    ppg_signals, ppg_info = nk.ppg_process(ts, sampling_rate=sf)
    return ppg_signals


def nk_extract_resp_signals(ts, sf):
    # extract respiration amplitude and frequency
    r_signals, r_info = nk.rsp_process(ts, sampling_rate=sf)
    # compute rvt = breathing frequency * breathing amplitude
    r_signals['RSP_RVT'] = r_signals['RSP_Rate'] * r_signals['RSP_Amplitude']
    r_signals['RSP_AMP_HILBERT'] = hilbert_resp_amplitude(ts, sf)
    return r_signals


def nk_extract_map(ts_dbp, ts_sbp, sf, bp_filt=(0.01, 0.1)):
    # filter to typical resting-state BOLD range
    ts_dbp_filt = butterworth_filter(ts_dbp, bp_filt[0], bp_filt[1], sf, 'bandpass')
    ts_sbp_filt = butterworth_filter(ts_sbp, bp_filt[0], bp_filt[1], sf, 'bandpass')
    # Calculate MAP
    ts_map = ts_dbp_filt*(2/3) + ts_sbp_filt*(1/3)
    return pd.DataFrame(ts_map, columns=['MAP'])


def nk_extract_physio(ts, phys_label, sf_physio, func_len, lowpass):
    # Neurokit physio preprocessing - output of nk extractor functions are Pandas dfs
    # If blood pressure, separate dbp and sbp
    if phys_label == 'bpp':
        # Assume dbp is first indx (double check this )
        ts_dbp, ts_sbp = ts[0], ts[1]
        phys_ts = nk_extract_map(ts_dbp, ts_sbp, sf_physio)
    elif phys_label == 'ecg':
        phys_ts = nk_extract_ecg_signals(ts, sf_physio)
    elif phys_label == 'resp':
        phys_ts = nk_extract_resp_signals(ts, sf_physio)
    elif phys_label == 'ppg':
        phys_ts = nk_extract_ppg_signals(ts, sf_physio)

    # resample physio to functional scans
    phys_ts_r = phys_ts.apply(filt_resample_physio_to_func, high_cut=lowpass, 
                              func_len=func_len, sf_physio=sf_physio, axis=0)

    return phys_ts_r


def trigger_extract_physio(ts, phys_label, func_scan_len, sf_func, 
                           sf_physio, win=3, despike=True):
    """
    Custom window-based physio preprocessing from Chang lab
    May. 2021
    Author: Karuna Gujar (modified by Taylor Bolt)
    PHYSIO LIBRARY: FMRI Regressors
    https://github.com/neurdylab/physio_scripts (private)
    """
    # Get sampling rate for func
    tr_func = 1 / sf_func
    tr_physio = 1 / sf_physio


    # If ppg, do peak detection
    if phys_label == 'ppg':
        ppg_ibi, ppg_t_ibi = interbeat_interval_ppg(ts, sf_physio)
    # If resp, filter to typical breathing frequency (0.05 - 3Hz)
    elif phys_label == 'resp':
        ts = pd.Series(butterworth_filter(ts, 0.05, 3, sf_physio, 'bandpass'))
    # If blood pressure, separated dbp and sbp and filter to typical resting-state fMRI range
    elif phys_label == 'bpp':
        # Assume dbp is first indx (double check)
        ts_dbp, ts_sbp = ts[0], ts[1]
        ts_dbp = pd.Series(butterworth_filter(ts_dbp, 0.01, 0.1, sf_physio, 'bandpass'))
        ts_sbp = pd.Series(butterworth_filter(ts_sbp, 0.01, 0.1, sf_physio, 'bandpass'))

    # Get time stamps of window centered on func TR
    t_fmri = trigger_window_center(tr_func, func_scan_len)

    t_phys = []
    # Loop through time points in func and calculate physio index within window
    for kk in range(func_scan_len):
        t = t_fmri[kk]
        if phys_label == 'ppg':
            t_p = trigger_extract_hr_rate(t, ppg_ibi, ppg_t_ibi, win, tr_func, func_scan_len)
        elif phys_label == 'resp':
            t_p = trigger_extract_resp_var(ts, t, win, tr_physio)
        elif phys_label == 'bpp':
            t_p = trigger_extract_map(ts_dbp, ts_sbp, t, win, tr_physio)

        t_phys.append(t_p)

    return t_phys
    


def trigger_extract_hr_rate(t, IBI, t_IBI, win, func_tr, scan_len, 
                            max_hr=100, min_hr=40):
    # Get bounds of window
    t1 = max(0, t - win*0.5)
    t2 = min(func_tr*scan_len, t + (win*0.5))

    inds1 = [idx for idx, element in enumerate(t_IBI) if element <= t2]
    inds2 = [idx for idx, element in enumerate(t_IBI) if element >= t1]

    inds = set(inds1) & set(inds2)
    
    isEmpty = (len(inds) == 0)
    if isEmpty:
        return np.nan
    else:   
        ids_ele = []
        for i in range(len(IBI)-1):
            if i in inds:
                ids_ele.append(IBI[i])

        t_hr = 60./np.median(ids_ele)
        # Check to ensure HR is within acceptable bounds
        if t_hr > max_hr:
            t_hr = np.nan
        elif t_hr < min_hr:
            t_hr = np.nan
        return t_hr


def trigger_extract_resp_var(resp_raw, t, win, physio_tr):
    n_p = resp_raw.size
    list_index = trigger_window_range(t, n_p, physio_tr, win)
    filtered_rv = resp_raw.loc[resp_raw.index.isin(list_index)]
    return filtered_rv.std()


def trigger_extract_map(dbp_ts, sbp_ts, t, win, physio_tr):
    n_p = dbp_ts.size
    list_index = trigger_window_range(t, n_p, physio_tr, win)
    dbp_ts_filt = dbp_ts.loc[dbp_ts.index.isin(list_index)]
    sbp_ts_filt = dbp_ts.loc[dbp_ts.index.isin(list_index)]
    # Extract mean arterial pressure (MAP) from diasytolic and systolic bp mean
    map_mean = dbp_ts_filt.mean()*(2/3) + sbp_ts_filt.mean()*(1/3)
    return map_mean


def trigger_window_center(tr_func, func_len):
    # Get time stamps of window centered on func TR
    up_limit = tr_func*func_len
    list_t = []
    x = 0
    while x <= up_limit:
        list_t.append(x)
        x += tr_func
    # Create time series index by func TR
    t_fmri = [x + (tr_func/2) for x in list_t]
    return t_fmri


def trigger_window_range(t, ts_max, physio_tr, win):
    # Get physio time stamps within window centered on functional TR
    i1 = max(1, int(np.floor((t - win*0.5)/physio_tr)))
    i2 = min(ts_max, int(np.floor((t + win*0.5)/physio_tr)))
    list_index = range(i1, i2)    
    return list_index




