import argparse
import os
import mne
import neurokit2 as nk
import numpy as np
import pandas as pd
import sys

from scipy.io import loadmat
from scipy.stats import zscore
from utils.signal.butterworth_filters import butterworth_filter
from utils.electrophys_utils import clip_spikes, eeg_band_amplitudes, \
nk_extract_physio, filt_resample_physio_to_func, \
trigger_extract_physio



# EEG bands of interest
band_freqs = [
    ('Delta', 1, 3),
    ('Theta', 4, 7),
    ('Alpha', 8, 12),
]

# physio signals
physio_labels = ['ECG']

eeg_chan_avg = [
    'P3', 
    'P4', 
    'Pz', 
    'O1',
    'O2',
    'Oz'                   
]

func_tr = 2.1 # Functional TR
sf_func = 1/2.1 # Functional sampling frequency
sf_resamp=100 # Resample frequency of EEG and physio time courses

def construct_mne_obj(eeg_mat, subj, output_mne=None):
    sf = eeg_mat['EEG']['srate'].item()
    n_chan = eeg_mat['EEG']['nbchan'].item()
    chan_labels = [chan[0] for chan in eeg_mat['EEG']['chanlocs'].item()]
    if subj == 'sub_0010':
        chan_types = ['eeg'] * 30 + ['ecg']
    else:
        chan_types = ['eeg'] * 31 + ['ecg']
    info = mne.create_info(chan_labels, ch_types=chan_types, sfreq=sf)
    data = np.vstack(eeg_mat['EEG']['data'].item())
    eeg_raw = mne.io.RawArray(data, info, verbose=False)
    # Resample from 250 to 100Hz
    eeg_resamp = eeg_raw.copy().resample(sf_resamp)
    if output_mne is not None:
        eeg_resamp.save(f'{output_mne}.raw.fif', overwrite=True)
    return eeg_resamp


def compute_vigilance_rms(eeg, sf, window_s=2, alpha_indx=2, delta_indx=0):
    # 2 second window
    win= window_s*sf
    # Freq ranges 
    delta_l, delta_h = band_freqs[delta_indx][1], band_freqs[delta_indx][2]
    alpha_l, alpha_h = band_freqs[alpha_indx][1], band_freqs[alpha_indx][2]
    # Filter to theta range
    eeg_freq = eeg.copy()
    eeg_freq.filter(delta_l, delta_h, l_trans_bandwidth=1, h_trans_bandwidth=1)
    delta = eeg_freq.get_data()
    # Filter to alpha range
    eeg_freq = eeg.copy()
    eeg_freq.filter(alpha_l, alpha_h, l_trans_bandwidth=1, h_trans_bandwidth=1)
    alpha = eeg_freq.get_data()
    # Create temporary dataframe to use Pandas rolling function
    df = pd.DataFrame({'delta': np.mean(delta, axis=0), 'alpha': np.mean(alpha, axis=0)})
    # Calculate vigilance
    alpha_rolling_rms = df.rolling(win)['alpha'].apply(lambda x: np.sqrt(np.sum(x**2)), raw=True)
    delta_rolling_rms = df.rolling(win)['delta'].apply(lambda x: np.sqrt(np.sum(x**2)), raw=True)
    vigilance = alpha_rolling_rms/delta_rolling_rms
    return vigilance.fillna(0)


def extract_eeg_bands(eeg, sf_resamp, func_len):
    eeg_freq = eeg.copy()
    # Pick EEG channels
    eeg_freq.pick(eeg_chan_avg)
    # Free up memory
    del eeg
    # Extract EEG Bands
    eeg_bands = eeg_band_amplitudes(eeg_freq, band_freqs, l_trans='auto', h_trans='auto')
    # Compute vigilance index
    vigilance = compute_vigilance_rms(eeg_freq, sf_resamp)
    vigilance_resamp = filt_resample_physio_to_func(vigilance, high_cut=0.15, func_len=func_len, sf=sf_resamp)
    # Get infraslow component (needs an extra downsampling step before filtering)
    infraslow = filter_infraslow(eeg_freq)
    infraslow_resamp = filt_resample_physio_to_func(infraslow, high_cut=0.15, func_len=func_len, sf=2)
    # Free up memory
    del eeg_freq
    # Run through bands and resample
    bands_ts = []
    band_labels = [b[0] for b in band_freqs]
    for i, (band, band_label) in enumerate(zip(eeg_bands, band_labels)):
        eeg_df = band.to_data_frame()
        eeg_df.drop(columns=['time'], inplace=True)
        eeg_df_resamp = eeg_df.apply(lambda x: filt_resample_physio_to_func(x, high_cut=0.15, func_len=func_len, sf=sf_resamp), axis=0)
        eeg_avg = eeg_df_resamp.mean(axis=1)
        eeg_avg = eeg_avg.rename(band_label)
        bands_ts.append(eeg_avg)

    bands_df = pd.concat(bands_ts, axis=1)
    bands_df['vigilance'] = vigilance_resamp
    bands_df['Infraslow'] = infraslow_resamp
    return bands_df


def filter_infraslow(eeg):
    eeg_low = eeg.copy().resample(2, npad=2000)
    eeg_array = eeg_low.get_data().T
    eeg_filt = []
    for i in range(eeg_array.shape[1]):
        eeg_filt.append(butterworth_filter(eeg_array[:,i], 0.01, 0.1, 2, 'bandpass'))
    eeg_filt = np.stack(eeg_filt, axis=1)
    eeg_infraslow = np.mean(eeg_filt, axis=1)
    return eeg_infraslow


def process_physio(physio_dict, sf_resamp, func_len):
    ## Preprocess ECG signals from ECG
    ecg = physio_dict['ecg'].T
    ecg_signals_nk = nk_extract_physio(ecg, 'ecg', sf_resamp, func_len, lowpass=0.15)
    
    ## Preprocess Resp signals
    resp = physio_dict['resp']
    # neurokit preprocessing
    resp_signals_nk = nk_extract_physio(resp, 'resp', sf_resamp, func_len, lowpass=0.15)
    # trigger-based preprocessing
    resp_signals_window = trigger_extract_physio(resp, 'resp', func_len, sf_func, sf_resamp)

    ## Preprocess PPG signals
    ppg = physio_dict['ppg']
    # neurokit preprocessing
    ppg_signals_nk = nk_extract_physio(ppg, 'ppg', sf_resamp, func_len, lowpass=0.15)
    # trigger-based preprocessing
    ppg_signals_window = trigger_extract_physio(ppg, 'ppg', func_len, sf_func, sf_resamp)


    # Create final physio dataframe
    physio_list = [ecg_signals_nk['ECG_Rate'].values.tolist(), resp_signals_nk['RSP_Rate'].values.tolist(), 
                   resp_signals_nk['RSP_Amplitude'].values.tolist(),resp_signals_nk['RSP_RVT'].values.tolist(), 
                   resp_signals_nk['RSP_AMP_HILBERT'].values.tolist(), resp_signals_nk['RSP_RV'],
                   resp_signals_window, ppg_signals_nk['PPG_Rate'].values.tolist(), ppg_signals_window] 
    physio_labels = ['ECG_HR_NK', 'RESP_RATE_NK', 'RESP_AMP_NK', 'RESP_RVT_NK', 'RESP_AMP_HILBERT', 'RSP_RV', 
                     'RESP_VAR_W', 'PPG_RATE_NK', 'PPG_RATE_W']
    physio_df = pd.DataFrame({label: col for label, col in zip(physio_labels, physio_list)})

    return physio_df


def load_eeg_mat(eeg_mat_fp):
    eeg_raw = loadmat(eeg_mat_fp, squeeze_me=True)
    subj = eeg_mat_fp.split('-')[0].split('/')[-1]
    return eeg_raw, subj


def load_physio_mat(eeg_mat_physio, sf_physio_resamp):
    physio_raw = loadmat(eeg_mat_physio, squeeze_me=True)
    sf_physio = 1/physio_raw['OUT_p']['dt_phys'].item()
    # Pull physio data into dict
    ppg = physio_raw['OUT_p']['card_dat'].item()
    resp = physio_raw['OUT_p']['resp'].item()['wave'].item()
    physio = {'ppg': ppg, 'resp': resp}
    # Downsample to 100Hz (2000Hz is unnecessary)
    physio_secs = len(physio['resp'])/sf_physio
    physio_resamp_n = int(physio_secs * sf_physio_resamp)
    physio = {l: nk.signal_resample(physio[l], desired_length=physio_resamp_n, method='FFT') 
              for l in physio.keys()}
    return physio


def run_main(eeg_mat_fp, physio_mat_fp, func_len, output_mne, output_eeg, output_physio):
    # Load polysomnography and hypynogram
    eeg, subj = load_eeg_mat(eeg_mat_fp)
    # Load raw physio data
    physio = load_physio_mat(physio_mat_fp, sf_resamp)
    # Construct MNE object from eeg .mat file (also pull back sampling frequency)
    eeg = construct_mne_obj(eeg, subj, output_mne)
    # Load physio data
    ecg_data = eeg.get_data(physio_labels)
    physio['ecg'] = ecg_data
    # Extract EEG band amplitudes
    band_df = extract_eeg_bands(eeg, sf_resamp, func_len)
    # Preprocess physio data
    physio_df = process_physio(physio, sf_resamp, func_len)
    write_output(band_df, physio_df, output_eeg, output_physio)


def write_output(band_df, physio_df, output_eeg, output_physio):
    band_df.to_csv(f'{output_eeg}.csv', index=False)
    for col in band_df.columns:
        np.savetxt(f'{output_eeg}_{col}.txt', band_df[col].values)
    physio_df.to_csv(f'{output_physio}.csv', index=False)
    for col in physio_df.columns:
        np.savetxt(f'{output_physio}_{col}.txt', physio_df[col].values)


    

if __name__ == '__main__':
    """Preprocess EEG and physio data from Chang dataset"""
    parser = argparse.ArgumentParser(description='Preprocess EEG and physio data from Chang dataset')
    parser.add_argument('-e', '--eeg_mat_fp',
                        help='<Required> File path to .mat file containing (artifact-corrected) EEGLAB data',
                        required=True,
                        type=str)
    parser.add_argument('-p', '--physio_mat_fp',
                        help='<Required> File path to .mat file containing physio raw data',
                        required=True,
                        type=str)
    parser.add_argument('-f', '--func_len',
                        help='<Required> Number of time samples in functional MRI data',
                        required=True,
                        type=int)
    parser.add_argument('-om', '--output_file_mne',
                        help='output file path for raw mne object',
                        required=False,
                        default=os.getcwd(),
                        type=str)
    parser.add_argument('-oe', '--output_file_eeg',
                        help='output file path for preprocessed eeg',
                        required=False,
                        default=os.getcwd(),
                        type=str)
    parser.add_argument('-op', '--output_file_physio',
                        help='output file path for physio time courses',
                        required=False,
                        default=os.getcwd(),
                        type=str)
    args_dict = vars(parser.parse_args())
    run_main(args_dict['eeg_mat_fp'], args_dict['physio_mat_fp'], args_dict['func_len'],
             args_dict['output_file_mne'], args_dict['output_file_eeg'], 
             args_dict['output_file_physio'])
