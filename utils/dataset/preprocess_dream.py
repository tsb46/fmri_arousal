import argparse
import os
import mne
import numpy as np
import pandas as pd
import sys

# from mne.filter import resample as mne_resample
from scipy.signal import resample, hilbert
from scipy.stats import zscore
from utils.signal.butterworth_filters import butterworth_filter
from utils.electrophys_utils import eeg_band_amplitudes, \
nk_extract_physio, filt_resample_physio_to_func, \
reject_amplitude


# global variables

 # sleep stage to numeric dict
annotation_desc_2_event_id = {
 '0': 5,  # sleep stage s4
 '1': 4,  # sleep stage s3
 '2': 3,  # sleep stage s2
 '3': 2,  # sleep stage s1
 '4': 6,  # REM
 '5': 1  # Wake
 }

# EEG bands of interest
band_freqs = [
    ('Slow_Delta', 0.5, 2),
    ('Delta', 1, 3),
    ('Theta', 4, 7),
    ('Alpha', 8, 12),
]

# eeg channels
eeg_chan = [
    'FP1-A2',
    'CZ-A1',
    'O1-A2'
]

# channel mapping dict
type_dict = {
    'FP1-A2': 'eeg',
    'CZ-A1': 'eeg',
    'O1-A2': 'eeg',
    'FP2-A1': 'eeg',
    'O2-A1': 'eeg',
    'CZ2-A1': 'eeg',
    'ECG': 'ecg',
    'EOG1': 'eog',
    'EOG2': 'eog',
    'EMG1': 'emg',
    'EMG2': 'emg',
    'EMG3': 'emg',
    'VTOT': 'resp',
    'NAF2P-A1': 'resp',
    'NAF1': 'resp',
    'VTH': 'resp',
    'VAB': 'resp',
    'SAO2': 'misc'
}

# Channel rename mapping for consistency across subjects
ch_rename = {
    'EOG1-A2': 'EOG1',
    'EOG2-A2': 'EOG2'
}

# physio signals
physio_labels = ['ECG', 'VTOT', 'NAF2P-A1',
                 'NAF1', 'VTH', 'VAB', 'SAO2', 
                 'EMG1', 'EMG2', 'EMG3']

# EMG measures
emg_cols = ['EMG1', 'EMG2', 'EMG3']

# respiration measures
resp_measures = ['VTOT', 'VTH', 'VAB', 'NAF1', 'NAF2P-A1']

# emg measures
emg_measures = ['EMG1', 'EMG2', 'EMG3']

# Sampling frequency of EEG
sf_eeg = 200

# Resample frequency
sf_resample = 30


def compute_vigilance_index(eeg_df, low_rng=(1,7), high_rng=(8,12)):
    low_amp_all = []
    high_amp_all = []
    for chan in eeg_chan:
        low_sig = butterworth_filter(eeg_df[chan], low_rng[0], low_rng[1], sf_resample, 'bandpass')
        low_amp = np.abs(hilbert(low_sig))
        low_amp_all.append(low_amp)

        high_sig = butterworth_filter(eeg_df[chan], high_rng[0], high_rng[1], sf_resample, 'bandpass')
        high_amp = np.abs(hilbert(high_sig))
        high_amp_all.append(high_amp)
    low_amp_avg = np.mean(low_amp_all, axis=0)
    high_amp_avg = np.mean(high_amp_all, axis=0)

    return high_amp_avg/low_amp_avg



def eeg_preprocess(eeg):
    eeg_sleep = eeg.copy()
    eeg_sleep.pick(eeg_chan)
    # Free up memory
    del eeg
    # Create sleep stage annotations
    sleep_events, _ = mne.events_from_annotations(eeg_sleep, chunk_duration=5, 
                                                  event_id=annotation_desc_2_event_id, 
                                                  verbose=False)
    tmax = 5 - 1. / eeg_sleep.info['sfreq']  # tmax in included
    # Identify high movement epochs
    print('identify high motion epochs...')
    reject_indx = reject_high_motion_epochs(eeg_sleep.copy(), sleep_events, tmax)

    # Epoch eeg based on sleep stages
    eeg_sleep, sleep_events_d = eeg_sleep.resample(sf_resample, events=sleep_events)
    eeg_sleep = mne.Epochs(raw=eeg_sleep.copy(), events=sleep_events_d,
                             tmin=0., tmax=tmax, reject=None,
                             baseline=None, verbose=False)
    eeg_df = eeg_sleep.to_data_frame()
    eeg_df = eeg_df.merge(reject_indx, right_on='epoch', left_on='epoch', how='left')

    return eeg_df


def process_physio(physio_dict, eeg_downsample_len):

    ## ECG signal
    # Neurokit physio signal extraction
    ecg_signals = nk_extract_physio(physio_dict['ECG'], 'ecg', sf_eeg, 
                                    eeg_downsample_len, lowpass=None)

    # Preprocess EMG signals
    emg_df_all = []
    for i, measure in enumerate(emg_measures):
        emg_signals = nk_extract_physio(physio_dict[measure], 'emg', sf_eeg, 
                                        eeg_downsample_len, lowpass=None)        
        emg_signals.rename(columns={'EMG_AMP': measure}, inplace=True)
        emg_df_all.append(emg_signals)

    emg_df_all = pd.concat(emg_df_all, axis=1)

    # Preprocess Resp signals
    resp_df_all = []
    for i, measure in enumerate(resp_measures):
        resp_signals = nk_extract_physio(physio_dict[measure], 'resp', sf_eeg, 
                                         eeg_downsample_len, lowpass=None)
        cols_select = ['RSP_Amplitude', 'RSP_Clean', 'RSP_Rate', 'RSP_RVT', 'RSP_AMP_HILBERT', 'RSP_RV']
        resp_signals = resp_signals[cols_select].copy()
        cols_m = {col: f'{col}_{measure}' for col in resp_signals.columns}
        resp_signals.rename(columns=cols_m, inplace=True)
        resp_df_all.append(resp_signals)

    resp_df_all = pd.concat(resp_df_all, axis=1)

    physio_df = pd.concat([resp_df_all, emg_df_all, ecg_signals], axis=1)

    return physio_df


def load_edf(poly_file):
    eeg_raw = mne.io.read_raw_edf(poly_file, verbose=False)
    return eeg_raw


def load_hypnogram(hypno_file):
    # Load sleep stage .csv into pandas dataframe
    hypno = np.loadtxt(hypno_file, skiprows=1)
    hypno = hypno.astype(int)
    onsets = np.arange(0, hypno.size*5, 5)
    # Create MNE annotation object
    annot = mne.Annotations(onsets, 5, hypno)
    return annot


def reject_high_motion_epochs(eeg, sleep_events, tmax, epoch_resample=100):
    # Epoch eeg based on sleep stages
    eeg_epoch, sleep_events_d = eeg.resample(epoch_resample, events=sleep_events)
    eeg_epoch = mne.Epochs(raw=eeg_epoch.copy(), events=sleep_events_d,
                             tmin=0., tmax=tmax, reject=None,
                             baseline=None, verbose=False)
    epoch_df = eeg_epoch.to_data_frame()
    eeg_ts = [epoch_df[chan] for chan in eeg_chan]
    epoch_df['reject'] = reject_amplitude(eeg_ts, epoch_df['epoch'], epoch_resample)
    reject_indx = epoch_df.groupby('epoch')['reject'].apply(lambda x: (x > 0).any()).reset_index()
    return reject_indx


def regress_eyeblinks(eeg_raw):
    eeg_proc, _ = mne.preprocessing.regress_artifact(eeg_raw, picks=eeg_chan, 
                                                     picks_artifact=['EOG1', 'EOG2'], 
                                                     betas=None, copy=True, verbose=False)
    return eeg_proc


def run_main(poly_file, hypno_file, output_eeg, output_physio):
    # Load polysomnography 
    print('load eeg data...')
    eeg = load_edf(poly_file)
    # Load Hypnogram
    annot = load_hypnogram(hypno_file)
    # Set annotations
    eeg.set_annotations(annot)
    # regress out eyeblinks
    eeg.load_data()
    eeg = regress_eyeblinks(eeg)
    # Load physio data
    physio_data = eeg.get_data(picks=physio_labels)
    physio_signals = {phys: physio_data[i,:] for i, phys in enumerate(physio_labels)}
    # Free up memory
    del physio_data
    # Extract EEG band amplitudes
    print('preprocess eeg data...')
    eeg_df = eeg_preprocess(eeg)
    print('compute vigilance index')
    eeg_df['vigilance'] = compute_vigilance_index(eeg_df)
    # Preprocess physio data
    print('preprocess physio data...')
    physio_df = process_physio(physio_signals, eeg_df.shape[0])
    write_output(eeg_df, physio_df, output_eeg, output_physio)


def write_output(eeg_df, physio_df, output_eeg, output_physio):
    eeg_df.to_csv(output_eeg, index=False)
    physio_df.to_csv(f'{output_physio}', index=False)

    

if __name__ == '__main__':
    """Preprocess polysomnography data from EDF data"""
    parser = argparse.ArgumentParser(description='Preprocess polysomnography data from EDF data')
    parser.add_argument('-p', '--polysomnography',
                        help='<Required> File path to polysomnography file in edf dataset',
                        required=True,
                        type=str)
    parser.add_argument('-g', '--hypnogram',
                        help='<Required> File path to hynogram file in edf dataset',
                        required=True,
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
    run_main(args_dict['polysomnography'], args_dict['hypnogram'], 
             args_dict['output_file_eeg'], args_dict['output_file_physio'])
