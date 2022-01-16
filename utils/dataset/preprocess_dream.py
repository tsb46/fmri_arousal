import argparse
import os
import mne
import numpy as np
import pandas as pd
import sys

from scipy.signal import resample
from scipy.stats import zscore
from utils.signal.butterworth_filters import butterworth_filter
from utils.electrophys_utils import compute_rvt, eeg_band_amplitudes, \
extract_ecg_signals, extract_resp_signals, \
extract_rr_interval, extract_rsa, find_interpolate_spikes, \
extract_ppg_signals, extract_emg_signals



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
    ('Beta', 13, 25),
    ('Gamma', 30, 45)
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

# respiration measures
resp_measures = ['VTOT', 'VTH', 'VAB', 'NAF1', 'NAF2P-A1']

# emg measures
emg_measures = ['EMG1', 'EMG2', 'EMG3']


def extract_eeg_bands(eeg):
    eeg_freq = eeg.copy()
    eeg_freq.pick(eeg_chan)
    # Free up memory
    del eeg
    slow_wave_freq = band_freqs.pop(0)
    eeg_bands = eeg_band_amplitudes(eeg_freq, band_freqs)
    slow_wave_band = eeg_band_amplitudes(eeg_freq, [slow_wave_freq], l_trans=0.3, h_trans=0.3)[0]
    eeg_bands.insert(0, slow_wave_band)
    band_freqs.insert(0, slow_wave_freq)
    sleep_events, _ = mne.events_from_annotations(eeg_freq, chunk_duration=5, 
                                                  event_id=annotation_desc_2_event_id, 
                                                  verbose=False)
    tmax = 5 - 1. / eeg_freq.info['sfreq']  # tmax in included
    # Free up memory
    del eeg_freq
    bands_ts = []
    band_labels = [b[0] for b in band_freqs]
    for i, (band, band_label) in enumerate(zip(eeg_bands, band_labels)):
        # Apply low-pass to frequency range of interest
        band.filter(None, 0.15)
        # Resample to 10 Hz
        band, sleep_events_d = band.resample(10, events=sleep_events)
        band_epochs = mne.Epochs(raw=band.copy(), events=sleep_events_d,
                                 tmin=0., tmax=tmax, 
                                 baseline=None, verbose=False)
        eeg_df = band_epochs.to_data_frame()
        if i > 0:
            eeg_df.drop(columns=['time', 'condition'], inplace=True)

        eeg_rename = {e: f'{e}-{band_label}' for e in eeg_chan}
        eeg_df.rename(columns=eeg_rename, inplace=True)
        bands_ts.append(eeg_df)

    bands_df = pd.concat(bands_ts, axis=1)

    return bands_df


def process_physio(physio_dict, eeg_downsample_len):
    # Get ECG signal
    ecg = physio_dict['ECG']
    # Preprocess ECG signals from ECG
    ecg_signals = extract_ecg_signals(ecg, 200)
    # Extract interbeat interval
    rr_interval_ecg = extract_rr_interval(ecg_signals, 200, 'ecg')

    # Preprocess EMG signals
    emg_df_all = []
    for i, measure in enumerate(emg_measures):
        emg_signals = extract_emg_signals(physio_dict[measure], 200)
        # emg_labels = ['EMG_Amplitude', 'EMG_Onsets', 'EMG_Offsets']
        # emg_labels_measure = {label: f'{label}_{measure}' for label in emg_labels}
        # emg_df = emg_signals[emg_labels].copy()
        # emg_df.rename(columns=emg_labels_measure, inplace=True)
        emg_df = pd.DataFrame({measure: emg_signals})
        emg_df_all.append(emg_df)

    # Preprocess Resp signals
    resp_df_all = []
    for i, measure in enumerate(resp_measures):
        resp_signals = extract_resp_signals(physio_dict[measure], 200)
        # Calculate Respiration volume per time
        resp_signals['RVT'] = compute_rvt(resp_signals)
        # Calculate RSA signal
        resp_signals['RSA'] = extract_rsa(ecg_signals, resp_signals)
        # Create one physio df
        resp_labels = ['RVT', 'RSA', 'RSP_Rate', 'RSP_Amplitude']
        resp_labels_measure = {label: f'{label}_{measure}' for label in resp_labels}
        resp_df = resp_signals[resp_labels].copy()
        resp_df.rename(columns = resp_labels_measure, inplace=True)
        resp_df_all.append(resp_df)

    physio_df = pd.concat(resp_df_all + emg_df_all, axis=1)
    physio_df['RR_Interval_ECG'] = rr_interval_ecg
    # Low pass physio signals below 0.15Hz (same as EEG) and downsample to EEG frequency (5Hz)
    physio = []
    for col in physio_df.columns:
        # cleaned_signal = find_interpolate_spikes(physio_df[col].values)
        filtered_signal = butterworth_filter(physio_df[col], None, 
                                             0.15, 200, 'lowpass')
        downsampled_signal = resample(filtered_signal, eeg_downsample_len)
        physio.append(downsampled_signal)

    physio_df = pd.DataFrame(np.array(physio).T, columns = physio_df.columns)
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


def regress_eyeblinks(eeg_raw):
    eeg_proc, _ = mne.preprocessing.regress_artifact(eeg_raw, picks=eeg_chan, 
                                                     picks_artifact=['EOG1', 'EOG2'], 
                                                     betas=None, copy=True, verbose=False)
    return eeg_proc


def run_main(poly_file, hypno_file, output_eeg, output_physio):
    # Load polysomnography and hypynogram
    eeg = load_edf(poly_file)
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
    band_df = extract_eeg_bands(eeg)
    # Preprocess physio data
    physio_df = process_physio(physio_signals, band_df.shape[0])
    write_output(band_df, physio_df, output_eeg, output_physio)


def write_output(band_df, physio_df, output_eeg, output_physio):
    band_df.to_csv(output_eeg, index=False)
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
