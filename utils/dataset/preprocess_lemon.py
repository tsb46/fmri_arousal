import argparse
import os
import mne
import numpy as np
import pandas as pd
import sys

from scipy.signal import resample
from scipy.stats import zscore
from utils.signal.butterworth_filters import butterworth_filter
from utils.electrophys_utils import nk_extract_physio, trigger_extract_physio, \
filt_resample_physio_to_func


# Sampling frequency of physio files
sf = 1000
# Sampling frequency of funcional scans (1/TR)
sf_func = 1/1.4
# Number of original time samples (before trim)
n_orig = 657
# Number of time samples (PPG and BP)
n_sub = 470

physio_prefix = ['bpp', 'ecg', 'ppg', 'resp']

subj_physio_str = lambda x, y: f'data/dataset_lemon/physio/raw/{x}_task_rest_physio_{y}.tsv.gz'


def process_physio(physio_signals):
    ## ECG signal
    # Neurokit physio signal extraction
    ecg_signals_nk = nk_extract_physio(physio_signals['ecg'][0], 'ecg', sf, 
                                       n_orig, lowpass=0.1)

    ## BP signals
    bp_dbp, bp_sbp = physio_signals['bpp'][0], physio_signals['bpp'][1]
    # Neurokit physio signal extraction
    map_signals_nk = nk_extract_physio([bp_dbp, bp_sbp], 'bpp', sf, 
                                       n_sub, lowpass=0.1)
    # window based physio signal extraction
    map_signals_window = trigger_extract_physio([bp_dbp, bp_sbp], 'bpp', 
                                                n_sub, sf_func, sf)
    
    ## Resp signals
    resp_signals_nk = nk_extract_physio(physio_signals['resp'][0], 'resp', sf, 
                                        n_orig, lowpass=0.1)
    resp_signals_window = trigger_extract_physio(physio_signals['resp'][0], 'resp', 
                                                 n_orig, sf_func, sf)
    
    ## PPG signals
    ppg_signals_nk = nk_extract_physio(physio_signals['ppg'][0], 'ppg', sf, 
                                       n_sub, lowpass=0.1)
    ppg_signals_window = trigger_extract_physio(physio_signals['ppg'][0], 'ppg', 
                                                n_sub, sf_func, sf)

    # Create df of PPG and MAP
    ppg_map_list = [map_signals_nk['MAP'].values.tolist(), map_signals_window, 
                    ppg_signals_nk['PPG_Rate'].values.tolist(), ppg_signals_window]
    ppg_map_labels = ['MAP_NK', 'MAP_W', 'PPG_HR_NK', 'PPG_HR_W']
    ppg_map_df = pd.DataFrame({label: col for label, col in zip(ppg_map_labels, ppg_map_list)})
    # Forward fill ppg (window) signal 
    ppg_map_df['PPG_HR_W'] = ppg_map_df['PPG_HR_W'].ffill()
    # Create df of ECG and Resp
    ecg_resp_list = [ecg_signals_nk['ECG_Rate'].values.tolist(), resp_signals_nk['RSP_Rate'].values.tolist(), 
                     resp_signals_nk['RSP_Amplitude'].values.tolist(),resp_signals_nk['RSP_RVT'].values.tolist(), 
                     resp_signals_nk['RSP_AMP_HILBERT'].values.tolist(), resp_signals_window] 
    ecg_resp_labels = ['ECG_HR_NK', 'RESP_RATE_NK', 'RESP_AMP_NK', 'RESP_RVT_NK', 'RESP_AMP_HILBERT', 'RESP_VAR_W']
    ecg_resp_df = pd.DataFrame({label: col for label, col in zip(ecg_resp_labels, ecg_resp_list)})
    return ppg_map_df, ecg_resp_df


def load_physio(subj):
    physio_signals = {}
    for phys_label in physio_prefix: 
        fp = subj_physio_str(subj, phys_label)
        phys_df = pd.read_csv(fp, compression='gzip', delimiter='\t', header=None)
        physio_signals[phys_label] = phys_df
    return physio_signals


def run_main(subj, trim, output_physio):
    # Load physio data
    physio_signals = load_physio(subj)
    ppg_map_df, ecg_resp_df = process_physio(physio_signals)
    write_output(ppg_map_df, ecg_resp_df, output_physio, trim)


def write_output(ppg_map_df, ecg_resp_df, output_physio, trim):
    if trim:
        ecg_resp_df = ecg_resp_df.iloc[:n_sub, :]

    ppg_map_df.to_csv(f'{output_physio}_map_ppg.csv', index=False)
    ecg_resp_df.to_csv(f'{output_physio}_ecg_resp.csv', index=False)
    for col in ppg_map_df.columns:
        np.savetxt(f'{output_physio}_{col}.txt', ppg_map_df[col].values)
    for col in ecg_resp_df.columns:
        np.savetxt(f'{output_physio}_{col}.txt', ecg_resp_df[col].values)


if __name__ == '__main__':
    """Preprocess physio data from LEMON data"""
    parser = argparse.ArgumentParser(description='Preprocess physio data from LEMON data')
    parser.add_argument('-s', '--subject',
                        help='<Required> subject number in LEMON dataset',
                        required=True,
                        type=str)
    parser.add_argument('-t', '--trim_data',
                        help='whether to trim ECG and RESP data to PPG and BP data length (470)',
                        required=False,
                        default=1,
                        type=int)
    parser.add_argument('-o', '--output_file_physio',
                        help='output file path for physio time courses',
                        required=False,
                        default=os.getcwd(),
                        type=str)
    args_dict = vars(parser.parse_args())
    run_main(args_dict['subject'], args_dict['trim_data'], 
             args_dict['output_file_physio'])
