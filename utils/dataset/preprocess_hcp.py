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
sf = 400
# Sampling frequency of funcional scans (1/TR)
sf_func = 1/0.72
# Number of time samples 
n_scan = 1200

physio_prefix = ['ppg', 'resp']


def process_physio(physio_signals):
    ## Resp signals
    resp_signals_nk = nk_extract_physio(physio_signals[:, 1], 'resp', sf, 
                                        n_scan, lowpass=0.1)
    resp_signals_window = trigger_extract_physio(physio_signals[:, 1], 'resp', 
                                                 n_scan, sf_func, sf)
    
    ## PPG signals
    ppg_signals_nk = nk_extract_physio(physio_signals[:, 2], 'ppg', sf, 
                                       n_scan, lowpass=0.1)
    ppg_signals_window = trigger_extract_physio(physio_signals[:, 2], 'ppg', 
                                                n_scan, sf_func, sf)

    # Create physio df
    physio_list = [ppg_signals_nk['PPG_Rate'].values.tolist(), ppg_signals_window, 
                   resp_signals_nk['RSP_Rate'].values.tolist(), 
                   resp_signals_nk['RSP_Amplitude'].values.tolist(),resp_signals_nk['RSP_RVT'].values.tolist(), 
                   resp_signals_nk['RSP_AMP_HILBERT'].values.tolist(), resp_signals_window]

    physio_labels = ['PPG_HR_NK', 'PPG_HR_W', 'RESP_RATE_NK', 'RESP_AMP_NK', 
                     'RESP_RVT_NK', 'RESP_AMP_HILBERT', 'RESP_VAR_W']
    physio_df = pd.DataFrame({label: col for label, col in zip(physio_labels, physio_list)})

    # Forward fill ppg (window) signal 
    physio_df['PPG_HR_W'] = physio_df['PPG_HR_W'].ffill()
    
    return physio_df


def load_physio(subj):
    physio_signals = np.loadtxt(subj)
    return physio_signals


def run_main(subj, output_physio):
    # Load physio data
    physio_signals = load_physio(subj)
    physio_df = process_physio(physio_signals)
    write_output(physio_df, output_physio)


def write_output(physio_df, output_physio):
    physio_df.to_csv(f'{output_physio}.csv', index=False)
    for col in physio_df.columns:
        np.savetxt(f'{output_physio}_{col}.txt', physio_df[col].values)


if __name__ == '__main__':
    """Preprocess physio data from LEMON data"""
    parser = argparse.ArgumentParser(description='Preprocess physio data from HCP data')
    parser.add_argument('-s', '--subject',
                        help='<Required> subject number in HCP dataset',
                        required=True,
                        type=str)
    parser.add_argument('-o', '--output_file_physio',
                        help='output file path for physio time courses',
                        required=False,
                        default=os.getcwd(),
                        type=str)
    args_dict = vars(parser.parse_args())
    run_main(args_dict['subject'], 
             args_dict['output_file_physio'])
