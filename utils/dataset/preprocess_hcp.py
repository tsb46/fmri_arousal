import argparse
import os
import matplotlib.pyplot as plt
import neurokit2 as nk
import numpy as np
import pandas as pd

from scipy.stats import zscore
from utils.electrophys_utils import nk_extract_physio, trigger_extract_physio, \
filt_resample_physio_to_func


# Sampling frequency of physio files
sf = 400
sf_resamp=100 # Resample frequency of EEG and physio time courses
# Sampling frequency of funcional scans (1/TR)
sf_func = 1/0.72
physio_prefix = ['ppg', 'resp']


def process_physio(physio_signals, n_scan):
    ## Resp signals
    resp_signals_nk = nk_extract_physio(physio_signals['resp'], 'resp', sf_resamp, 
                                        n_scan)
    resp_signals_window = trigger_extract_physio(physio_signals['resp'], 'resp', 
                                                 n_scan, sf_func, sf_resamp)

    resp_raw = nk.signal_resample(physio_signals['resp'], desired_length=n_scan, method='FFT')
    
    ## PPG signals
    ppg_signals_nk = nk_extract_physio(physio_signals['ppg'], 'ppg', sf_resamp, 
                                       n_scan)
    ppg_signals_window = trigger_extract_physio(physio_signals['ppg'], 'ppg', 
                                                n_scan, sf_func, sf_resamp)


    # Create physio df
    physio_list = [ppg_signals_nk['PPG_Rate'].values.tolist(), ppg_signals_window, 
                   ppg_signals_nk['PPG_LOW'].values.tolist(), ppg_signals_nk['PPG_RMS_AMP'].values.tolist(), 
                   ppg_signals_nk['PPG_PEAK_AMP'].values.tolist(), resp_signals_nk['RSP_Rate'].values.tolist(), 
                   resp_signals_nk['RSP_Amplitude'].values.tolist(),resp_signals_nk['RSP_RVT'].values.tolist(), 
                   resp_signals_nk['RSP_AMP_HILBERT'].values.tolist(), resp_signals_window, 
                   resp_raw.tolist()]

    physio_labels = ['PPG_HR_NK', 'PPG_HR_W', 'PPG_LOW_NK', 'PPG_RMS_AMP', 'PPG_PEAK_AMP', 
                     'RESP_RATE_NK', 'RESP_AMP_NK', 'RESP_RVT_NK', 'RESP_AMP_HILBERT', 
                     'RESP_VAR_W', 'RESP_RAW']
    physio_df = pd.DataFrame({label: col for label, col in zip(physio_labels, physio_list)})

    # Forward fill ppg (window) signal 
    physio_df['PPG_HR_W'] = physio_df['PPG_HR_W'].ffill()
    
    return physio_df


def load_physio(subj):
    physio_signals = np.loadtxt(subj)
    physio_signals_df = pd.DataFrame({ 
     'resp': physio_signals[:,1],
     'ppg': physio_signals[:,2]
     })
    # Downsample to 100Hz (400Hz is unnecessary)
    physio_secs = physio_signals_df.shape[0]/sf
    physio_resamp_n = int(physio_secs * sf_resamp)
    physio_signals_resamp = physio_signals_df.apply(
        nk.signal_resample, desired_length=physio_resamp_n, method='FFT', axis=0
    )
    return physio_signals_resamp


def run_main(subj, output_physio, dataset):
    if dataset == 'rest':
        n_scan = 1200
    elif dataset == 'rel':
        n_scan = 232
    elif dataset == 'wm':
        n_scan = 405

    # Load physio data
    physio_signals = load_physio(subj)
    physio_df = process_physio(physio_signals, n_scan)
    write_output(physio_df, output_physio)


def write_output(physio_df, output_physio):
    physio_df.to_csv(f'{output_physio}.csv', index=False)
    # Save physio plots
    fig, ax = plt.subplots(figsize=(15,5))
    physio_df[['RESP_AMP_HILBERT', 'PPG_HR_NK']].apply(zscore, axis=0).plot(ax=ax)
    plt.savefig(f'{output_physio}.png')

    for col in physio_df.columns:
        np.savetxt(f'{output_physio}_{col}.txt', physio_df[col].values)


if __name__ == '__main__':
    """Preprocess physio data from HCP data"""
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
    parser.add_argument('-d', '--dataset',
                        help='hcp dataset (rest, wm, rel)',
                        required=True,
                        type=str)
    args_dict = vars(parser.parse_args())
    run_main(args_dict['subject'], args_dict['output_file_physio'], 
             args_dict['dataset'])
