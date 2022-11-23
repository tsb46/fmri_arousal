import argparse
import json
import os
import neurokit2 as nk
import numpy as np
import pandas as pd

from scipy.stats import zscore
from utils.electrophys_utils import nk_extract_physio, trigger_extract_physio, \
filt_resample_physio_to_func


# Sampling frequency of funcional scans (1/TR)
tr_func = 3.0
sf_func = 1/tr_func
# Number of time samples 
n_scan = 200

# Number of time samples trimmed off original functional scan (4)
n_trim = 4

physio_prefix = ['ppg', 'resp']


def process_physio(physio_signals, sf_physio):
    ## Resp signals
    resp_signals_nk = nk_extract_physio(physio_signals['resp'], 'resp', sf_physio, n_scan)
    resp_signals_window = trigger_extract_physio(physio_signals['resp'], 'resp', 
                                                 n_scan, sf_func, sf_physio)
    
    ## PPG signals
    ppg_signals_nk = nk_extract_physio(physio_signals['ppg'], 'ppg', sf_physio, 
                                       n_scan)
    ppg_signals_window = trigger_extract_physio(physio_signals['ppg'], 'ppg', 
                                                n_scan, sf_func, sf_physio)

    # Create physio df
    physio_list = [ppg_signals_nk['PPG_Rate'].values.tolist(), ppg_signals_window, 
                   ppg_signals_nk['PPG_LOW'].values.tolist(), ppg_signals_nk['PPG_RMS_AMP'].values.tolist(),
                   ppg_signals_nk['PPG_PEAK_AMP'].values.tolist(),resp_signals_nk['RSP_Rate'].values.tolist(), 
                   resp_signals_nk['RSP_Amplitude'].values.tolist(),resp_signals_nk['RSP_RVT'].values.tolist(), 
                   resp_signals_nk['RSP_AMP_HILBERT'].values.tolist(), resp_signals_window,
                   resp_signals_nk['RSP_RVT_IF'].values.tolist(), resp_signals_nk['RSP_RVT_AMP'].values.tolist()]

    physio_labels = ['PPG_HR_NK', 'PPG_HR_W', 'PPG_LOW_NK', 'PPG_RMS_AMP', 'PPG_PEAK_AMP', 
                     'RESP_RATE_NK', 'RESP_AMP_NK', 'RESP_RVT_NK', 'RESP_AMP_HILBERT', 'RESP_VAR_W',
                     'RESP_RVT_IF_NK', 'RESP_RVT_AMP_NK']
    physio_df = pd.DataFrame({label: col for label, col in zip(physio_labels, physio_list)})

    # Forward fill ppg (window) signal 
    physio_df['PPG_HR_W'] = physio_df['PPG_HR_W'].ffill()
    
    return physio_df


def load_physio(subj):
    physio_df = pd.read_csv(subj, compression='gzip', sep='\t', header=None)
    subj_base = subj.rsplit('.tsv')[0]
    physio_json = json.load(open(f'{subj_base}.json'))
    physio_df.columns = physio_json['Columns']
    physio_df.rename(columns = {'cardiac': 'ppg', 'respiratory': 'resp'}, inplace=True)
    sf_physio = physio_json['SamplingFrequency']
    # Trim the beginning of the physio (12sec) to account for the removal of the first four functional volumes
    n_trim_sec = n_trim*tr_func # (TR * num of time samples)
    trim_indx = int(n_trim_sec*sf_physio)
    physio_df_trim = physio_df.iloc[trim_indx:, :].reset_index()
    return physio_df_trim, sf_physio


def run_main(subj, output_physio):
    # Load physio data
    physio_signals, sf_physio = load_physio(subj)
    physio_df = process_physio(physio_signals, sf_physio)
    write_output(physio_df, output_physio)


def write_output(physio_df, output_physio):
    physio_df.to_csv(f'{output_physio}.csv', index=False)
    for col in physio_df.columns:
        np.savetxt(f'{output_physio}_{col}.txt', physio_df[col].values)


if __name__ == '__main__':
    """Preprocess physio data from NKI data"""
    parser = argparse.ArgumentParser(description='Preprocess physio data from HCP data')
    parser.add_argument('-s', '--subject',
                        help='<Required> subject number in Spreng dataset',
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
