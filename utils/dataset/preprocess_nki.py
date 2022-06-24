import argparse
import json
import os
import matplotlib.pyplot as plt
import neurokit2 as nk
import numpy as np
import pandas as pd

from scipy.stats import zscore
from utils.electrophys_utils import nk_extract_physio, trigger_extract_physio, \
filt_resample_physio_to_func


# Sampling frequency of physio files
sf = 62.5
# Sampling frequency of funcional scans (1/TR)
sf_func = 1/1.4
# Number of time samples 
n_scan = 186

physio_prefix = ['ppg', 'resp']


def process_physio(physio_signals):
    ## Resp signals
    resp_signals_nk = nk_extract_physio(physio_signals['resp'], 'resp', sf, 
                                        n_scan, lowpass=0.15)
    resp_signals_window = trigger_extract_physio(physio_signals['resp'], 'resp', 
                                                 n_scan, sf_func, sf)
    
    ## PPG signals
    ppg_signals_nk = nk_extract_physio(physio_signals['ppg'], 'ppg', sf, 
                                       n_scan, lowpass=0.15)
    ppg_signals_window = trigger_extract_physio(physio_signals['ppg'], 'ppg', 
                                                n_scan, sf_func, sf)

    # Create physio df
    physio_list = [ppg_signals_nk['PPG_Rate'].values.tolist(), ppg_signals_window, 
                   ppg_signals_nk['PPG_LOW'].values.tolist(), ppg_signals_nk['PPG_RMS_AMP'].values.tolist(),
                   ppg_signals_nk['PPG_PEAK_AMP'].values.tolist(),resp_signals_nk['RSP_Rate'].values.tolist(), 
                   resp_signals_nk['RSP_Amplitude'].values.tolist(),resp_signals_nk['RSP_RVT'].values.tolist(), 
                   resp_signals_nk['RSP_AMP_HILBERT'].values.tolist(), resp_signals_window]

    physio_labels = ['PPG_HR_NK', 'PPG_HR_W', 'PPG_LOW_NK', 'PPG_RMS_AMP', 'PPG_PEAK_AMP', 
                     'RESP_RATE_NK', 'RESP_AMP_NK', 'RESP_RVT_NK', 'RESP_AMP_HILBERT', 'RESP_VAR_W']
    physio_df = pd.DataFrame({label: col for label, col in zip(physio_labels, physio_list)})

    # Forward fill ppg (window) signal 
    physio_df['PPG_HR_W'] = physio_df['PPG_HR_W'].ffill()
    
    return physio_df


def load_physio(subj):
    physio_df = pd.read_csv(subj, compression='gzip', sep='\t', header=None)
    physio_df = physio_df.dropna(axis=1, how='all')
    subj_base = subj.rsplit('.tsv')[0]
    physio_json = json.load(open(f'{subj_base}.json'))
    physio_df.columns = physio_json['Columns']
    physio_df.rename(columns = {'cardiac': 'ppg', 'respiratory': 'resp'}, inplace=True)
    # Trim the end of the physio to the last functional trigger (5 is trigger indx)
    trim_n = physio_df.loc[physio_df.trigger == 5].tail(1).index[0]
    physio_df_trim = physio_df.iloc[:trim_n, :].copy()
    return physio_df_trim


def run_main(subj, output_physio):
    # Load physio data
    physio_signals = load_physio(subj)
    physio_df = process_physio(physio_signals)
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
    """Preprocess physio data from NKI data"""
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
