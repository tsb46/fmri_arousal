import argparse
import os
import mne
import numpy as np
import pandas as pd

from neurokit2 import signal_findpeaks, rsp_peaks, rsp_amplitude, rsp_rate
from scipy.signal import detrend, hilbert
from scipy.stats import zscore
from utils.signal.butterworth_filters import butterworth_filter

# Global variables
fs = 3 # eke sampling rate

def load_data(data_file, condition_file):
    # Load data file
    df = pd.read_csv(data_file.replace('\r', ''), encoding="latin1")
    df.columns = [col_name.replace('"', '') for col_name in df.columns]
    # Encoding of csv pulls a lot of extra columns (data columns end at 53)
    df = df.iloc[:,:53].copy()
    # Load condition file
    df_cond = pd.read_csv(condition_file.replace('\r', ''), encoding='latin1', header=None)
    return df, df_cond


def preprocess_capno(resp_raw, lowcut=0.05, highcut=0.4):
    resp_filt = butterworth_filter(resp_raw, lowcut, highcut, fs, 'bandpass')
    resp_peaks, info = rsp_peaks(resp_filt, sampling_rate=fs) 
    resp_amp = rsp_amplitude(resp_filt, resp_peaks)
    resp_rate = rsp_rate(resp_filt, sampling_rate=3)
    resp_amp_hilbert = np.abs(hilbert(resp_filt))
    return resp_amp, resp_rate, resp_amp_hilbert


def run_main(data_file, condition_file, output_data):
    # Load preprocessed time series data
    df, df_cond = load_data(data_file, condition_file)
    time_indx = df.pop('Time')
    # Linear detrend physio data
    df = df.apply(lambda x: detrend(x), axis=0)
    resp_amp, resp_rate, resp_amp_hilbert = preprocess_capno(df['resp_uncalibrated'])
    df['resp_amp'] = resp_amp
    df['resp_rate'] = resp_rate
    df['resp_amp_hilbert'] = resp_amp_hilbert
    # Butterworth low-pass filter
    df = df.apply(lambda x: butterworth_filter(x, None, 0.15, fs, 'lowpass'), axis=0)
    # subset to resting-state period (~ first 30min)
    df['Time'] = time_indx.copy()
    rest_end_t = df_cond.loc[df_cond[0] == 1, 1].values[0]
    df_rest = df.loc[df.Time < rest_end_t]
    write_output(df_rest, output_data)


def write_output(df_preprocessed, output_data):
    df_preprocessed.to_csv(output_data.replace('\r', ''), index=False)


if __name__ == '__main__':
    """Preprocess data from EKE data"""
    parser = argparse.ArgumentParser(description='Preprocess data from EKE data')
    parser.add_argument('-d', '--data_file',
                        help='<Required> File path to data file in eke dataset',
                        required=True,
                        type=str)
    parser.add_argument('-c', '--condition_file',
                        help='<Required> File path to condition timing file in eke dataset',
                        required=True,
                        type=str)
    parser.add_argument('-o', '--output_file',
                        help='output file path for preprocessed time courses',
                        required=True,
                        default=os.getcwd(),
                        type=str)
    args_dict = vars(parser.parse_args())
    run_main(args_dict['data_file'], args_dict['condition_file'], args_dict['output_file'])
