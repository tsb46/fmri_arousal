import argparse
import matplotlib.pyplot as plt
import nibabel as nb
import numpy as np
import pandas as pd
import os

from scipy.interpolate import interp1d
from scipy.signal import butter, find_peaks, filtfilt
from scipy.stats import zscore, iqr

"""
May. 2021
Author: Karuna Gujar (modified by Taylor Bolt)
PHYSIO LIBRARY: FMRI Regressors
https://github.com/neurdylab/physio_scripts (private)
"""

def extract_hr(IBI, t_IBI, t, win, func_tr, scan_len):
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

        return 60./np.median(ids_ele)


def extract_rv(resp_raw, t, win, physio_tr):
    n_p = resp_raw.size
    i1 = max(1, int(np.floor((t - win*0.5)/physio_tr)))
    i2 = min(n_p, int(np.floor((t + win*0.5)/physio_tr)))
    list_index = range(i1, i2)
    filtered_rv = resp_raw.loc[resp_raw.index.isin(list_index)]
    return filtered_rv.std()


def get_peaks_ppg(ppg_sig, ppg_rng, fs_physio):
    minHeight = 0.05*ppg_rng
    minDist = (fs_physio*.1)+(fs_physio/2)
    locs, pks_d = find_peaks(ppg_sig, height = minHeight, distance = minDist)
    return locs, pks_d 


def load_physio(file, dataset, resp_col, ppg_col, trigger_col):
    if dataset == 'nki':
        physio_df = pd.read_csv(file, compression='gzip', sep='\t', header=None)
        resp = physio_df.iloc[:,resp_col].copy()
        ppg = physio_df.iloc[:, ppg_col].copy()
        trigger = physio_df.iloc[:, trigger_col].copy()
    elif dataset == 'hcp':
        physio_mat = np.loadtxt(file)
        physio_df = pd.DataFrame(physio_mat)
        resp = physio_df.iloc[:, resp_col].copy()
        ppg = physio_df.iloc[:, ppg_col].copy()
        trigger = physio_df.iloc[:, trigger_col].copy()
    return resp, ppg, trigger


def run_main(dataset, file, resp_col, ppg_col, trigger_col, func_tr, scan_len, 
             fs_physio, output_file, bp_filt=[0.5,2], win=3, max_hr=100, min_hr=40):
    # Get sampling rate for func
    func_fs = 1 / func_tr
    # Get physio time interval
    physio_tr = 1 / fs_physio
    Fn = fs_physio/2  # physio nyquist freq

    # Load physio file
    resp_raw, ppg_raw, trigger = load_physio(file, dataset, resp_col, ppg_col, trigger_col)

    # Filter ppg data
    Wn = [x / Fn for x in bp_filt]
    Nb = 2
    b, a = butter(Nb,Wn, btype='band')
    ppg_bpf =  filtfilt(b,a, ppg_raw)


    # Get IQR OF PPG (25%, 75%)
    ppg_rng = iqr(ppg_bpf)

    # Get peaks of PPG
    locs, pk_d = get_peaks_ppg(ppg_bpf, ppg_rng, fs_physio)
    ppg_samples = locs
    ppg_trig_times = ppg_samples*physio_tr

    # Calculate interbeat interval
    IBI = (np.diff(ppg_trig_times))
    t_IBI = 0.5*(ppg_trig_times[1:] + ppg_trig_times[0:-1])


    up_limit = func_tr*scan_len
    list_t = []
    x = 0
    while x <= up_limit:
        list_t.append(x)
        x += func_tr

    # Create time series index by func TR
    t_fmri = [x + (func_tr/2) for x in list_t]

    hr = []
    rv = []
    # Loop through time points in func and calculate RV and HR within window
    for kk in range(scan_len):
        t = t_fmri[kk]
        # Extract HR within window
        t_hr = extract_hr(IBI, t_IBI, t, win, func_tr, scan_len)
        # Check to ensure HR is within acceptable bounds
        if t_hr > max_hr:
            t_hr = np.nan
        elif t_hr < min_hr:
            t_hr = np.nan

        hr.append(t_hr)
        # Extract RVT within window
        t_rv = extract_rv(resp_raw, t, win, physio_tr)
        if np.isnan(t_rv):
            t_rv = np.nan
        rv.append(t_rv)

    # Interpolate nans in HR and RV time series with cubic spline
    signal_interp = []
    for signal in [hr, rv]:
        signal = np.array(signal)
        not_nan = np.logical_not(np.isnan(signal))
        indices = np.arange(len(signal))
        interp = interp1d(indices[not_nan], signal[not_nan], kind='cubic', 
                          fill_value="extrapolate")
        signal_interp.append(interp(indices))

    hr, rv = signal_interp

    fig, ax = plt.subplots()
    ax.plot(zscore(rv), label='rv')
    ax.plot(zscore(hr), label='hr')
    ax.legend()
    fig.savefig(f'{output_file}.png')
    np.savetxt(f'{output_file}_rv.txt', rv)
    np.savetxt(f'{output_file}_hr.txt', hr)


if __name__ == '__main__':
    """Preprocess respiration and PPG data """
    parser = argparse.ArgumentParser(description='Preprocess PPG and respiration data')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to run analysis on - only option is NKI',
                        choices=['nki', 'hcp'], 
                        required=True,
                        type=str)
    parser.add_argument('-f', '--file_path_bids',
                        help='<Required> path to physio file (in BIDS format)',
                        required=True,
                        type=str)
    parser.add_argument('-r', '--resp',
                        help='column number of respiration time series (starts at 0)',
                        required=True,
                        type=int)
    parser.add_argument('-p', '--ppg',
                        help='column number of ppg time series (starts at 0)',
                        required=True,
                        type=int)
    parser.add_argument('-g', '--trigger',
                        help='column number of trigger time series (starts at 0)',
                        required=True,
                        type=int)
    parser.add_argument('-t', '--tr',
                        help='the repetition time of the data',
                        required=True,
                        type=float)
    parser.add_argument('-n', '--scan_length',
                        help='number of volumes in scan',
                        required=True,
                        type=int)
    parser.add_argument('-s', '--sampling_rate',
                        help='sampling rate of physio data',
                        required=True,
                        type=float)
    parser.add_argument('-o', '--output_file',
                        help='output file path',
                        required=False,
                        default=os.getcwd(),
                        type=str)
    args_dict = vars(parser.parse_args())
    run_main(args_dict['dataset'], args_dict['file_path_bids'], args_dict['resp'], args_dict['ppg'],
             args_dict['trigger'], args_dict['tr'], args_dict['scan_length'], 
             args_dict['sampling_rate'], args_dict['output_file'])
