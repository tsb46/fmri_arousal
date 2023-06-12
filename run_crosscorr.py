import argparse
import numpy as np
import pandas as pd
import pickle

from scipy.interpolate import interp1d
from scipy.stats import zscore
from utils.load_write import load_data, write_nifti
from utils.glm_utils import xcorr


def cross_corr_maps(func_data, physio, max_lag, n_eval=100):
    # calculate cross-correlation functions across lags and 
    # interpolate with cubic splines
    cc_map_max = []
    cc_map_min = []
    x_interp=np.linspace(-max_lag, max_lag, n_eval)
    for i in range(func_data.shape[1]):
        lags, cc = xcorr(func_data[:,i], physio, max_lag)
        f_interp = interp1d(lags, cc, kind='cubic')
        cc_interp = f_interp(x_interp)
        cc_map_max.append(x_interp[np.argsort(cc_interp)[-1]])
        cc_map_min.append(x_interp[np.argsort(-cc_interp)[-1]])
    return np.array(cc_map_max), np.array(cc_map_min)


def run_crosscorr(dataset, physio, max_lag, out_dir=None):
    # create cross-correlation max/min maps from cross-correlation
    # b/w physio signal and BOLD signals
    # Load data
    func_data, physio_sig, zero_mask, n_vert = load_data(dataset, physio=physio) 
    # squeeze out extra dimension
    physio_sig = np.squeeze(physio_sig)
    # Get cross-correlation peak maps
    cc_map_max, cc_map_min = cross_corr_maps(func_data, physio_sig, max_lag)
    # Write out results
    write_results(dataset, physio, cc_map_max, cc_map_min, zero_mask, n_vert, 
                  out_dir)


def write_results(dataset, term, cc_max, cc_min, zero_mask, n_vert, out_dir):
    # write out results
    if out_dir is not None:
        analysis_str = f'{out_dir}/{dataset}_cc_{term}'
    else:
        analysis_str = f'{dataset}_cc_{term}'
    # Minimum (negative) cross-correlation
    write_nifti(cc_min[np.newaxis, :], f'{analysis_str}_min.nii', zero_mask, n_vert)
    # Maximum (positive cross-correlation)
    write_nifti(cc_max[np.newaxis, :], f'{analysis_str}_max.nii', zero_mask, n_vert)


if __name__ == '__main__':
    """Run main analysis"""
    parser = argparse.ArgumentParser(description='Run Cross-Correlation Lag Analysis')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to run analysis on',
                        choices=['chang', 'chang_bh', 'chang_cue', 
                                 'nki', 'hcp', 'yale', 'spreng', 
                                 'natview'], 
                        required=True,
                        type=str)
    parser.add_argument('-p', '--physio',
                        help='select physio - can only provide one',
                        required=True,
                        choices=['PPG_HR', 'ECG_HR', 'PPG_PEAK_AMP', 'PPG_LOW', 
                                 'RSP_RVT', 'GSR', 'ALPHA', 'THETA', 'DELTA', 
                                 'PUPIL'],
                        type=str) 
    parser.add_argument('-c', '--max_lag',
                        help='max lag to set for evaluating cross-correlation',
                        default=20,
                        required=False,
                        type=int)


    args_dict = vars(parser.parse_args())
    run_crosscorr(args_dict['dataset'], args_dict['physio'], 
                  args_dict['max_lag']
                  )
