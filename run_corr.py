import argparse
import numpy as np
import pandas as pd
import pickle

from scipy.interpolate import interp1d
from scipy.stats import zscore
from utils.load_write import load_data, write_nifti
from utils.glm_utils import xcorr


def cross_corr_maps(func_data, physio, max_lag, n_samples):
    # calculate cross-correlation functions across lags and 
    # interpolate with cubic splines
    cc_map = []
    x_interp=np.linspace(-max_lag, max_lag, n_samples)
    for i in range(func_data.shape[1]):
        lags, cc = xcorr(func_data[:,i], physio, max_lag)
        f_interp = interp1d(lags, cc, kind='cubic')
        cc_interp = f_interp(x_interp)
        cc_map.append(cc_interp)
    return np.array(cc_map).T


def run_crosscorr(dataset, physio, max_lag, n_samples, m_param, out_dir=None):
    # create cross-correlation max/min maps from cross-correlation
    # b/w physio signal and BOLD signals
    # Load data
    func_data, physio_sig, params = load_data(dataset, physio=physio, 
                                              multiecho=m_param) 
    # squeeze out extra dimension
    physio_sig = np.squeeze(physio_sig)
    # Get cross-correlation peak maps
    cc_map = cross_corr_maps(func_data, physio_sig, 
                             max_lag, n_samples)
    # Write out results
    write_results(dataset, physio, cc_map, params, out_dir)


def write_results(dataset, term, cc_map, params, out_dir):
    # write out results
    if out_dir is not None:
        analysis_str = f'{out_dir}/{dataset}_cc_{term}'
    else:
        analysis_str = f'{dataset}_cc_{term}'
    # Write cross-correlation lag maps
    write_nifti(cc_map, f'{analysis_str}.nii', params)


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
                        help='max lag to set for evaluating cross-correlation (in volumes)',
                        default=10,
                        required=False,
                        type=int)
    parser.add_argument('-n', '--n_samp',
                        help='cross correlations are cubic interpolated. Specify number of '
                        'equally spaced interpolation points between -max_lag and max_lag',
                        default=20,
                        required=False,
                        type=int)
    parser.add_argument('-m', '--m_param',
                        help='For multiecho data, specify multiecho parameter - kappa or rho',
                        choices=['kappa', 'rho'], 
                        required=False,
                        default=None,
                        type=str)

    args_dict = vars(parser.parse_args())
    run_crosscorr(args_dict['dataset'], args_dict['physio'], 
                  args_dict['max_lag'], args_dict['n_samp'],
                  args_dict['m_param']
                  )
