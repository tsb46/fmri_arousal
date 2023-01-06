import argparse
import numpy as np
import pandas as pd
import pickle

from scipy.interpolate import interp1d
from scipy.stats import zscore
from utils.load_write import load_data, write_nifti
from utils.glm_utils import xcorr


def cross_corr_maps(func_data, physio, max_lags, n_eval=100):
    cc_map_max = []
    cc_map_min = []
    x_interp=np.linspace(-max_lags, max_lags, n_eval)
    for i in range(func_data.shape[1]):
        lags, cc = xcorr(func_data[:,i], physio, max_lags)
        f_interp = interp1d(lags, cc, kind='cubic')
        cc_interp = f_interp(x_interp)
        cc_map_max.append(x_interp[np.argsort(cc_interp)[-1]])
        cc_map_min.append(x_interp[np.argsort(-cc_interp)[-1]])
    return np.array(cc_map_max), np.array(cc_map_min)


def run_main(dataset, physio, cross_corr_max_lag, regress_global_sig):
    # Load data
    func_data, physio_sig, physio_label, zero_mask, n_vert, analysis_params = \
    load_data(dataset, physio=[physio], load_physio=True, regress_global=regress_global_sig) 

    physio_data = np.squeeze(np.stack(physio_sig, axis=1))
    print('cross corr')
    cc_map_max, cc_map_min = cross_corr_maps(func_data, physio_data, cross_corr_max_lag)
    write_results(dataset, physio, cc_map_max, cc_map_min, zero_mask, n_vert, analysis_params)

def write_results(dataset, term, cc_max, cc_min, zero_mask, n_vert, analysis_params):
    analysis_str = f'{dataset}_cc_{term}_group'
    # Minimum (negative) cross-correlation
    write_nifti(cc_min[np.newaxis, :], f'{analysis_str}_min.nii', zero_mask, n_vert, analysis_params['mask'])
    # Maximum (positive cross-correlation)
    write_nifti(cc_max[np.newaxis, :], f'{analysis_str}_max.nii', zero_mask, n_vert, analysis_params['mask'])


if __name__ == '__main__':
    """Run main analysis"""
    parser = argparse.ArgumentParser(description='Run Cross-Correlation Lag Analysis')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to run analysis on',
                        choices=['chang', 'nki', 'hcp', 'yale', 'spreng'], 
                        required=True,
                        type=str)
    parser.add_argument('-p', '--physio',
                        help='select physio - can only provide one',
                        required=True,
                        choices=['hr', 'rv', 'alpha', 'delta', 'ppg_low', 
                                 'vigilance', 'pupil'],
                        type=str) 
    parser.add_argument('-c', '--cross_corr_max_lag',
                        help='max lag to set for evaluating cross-correlation',
                        default=20,
                        required=False,
                        type=int)
    parser.add_argument('-g', '--regress_global_sig',
                        help='Whether to regress out global signal from functional data',
                        default=0,
                        required=False,
                        type=int)


    args_dict = vars(parser.parse_args())
    run_main(args_dict['dataset'], args_dict['physio'], 
             args_dict['cross_corr_max_lag'],
             args_dict['regress_global_sig']
             )
