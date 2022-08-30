import argparse
import numpy as np
import pandas as pd
import pickle

from patsy import dmatrix
from scipy.interpolate import interp1d
from scipy.stats import zscore
from scipy.signal import hilbert
from utils.load_write import load_data, write_nifti
from utils.glm_utils import linear_regression, construct_lag_splines, xcorr


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

def evaluate_model(lag_vec, model, spline_basis, var_eval):
    # Create basis from model evaluation using previously defined design matrix (for model fit)
    basis_lag_pred = dmatrix(spline_basis.design_info, {'x': lag_vec}, return_type='dataframe')
     # Intialize basis matrix
    pred_list = [var_eval * basis_lag_pred.iloc[:, l].values for l in range(basis_lag_pred.shape[1])]
    pred_mat = np.vstack(pred_list).T
    # Get predictions from model
    pred_lags = model.predict(pred_mat)
    return pred_lags


def run_main(dataset, physio, p_n_lags, n_n_lags, nknots, 
             regress_global_sig, n_eval, cross_corr, 
             cross_corr_max_lag, var_eval=1):
    # Load data
    func_data, physio_sig, physio_label, zero_mask, n_vert, analysis_params = \
    load_data(dataset, 'group', physio=[physio], load_physio=True, regress_global=regress_global_sig) 

    physio_data = np.squeeze(np.stack(physio_sig, axis=1))
    if cross_corr:
        print('cross corr')
        cc_map_max, cc_map_min = cross_corr_maps(func_data, physio_data, cross_corr_max_lag)
        pred_maps = None
        pred_lag_vec = None
        write_results(dataset, physio, pred_maps, pred_lag_vec, cc_map_max, cc_map_min,
                  zero_mask, n_vert, analysis_params)
    else:
        physio_df = pd.DataFrame(physio_data, columns=physio_label)
        # Construct Design matrix using patsy style formula
        print('construct spline matrix')
        design_mat, spline_basis = construct_lag_splines(physio_df, p_n_lags, n_n_lags, nknots)
        # Lag introduces null values - trim beginning of predictor matrix
        na_indx = ~(np.isnan(design_mat).any(axis=1))
        func_data = func_data[na_indx, :]
        design_mat = design_mat[na_indx, :]
        print('run regression')
        lin_reg = linear_regression(design_mat, func_data, return_model=True, 
                                    intercept=False, norm=False)
        # Get predicted maps at all time lags
        print('get predicted respone at time lags')
        pred_lag_vec = np.linspace(-n_n_lags, p_n_lags, n_eval)
        pred_maps = evaluate_model(pred_lag_vec, lin_reg, spline_basis, var_eval)

        cc_map_max = None
        cc_map_min = None
        # Write out results
        write_results(dataset, physio, pred_maps, pred_lag_vec, cc_map_max, cc_map_min,
                      zero_mask, n_vert, analysis_params)

def write_results(dataset, term, pred_maps, pred_lag_vec, cc_max, cc_min, zero_mask, n_vert, analysis_params):
    analysis_str = f'{dataset}_physio_reg_group_{term}'
    # if cross corr maps specified, write out cross corr maps
    if cc_max is not None:
        # Minimum (negative) cross-correlation
        write_nifti(cc_min[np.newaxis, :], f'{analysis_str}_min_cc.nii', zero_mask, n_vert, analysis_params['mask'])
        # Maximum (positive cross-correlation)
        write_nifti(cc_max[np.newaxis, :], f'{analysis_str}_max_cc.nii', zero_mask, n_vert, analysis_params['mask'])
    else:
        write_nifti(pred_maps, analysis_str, zero_mask, n_vert, analysis_params['mask'])



if __name__ == '__main__':
    """Run main analysis"""
    parser = argparse.ArgumentParser(description='Run Physio GLM w/ Time-lag Spline Regressors')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to run analysis on',
                        choices=['chang', 'nki', 'hcp', 'hcp_fix', 
                                 'hcp_rel', 'hcp_wm', 'spreng'], 
                        required=True,
                        type=str)
    parser.add_argument('-p', '--physio',
                        help='select physio - can only provide one',
                        required=True,
                        choices=['hr', 'rv', 'alpha', 'delta', 'infraslow', 'ppg_low'],
                        type=str)
    parser.add_argument('-pl', '--p_nlags',
                        help='Number of lags (TRs) of physio signal in the positive (forward) direction',
                        required=False,
                        default=0,
                        type=int) 
    parser.add_argument('-nl', '--n_nlags',
                        help='Number of lags (TRs) of physio signal in the negative (backward) direction',
                        required=False,
                        default=0,
                        type=int)   
    parser.add_argument('-k', '--nknots',
                        help='Number of knots in spline basis',
                        default=3, 
                        required=False,
                        type=int)    
    parser.add_argument('-g', '--regress_global_sig',
                        help='Whether to regress out global signal from functional data',
                        default=0,
                        required=False,
                        type=int)
    parser.add_argument('-e', '--n_eval',
                        help='Number of equally spaced samples to evaluate model response between time lags',
                        default=20, 
                        required=False,
                        type=int) 
    parser.add_argument('-t', '--cross_corr_maps',
                        help='Whether to create cross correlation maps between voxel time series and regressor',
                        default=0,
                        required=False,
                        type=int)
    parser.add_argument('-cl', '--cross_corr_max_lag',
                        help='max lag to set for evaluating cross-correlation',
                        default=20,
                        required=False,
                        type=int)

    args_dict = vars(parser.parse_args())
    run_main(args_dict['dataset'], args_dict['physio'], 
             args_dict['p_nlags'], args_dict['n_nlags'], 
             args_dict['nknots'], args_dict['regress_global_sig'], 
             args_dict['n_eval'], args_dict['cross_corr_maps'], 
             args_dict['cross_corr_max_lag'])
