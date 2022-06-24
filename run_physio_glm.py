import argparse
import numpy as np
import fbpca
import pandas as pd
import pickle

from numpy.linalg import pinv
from patsy import dmatrix
from run_pca import pca, rotation
from scipy.stats import zscore
from scipy.signal import hilbert
from utils.load_write import load_data, write_nifti
from utils.glm_utils import linear_regression


def construct_lag_splines(physio_var, p_nlags, n_n_lags, nknots):
    # Define model formula
    # Create lag sequence array (include lag of 0!)
    seq_lag = np.arange(-n_n_lags, p_nlags+1)
    # Create lag splines
    lag_splines = dmatrix("cr(x, df=nknots) - 1", 
                          {"x": seq_lag}, return_type='dataframe')
    # Create design matrix 
    basis_lag = np.zeros((len(physio_var), lag_splines.shape[1]))
    # Loop through lag bases and multiply column pairs
    lag_mat = pd.concat([physio_var.shift(l) for l in seq_lag], axis=1)
    for l in np.arange(lag_splines.shape[1]):
        basis_lag[:, l] = np.dot(lag_mat.values, lag_splines.iloc[:,l].values)
    return basis_lag, lag_splines


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
             regress_global_sig, n_eval, time_lag_maps, var_eval=1):
    # Load data
    func_data, physio_sig, physio_label, zero_mask, n_vert, _ = \
    load_data(dataset, 'group', physio=[physio], load_physio=True, regress_global=regress_global_sig) 

    physio_data = np.squeeze(np.stack(physio_sig, axis=1))
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
    # Write out results
    write_results(dataset, physio, pred_maps, time_lag_maps, pred_lag_vec, 
                  zero_mask, n_vert)

def write_results(dataset, term, pred_maps, time_lag_maps, pred_lag_vec, zero_mask, n_vert):
    analysis_str = f'{dataset}_physio_reg_group_{term}'
    # if time-lag maps specified, get lag of maximum/minimum cross-correlation of each voxel.
    if time_lag_maps:
        # Minimum (negative) cross-correlation
        neg_indx = np.argmin(pred_maps, axis=0)
        neg_lag_indx = np.array([pred_lag_vec[i] for i in neg_indx])
        write_nifti(neg_lag_indx[np.newaxis, :], f'{analysis_str}_min_time_lag.nii', zero_mask, n_vert)
        # Maximum (positive cross-correlation)
        pos_indx = np.argmax(pred_maps, axis=0)
        pos_lag_indx = np.array([pred_lag_vec[i] for i in pos_indx])
        write_nifti(pos_lag_indx[np.newaxis, :], f'{analysis_str}_max_time_lag.nii', zero_mask, n_vert)
    else:
        write_nifti(pred_maps, analysis_str, zero_mask, n_vert)

        





if __name__ == '__main__':
    """Run main analysis"""
    parser = argparse.ArgumentParser(description='Run Physio GLM w/ Time-lag Spline Regressors')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to run analysis on',
                        choices=['chang', 'nki', 'yale', 'hcp', 'hcp_fix'], 
                        required=True,
                        type=str)
    parser.add_argument('-p', '--physio',
                        help='select physio - can only provide one',
                        required=True,
                        choices=['hr', 'rv', 'alpha', 'delta', 'infraslow', 'ppg_low'],
                        type=str)
    parser.add_argument('-pl', '--p_nlags',
                        help='Number of lags (TRs) of physio signal in the positive (backward) direction',
                        required=False,
                        default=0,
                        type=int) 
    parser.add_argument('-nl', '--n_nlags',
                        help='Number of lags (TRs) of physio signal in the negative (forward) direction',
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
    parser.add_argument('-t', '--time_lag_maps',
                        help='Whether to create time lag maps between voxel time series and regressor',
                        default=0,
                        required=False,
                        type=int)

    args_dict = vars(parser.parse_args())
    run_main(args_dict['dataset'], args_dict['physio'], 
             args_dict['p_nlags'], args_dict['n_nlags'], 
             args_dict['nknots'], args_dict['regress_global_sig'], 
             args_dict['n_eval'], args_dict['time_lag_maps'])
