import argparse
import numpy as np
import pandas as pd
import pickle

from patsy import dmatrix
from scipy.stats import zscore
from utils.glm_utils import linear_regression, mask_voxels
from utils.load_utils import load_subject_list
from utils.load_write import load_data, write_nifti


def construct_crossbasis(pvar, p_nlags, n_nlags, var_df, lag_df):
    # Create lag sequence array (include lag of 0!)
    seq_lag = np.arange(-n_nlags, p_nlags+1)
    # Create Cubic B-spline basis for predictor and lag
    basis_var = dmatrix("cr(x, df=var_df) - 1", {"x": pvar}, return_type='dataframe')
    basis_lag = dmatrix("cr(x, df=lag_df) - 1", {"x": seq_lag}, return_type='dataframe')
    # Intialize crossbasis matrix
    crossbasis = np.zeros((len(pvar), basis_var.shape[1]*basis_lag.shape[1]))
    # Loop through predictor and lag bases and multiply column pairs
    indx = 0
    for v in np.arange(basis_var.shape[1]):
        lag_mat = pd.concat([basis_var.iloc[:,v].shift(i) for i in seq_lag], axis=1)
        for l in np.arange(basis_lag.shape[1]):
            crossbasis[:, indx] = np.dot(lag_mat.values, basis_lag.iloc[:,l].values)
            indx+=1
    return crossbasis, basis_var, basis_lag
    

def evaluate_model(pvar, lin_reg, basis_var, basis_lag, lag_eval, p_nlags, n_nlags, n_voxels):
    # Select z-score units to evaluate physio variables
    var_pred_array = [0.5, 0.75, 1, 1.5, 2, 2.5]
    # Create lag sequence array (include lag of 0!)
    seq_lag = np.linspace(-n_nlags, p_nlags, lag_eval)
    # Create repeated values of pred and lag array for all possible pairwise combos
    varvec = np.tile(var_pred_array, len(seq_lag))
    lagvec = np.repeat(seq_lag,len(var_pred_array))
    # Define length of pred and lag array
    n_var = len(var_pred_array)
    n_lag = len(seq_lag)
    # Create basis from model evaluation using previously defined design matrix (for model fit)
    basis_var_pred = dmatrix(basis_var.design_info, {'x': varvec}, return_type='dataframe')
    basis_lag_pred = dmatrix(basis_lag.design_info, {'x': lagvec}, return_type='dataframe')

    # We must center our predictions around a reference value (set at 50th percentile/median)
    cen = pvar.quantile(0.5) 
    # Rescale 
    basis_cen = dmatrix(basis_var.design_info, {'x': cen}, return_type='dataframe')
    basis_var_pred = basis_var_pred.subtract(basis_cen.values, axis=1)

    v_len = basis_var_pred.shape[1]
    l_len = basis_lag_pred.shape[1]
    # Row-wise kronecker product between predicted bases to generate prediction matrix
    xpred_list = [basis_var_pred.iloc[:,v]*basis_lag_pred.iloc[:,l] 
                  for v in range(v_len) for l in range(l_len)]
    xpred = pd.concat(xpred_list, axis=1)
    # Get predictions from model
    pred_mat = lin_reg.predict(xpred)
    # v_pred_mat = pd.DataFrame(pred_mat.reshape(n_var, n_lag, order='F'), index=var_pred_array, columns=seq_lag)
    pred_all = reshape_output(pred_mat, n_var, n_lag)
    return pred_all, seq_lag, var_pred_array, 


def reshape_output(pred_mat, n_var, n_lag):
    # Loop through voxels and reshape predictions into physio val by lag matrix
    pred_list = []
    for i in range(pred_mat.shape[1]):
        v_pred_mat = pred_mat[:,i].reshape(n_var, n_lag, order='F')
        pred_list.append(v_pred_mat)
    return pred_list

def save_func_pred(design_mat, lin_reg, dataset):
    func_data_pred = lin_reg.predict(design_mat)
    analysis_str = f'{dataset}_dlnm_group_func_pred'
    pickle.dump(func_data_pred, open(analysis_str, 'wb'))


def write_results(dataset, term, beta_map, level, subj_n, scan, zero_mask, n_vert, params):
    if level == 'group':
        analysis_str = f'{dataset}_dlnm_group_{term}'
    write_nifti(beta_map, analysis_str, zero_mask, n_vert, params['mask'])


def run_main(dataset, physio_var, p_nlags, n_nlags, var_df, lag_df, 
             lag_eval, save_pred):
    func_data, physio_sig, physio_labels, zero_mask, n_vert, params = load_data(dataset, physio=[physio_var],
                                                                                load_physio=True) 
    # Create dataframe of physio signals
    physio_sig = pd.DataFrame(np.squeeze(np.stack(physio_sig,axis=1)), 
                              columns=physio_labels)

    # Construct Design matrix using patsy style formula
    crossbasis, basis_var, basis_lag = construct_crossbasis(physio_sig[physio_var], p_nlags, n_nlags, var_df, lag_df)
    # Lag introduces null values - trim beginning of predictor matrix
    na_indx = ~(np.isnan(crossbasis).any(axis=1))
    func_data = func_data[na_indx, :]

    lin_reg = linear_regression(crossbasis[na_indx,:], func_data, return_model=True, 
                                intercept=False, norm=False)
    if save_pred:
        del func_data # free up memory
        save_pred(cross_basis[na_indx,:], lin_reg, dataset)

    pred_maps, lag_eval, eval_points = evaluate_model(physio_sig[physio_var], lin_reg, basis_var, basis_lag, 
                                                      lag_eval, p_nlags, n_nlags, func_data.shape[1])
    pred_maps = np.stack(pred_maps,axis=2)

    pickle.dump([lin_reg, lag_eval, eval_points], open(f'{dataset}_dlnm_group_results.pkl', 'wb'))
    for i in range(pred_maps.shape[0]):
        write_results(dataset, f'eval_{i}', pred_maps[i,:,:], 'group', None, None, zero_mask, n_vert, params)


if __name__ == '__main__':
    """Run main analysis"""
    parser = argparse.ArgumentParser(description='Regress functional data on '
                                     'physio time courses using dynamic lag non linear model')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to run analysis on',
                        choices=['chang', 'chang_bh', 'nki', 'yale',
                                  'hcp', 'spreng'], 
                        required=True,
                        type=str)
    parser.add_argument('-p', '--physio_var',
                        help='<Required> physio predictor var',
                        required=True,
                        type=str)
    parser.add_argument('-pl', '--p_nlags',
                        help='Number of lags (TRs) of physio signal in the positive (forward) direction',
                        required=False,
                        default=15,
                        type=int) 
    parser.add_argument('-nl', '--n_nlags',
                        help='Number of lags (TRs) of physio signal in the negative (backward) direction',
                        required=False,
                        default=0,
                        type=int)   
    parser.add_argument('-vd', '--var_df',
                        help='Number of degrees of freedom in spline basis for physio var. '
                        'Knots are placed at equally spaced quantiles based on N knots',
                        default=6, 
                        required=False,
                        type=int)  
    parser.add_argument('-ld', '--lag_df',
                        help='Number of degrees of freedom in spline basis for lag. '
                        'Knots are placed along equally spaced quantiles across lag range',
                        default=5, 
                        required=False,
                        type=int)  
    parser.add_argument('-le', '--lag_evaluate',
                        help='for model predictions, what # of equally spaced values between '
                        ' min and max lag to evaluate on the fitted model.',
                        required=False,
                        default=40,
                        type=int)
    parser.add_argument('-s', '--save_func_pred',
                        help='whether to save the predicted functional BOLD data (will be very large)',
                        required=False,
                        default=0,
                        type=int)


    args_dict = vars(parser.parse_args())
    run_main(args_dict['dataset'], args_dict['physio_var'], args_dict['p_nlags'], 
             args_dict['n_nlags'], args_dict['var_df'], args_dict['lag_df'], 
             args_dict['lag_evaluate'], args_dict['save_func_pred'])
