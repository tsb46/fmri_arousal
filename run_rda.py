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


def construct_lag_splines(physio_vars, nlags, nknots):
    # Define model formula
    # Create lag sequence array (include lag of 0!)
    seq_lag = np.arange(nlags+1)
    # Create lag splines
    lag_splines = dmatrix("cr(x, df=nknots) - 1", 
                          {"x": seq_lag}, return_type='dataframe')
    # Create design matrix 
    design_mat = []
    for var_label in physio_vars.columns:
        var = physio_vars[var_label].copy()
        basis_lag = np.zeros((len(physio_vars), lag_splines.shape[1]))
        # Loop through lag bases and multiply column pairs
        lag_mat = pd.concat([var.shift(l) for l in seq_lag], axis=1)
        for l in np.arange(lag_splines.shape[1]):
            basis_lag[:, l] = np.dot(lag_mat.values, lag_splines.iloc[:,l].values)
        design_mat.append(basis_lag)
    design_mat = np.hstack(design_mat)
    return design_mat


def run_main(dataset, n_comps, physio, n_lags, nknots, rotate, regress_global_sig):
    # Load data
    func_data, physio_sig, physio_labels, zero_mask, n_vert, _ = \
    load_data(dataset, 'group', physio=physio, load_physio=True, regress_global=regress_global_sig) 

    physio_data = np.squeeze(np.stack(physio_sig, axis=1))
    physio_df = pd.DataFrame(physio_data, columns=physio_labels)
    # Construct Design matrix using patsy style formula
    print('construct spline matrix')
    design_mat = construct_lag_splines(physio_df, n_lags, nknots)
    # Lag introduces null values - trim beginning of predictor matrix
    na_indx = ~(np.isnan(design_mat).any(axis=1))
    func_data = func_data[na_indx]
    design_mat = design_mat[na_indx, :]
    print('run pca')
    # Run ordinary PCA 
    pca_res = pca(func_data, n_comps)
    # If rotation is chosen, rotate pca loadings
    if rotate is not None:
        pca_res = rotation(pca_res, func_data, rotate)
    print('run regression')
    lin_reg = linear_regression(design_mat, func_data, return_model=True, 
                                intercept=False, norm=False)
    print('regression projection')
    # Predict individual voxel ts from design mat
    func_data_pred = lin_reg.predict(design_mat)
    # Run PCA on predicted values
    print('run pca on predicted values')
    rda = pca(func_data_pred, n_comps)
    # If rotation is chosen, rotate rda loadings
    if rotate is not None:
        rda = rotation(rda, func_data_pred, rotate)
        rda_proj = rda['pc_scores']
        rda_proj_y = func_data @ pinv(rda['loadings'].T).T
    else:
        # Project predicted and observed functional data on RDA axes
        rda_proj = np.dot(func_data_pred, rda['Va'].T)
        rda_proj_y = np.dot(func_data, rda['Va'].T)

    # Run PCA on residuals (apply in-place to original func data to not overload memory)
    func_data -= func_data_pred
    # Free up memory
    del func_data_pred
    # Run pca on residuals
    print('run pca on residuals')
    rda_resid = pca(func_data, n_comps)
    # If rotation is chosen, rotate residual rda loadings
    if rotate is not None:
        rda_resid = rotation(rda_resid, func_data, rotate)
        rda_proj_resid = rda_resid['pc_scores']
    else:
        # Project residual functional data on RDA axes
        rda_proj_resid = np.dot(func_data, rda_resid['Va'].T)
    res_list = [rda, rda_proj, rda_proj_y, rda_resid, rda_proj_resid, pca_res]
    res_labels = ['rda_comps', 'rda_proj', 'rda_proj_y', 'rda_comps_residuals', 
                  'rda_proj_residuals','pca_comps']
    write_results(dataset, res_list, res_labels,
                  rotate, zero_mask, n_vert)


def write_results(dataset, res_list, res_labels, rotate, zero_mask, n_vert):
    analysis_str = f'{dataset}_rda'
    if rotate is not None:
        analysis_str += f'_{rotate}'
    results_dict = {l: r for r,l in zip(res_list, res_labels)}
    write_nifti(results_dict['rda_comps']['loadings'], analysis_str, zero_mask, n_vert)
    write_nifti(results_dict['rda_comps_residuals']['loadings'], f'{analysis_str}_resid', zero_mask, n_vert)
    pickle.dump(results_dict, open(f'{analysis_str}_results.pkl', 'wb'))



if __name__ == '__main__':
    """Run main analysis"""
    parser = argparse.ArgumentParser(description='Run redundancy analysis')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to run analysis on',
                        choices=['chang', 'nki', 'yale', 'hcp', 'hcp_fix'], 
                        required=True,
                        type=str)
    parser.add_argument('-n', '--n_comps',
                        help='<Required> Number of components from RDA',
                        required=True,
                        type=int)
    parser.add_argument('-p', '--physio',
                        help='select physio - can provide multiple (separated by space)',
                        required=False,
                        default=None,
                        action='append',
                        type=str)
    parser.add_argument('-l', '--nlags',
                        help='Number of lags',
                        default=10, 
                        required=False,
                        type=int)  
    parser.add_argument('-k', '--nknots',
                        help='Number of knots in spline basis',
                        default=3, 
                        required=False,
                        type=int)    
    parser.add_argument('-r', '--rotate',
                        help='Whether to rotate pca weights',
                        default=None,
                        required=False,
                        choices=['varimax', 'promax'],
                        type=str)
    parser.add_argument('-g', '--regress_global_sig',
                        help='Whether to regress out global signal from functional data',
                        default=0,
                        required=False,
                        type=int)

    args_dict = vars(parser.parse_args())
    run_main(args_dict['dataset'], args_dict['n_comps'], args_dict['physio'], 
             args_dict['nlags'], args_dict['nknots'], args_dict['rotate'], 
             args_dict['regress_global_sig'])
