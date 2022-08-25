import argparse
import numpy as np
import fbpca
import pandas as pd
import pickle

from numpy.linalg import pinv
from patsy import dmatrix
from run_pca import pca, rotation
from run_physio_glm_dlnlm import construct_crossbasis
from scipy.stats import zscore
from scipy.signal import hilbert
from utils.load_write import load_data, write_nifti
from utils.glm_utils import linear_regression


def run_main(dataset, n_comps, physio, n_lags, lag_nknots, var_nknots, rotate, regress_global_sig):
    # Load data
    func_data, physio_sig, physio_labels, zero_mask, n_vert, params = \
    load_data(dataset, 'group', physio=[physio], load_physio=True, regress_global=regress_global_sig) 

    # Create dataframe of physio signals
    physio_sig = pd.DataFrame(np.squeeze(np.stack(physio_sig,axis=1)), 
                              columns=physio_labels)

    # Construct Design matrix using patsy style formula
    print('construct spline matrix')
    # Construct Design matrix using patsy style formula
    design_mat, basis_var, basis_lag = construct_crossbasis(physio_sig[physio], 
                                                            n_lags, var_nknots, lag_nknots)
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
                  rotate, zero_mask, n_vert, params)


def write_results(dataset, res_list, res_labels, rotate, zero_mask, n_vert, params):
    analysis_str = f'{dataset}_rda'
    if rotate is not None:
        analysis_str += f'_{rotate}'
    results_dict = {l: r for r,l in zip(res_list, res_labels)}
    write_nifti(results_dict['rda_comps']['loadings'], analysis_str, zero_mask, n_vert, params['mask'])
    write_nifti(results_dict['rda_comps_residuals']['loadings'], f'{analysis_str}_resid', zero_mask, n_vert, params['mask'])
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
                        help='select physio',
                        required=False,
                        default=None,
                        type=str)
    parser.add_argument('-l', '--nlags',
                        help='Number of lags',
                        default=15, 
                        required=False,
                        type=int)  
    parser.add_argument('-lk', '--lag_nknots',
                        help='Number of knots in spline basis for lag. '
                        'Knots are placed along the range of lag values at a log '
                        'scale (more resolution at earlier lags)',
                        default=3, 
                        required=False,
                        type=int)   
    parser.add_argument('-vk', '--var_nknots',
                        help='Number of knots in spline basis for physio var. '
                        'Knots are placed at equally spaced quantiles based on N knots',
                        default=5, 
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
             args_dict['nlags'], args_dict['lag_nknots'], args_dict['var_nknots'],
             args_dict['rotate'], args_dict['regress_global_sig'])
