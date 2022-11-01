import argparse
import nibabel as nb
import numpy as np
import fbpca
import pandas as pd
import pickle

from numpy.linalg import pinv
from patsy import dmatrix
from run_physio_glm_dlnlm import construct_crossbasis
from scipy.stats import zscore
from scipy.signal import hilbert
from utils.load_write import convert_2d, load_data, write_nifti
from utils.glm_utils import linear_regression


def load_parcellation(parcellation_fp, mask, data_indx):
    nifti = nb.load(parcellation_fp, keep_file_open = True)
    # Load mask
    mask = nb.load(mask).get_fdata() > 0
    nifti_data = nifti.get_fdata()
    nifti.uncache()
    nifti_data = convert_2d(mask, nifti_data)
    return nifti_data[data_indx]


def run_main(dataset, parcellation, physio, n_lags, lag_nknots, var_nknots, regress_global_sig):
    # Load data
    func_data, physio_sig, physio_labels, zero_mask, n_vert, params = \
    load_data(dataset, 'group', physio=[physio], load_physio=True, regress_global=regress_global_sig) 

    # Create dataframe of physio signals
    physio_sig = pd.DataFrame(np.squeeze(np.stack(physio_sig,axis=1)), 
                              columns=physio_labels)

    parcel_data = load_parcellation(parcellation, params['mask'], zero_mask)
    parcel_indx = np.unique(parcel_data)

    # Construct FC matrix for original data
    parcel_ts_orig = np.array([func_data[:, parcel_data == i].mean(axis=1) for i in parcel_indx])
    fc_orig = np.corrcoef(parcel_ts_orig)


    # Construct Design matrix using patsy style formula
    design_mat, basis_var, basis_lag = construct_crossbasis(physio_sig[physio], 
                                                            n_lags, var_nknots, lag_nknots)

    # Lag introduces null values - trim beginning of predictor matrix
    na_indx = ~(np.isnan(design_mat).any(axis=1))
    func_data = func_data[na_indx]
    design_mat = design_mat[na_indx, :]
    lin_reg = linear_regression(design_mat, func_data, return_model=True, 
                                intercept=False, norm=False)

    # Free up memory
    del func_data, parcel_ts_orig
    
    # Predict individual voxel ts from design mat
    func_data_pred = lin_reg.predict(design_mat)

    # Construct FC matrix for predicted data
    parcel_ts_pred = np.array([func_data_pred[:, parcel_data == i].mean(axis=1) for i in parcel_indx])
    fc_pred = np.corrcoef(parcel_ts_pred)

    write_results(dataset, physio, [fc_orig, fc_pred])


def write_results(dataset, physio, results):
    analysis_str = f'{dataset}_fc_compare_{physio}'
    pickle.dump(results, open(f'{analysis_str}_results.pkl', 'wb'))


if __name__ == '__main__':
    """Run main analysis"""
    parser = argparse.ArgumentParser(description='Run redundancy analysis')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to run analysis on',
                        choices=['chang', 'nki', 'yale', 'hcp', 'hcp_fix'], 
                        required=True,
                        type=str)
    parser.add_argument('-pa', '--parcellation',
                       help='filepath to parcellation nifti for extracting time courses and '
                       'constructing FC matrix',
                       required=True,
                       type=str)
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
    parser.add_argument('-g', '--regress_global_sig',
                        help='Whether to regress out global signal from functional data',
                        default=0,
                        required=False,
                        type=int)


    args_dict = vars(parser.parse_args())
    run_main(args_dict['dataset'], args_dict['parcellation'], 
             args_dict['physio'], args_dict['nlags'], args_dict['lag_nknots'], 
             args_dict['var_nknots'],args_dict['regress_global_sig'])
