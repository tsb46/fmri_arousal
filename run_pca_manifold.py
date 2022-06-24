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


def construct_spline_basis(pca_ts, nknots):
    # Define model formula
    # Create lag splines
    spline_basis = dmatrix("te(cr(x, df=nknots), cr(y, df=nknots)) - 1", 
                          {"x": pca_ts[:,0], 'y': pca_ts[:,1]}, 
                          return_type='dataframe')
    return spline_basis


def evaluate_model(p1_eval, p2_eval, model, spline_basis):
    pred_maps = []
    for p1_e in p1_eval:
        p1_pred_maps = []
        for p2_e in p2_eval:
            # Create basis from model evaluation using previously defined design matrix (for model fit)
            pred_mat = dmatrix(spline_basis.design_info, {'x': p1_e, 'y': p2_e}, 
                                     return_type='dataframe')
            # Get predictions from model
            pred_map = model.predict(pred_mat)
            p1_pred_maps.append(pred_map)
        pred_maps.append(np.vstack(p1_pred_maps))
    return pred_maps


def run_main(dataset, pca_res_fp, p1, p2, nknots, regress_global_sig, n_eval):
    # Load data
    func_data, _, _, zero_mask, n_vert, _ = load_data(dataset, 'group', physio=None, load_physio=False, 
                                                      regress_global=regress_global_sig) 

    # Load pca results
    pca_res = pickle.load(open(pca_res_fp, 'rb'))
    pca_ts = pca_res['pc_scores'][:, [p1, p2]]

    # Construct Design matrix using patsy style formula
    print('construct spline matrix')
    spline_basis = construct_spline_basis(pca_ts, nknots)

    # Run regression
    print('run regression')
    lin_reg = linear_regression(spline_basis, func_data, return_model=True, 
                                intercept=False, norm=False)
    # Get predicted maps at all time lags at equally spaced percentiles
    print('get predicted maps at different levels of p1 and p2')
    p1_5, p1_95 = np.percentile(pca_ts[:,0], [5, 95])
    p2_5, p2_95 = np.percentile(pca_ts[:,1], [5, 95])
    p1_eval = np.linspace(p1_5, p1_95, n_eval)
    p2_eval = np.linspace(p2_5, p2_95, n_eval)
    pred_maps = evaluate_model(p1_eval, p2_eval, lin_reg, spline_basis)
    # Write out results
    write_results(dataset, pred_maps, [p1_eval, p2_eval], zero_mask, n_vert)


def write_results(dataset, pred_maps, pred_vec, zero_mask, n_vert):
    analysis_str = f'{dataset}_pc_manifold_group'
    # if time-lag maps specified, get lag of maximum/minimum cross-correlation of each voxel.
    for i, maps in enumerate(pred_maps):
        write_nifti(maps, f'{analysis_str}_{i}', zero_mask, n_vert)

        





if __name__ == '__main__':
    """Run main analysis"""
    parser = argparse.ArgumentParser(description='Run Physio GLM w/ Time-lag Spline Regressors')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to run analysis on',
                        choices=['chang', 'nki', 'yale', 'hcp', 'chang_bh'], 
                        required=True,
                        type=str)
    parser.add_argument('-r', '--pca_results',
                        help='file path to pca results pickle dictionary',
                        required=True,
                        type=str)
    parser.add_argument('-p1', '--pc_index1',
                        help='numeric index of principal component (e.g. 0 for first component)',
                        required=True,
                        choices=[0,1,2],
                        type=int)
    parser.add_argument('-p2', '--pc_index2',
                        help='numeric index of principal component',
                        required=True,
                        choices=[0,1,2],
                        type=int)
    parser.add_argument('-k', '--nknots',
                        help='Number of knots in spline basis',
                        default=5, 
                        required=False,
                        type=int)    
    parser.add_argument('-g', '--regress_global_sig',
                        help='Whether to regress out global signal from functional data',
                        default=0,
                        required=False,
                        type=int)
    parser.add_argument('-e', '--n_eval',
                        help='Number of equally spaced samples to evaluate model response',
                        default=5, 
                        required=False,
                        type=int) 

    args_dict = vars(parser.parse_args())
    run_main(args_dict['dataset'], args_dict['pca_results'], args_dict['pc_index1'],
             args_dict['pc_index2'], args_dict['nknots'], 
             args_dict['regress_global_sig'], args_dict['n_eval'])
