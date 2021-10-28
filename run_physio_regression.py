import argparse
import numpy as np
import pandas as pd
import pickle

from itertools import zip_longest
from utils.glm_utils import construct_design_matrix, \
lag_and_convolve_physio, linear_regression, create_interaction_maps, \
get_interaction_map, mask_voxels
from utils.load_utils import load_subject_list
from patsy import dmatrix
from scipy.stats import zscore
from utils.load_write import load_data, write_nifti


def write_results(dataset, term, beta_map, level, subj_n, scan, zero_mask, n_vert):
    if level == 'group':
        analysis_str = f'{dataset}_physio_reg_group_{term}'
    elif level == 'subject':
        analysis_str = f'{dataset}_physio_reg_s{subj_n}'
        if scan is not None:
            analysis_str += f'_{scan}_{term}'
    write_nifti(beta_map, analysis_str, zero_mask, n_vert)


def run_main(dataset, model_formula, interaction_map, time_lag, 
             physio, convolve):
    
    subject_df = load_subject_list(dataset)
    if dataset == 'chang':
        subj_list = subject_df.subject
        scan_list = subject_df.scan
    else:
        subj_list = subject_df.subject
        scan_list = [None]

    subj_beta_maps = []
    # Load through subjects and fit first-level subject maps
    for subj, scan in zip_longest(subj_list, scan_list):
        print(subj)
        func_data, physio_sig, physio_labels, zero_mask, n_vert, params = load_data(dataset, 'subject', physio=physio,
                                                                                    load_physio=True, subj_n=subj, 
                                                                                    scan_n=scan, verbose=False, 
                                                                                    group_method='list', 
                                                                                    filter_nan_voxels=False) 

        # Create lagged physio variables (if non-zero lag) and convolve (if specified)
        physio_sig_proc = lag_and_convolve_physio(physio_sig, physio_labels, 1, time_lag, convolve, params['tr'])
        design_mat = construct_design_matrix(model_formula, pd.DataFrame(physio_sig_proc[0], columns=physio_labels))
        func_data, mask, beta_empty = mask_voxels(func_data[0])
        beta_reg = linear_regression(design_mat, func_data, design_mat.columns)
        beta_maps = []
        for b in beta_reg:
            b_all = beta_empty.copy()
            b_all[mask] = b
            beta_maps.append(b_all)
        subj_beta_maps.append(beta_maps)

    # Parse interaction string (if specified) and identify beta map associated with the interaction term
    if interaction_map is not None:
        avg_beta_inter, v1_i, v2_i = get_interaction_map(interaction_map, design_mat.columns, 
                                                         subj_beta_maps)
        
        
    # Write out individual maps for each covariate in model
    for i, term in enumerate(design_mat.columns):
        avg_beta = np.mean([bmap[i] for bmap in subj_beta_maps], axis=0) 
        write_results(dataset, term, avg_beta[np.newaxis, :], 'group', None, None, zero_mask, n_vert)
        # If specified, and covariate is part of an interaction, create an interaction map
        if (interaction_map is not None) and ((term == v1_i) | (term == v2_i)):
            interaction_beta_maps = create_interaction_maps(avg_beta, avg_beta_inter)
            write_results(dataset, f'{term}_interaction_map', 
                          interaction_beta_maps, 'group', None, None, zero_mask, n_vert)





if __name__ == '__main__':
    """Run main analysis"""
    parser = argparse.ArgumentParser(description='Regress functional data on '
                                     'physio time series')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to run analysis on',
                        choices=['chang', 'nki', 'yale', 'hcp'], 
                        required=True,
                        type=str)
    parser.add_argument('-f', '--model_formula',
                        help='model formula string using patsy-style formula',
                        required=True,
                        type=str)
    parser.add_argument('-i', '--interaction_map',
                        help='interaction effect string from model formula (option -f) that signifies the user would '
                        'like an interaction map created. This string should match the interaction effect string used '
                        'in the patsy formula (option -f) This should only be used if an interaction effect was specified'
                        'the model',
                        default=None,
                        required=False,
                        type=str)
    parser.add_argument('-lag', '--time_lag', 
                        help='choice of lag (positive - shift to the right, negative - shift to the left) for physio time series',
                        default=0,
                        type=int)
    parser.add_argument('-p', '--physio',
                        help='select physio - can provide multiple (separated by space)',
                        required=False,
                        default=None,
                        action='append',
                        type=str)
    parser.add_argument('-c', '--convolve',
                        help='whether to convolve physio time series with canonical '
                        'hemodynamic response function',
                        default=0,
                        choices=[0,1],
                        type=int)

    args_dict = vars(parser.parse_args())
    run_main(args_dict['dataset'], args_dict['model_formula'], args_dict['interaction_map'],
              args_dict['time_lag'], args_dict['physio'], args_dict['convolve'])

