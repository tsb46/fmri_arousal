import argparse
import numpy as np
import pandas as pd
import pickle

from itertools import zip_longest
from utils.glm_utils import construct_design_matrix, \
lag_and_convolve_physio, linear_regression, create_interaction_maps, \
get_interaction_map, get_quadratic_map, create_quadratic_maps, \
parse_quadratic_string, mask_voxels
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


def run_main(dataset, model_formula, interaction_map, quadratic_map,
             time_lag, physio, convolve):
    
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
        beta_reg = linear_regression(design_mat, func_data, design_mat.columns, intercept=False, norm=False)
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
    # Parse quadratic string (if specified) and identify beta map associated with the quadratic term
    if quadratic_map is not None:
        avg_beta_quad, v1 = get_quadratic_map(quadratic_map, design_mat.columns, 
                                              subj_beta_maps)

        
    # Write out group map for each covariate in model
    for i, term in enumerate(design_mat.columns):
        # Modify quadratic term strings for output file path
        if '**' in term:
            var = parse_quadratic_string(term)
            term = f'{var.strip()}_2'
        # Create fixed effect map (averaged across subjects)
        avg_beta = np.mean([bmap[i] for bmap in subj_beta_maps], axis=0) 
        write_results(dataset, term, avg_beta[np.newaxis, :], 'group', None, None, zero_mask, n_vert)
        # If specified, and covariate is part of an interaction, create simple effect maps
        if (interaction_map is not None) and ((term == v1_i) | (term == v2_i)):
            interaction_beta_maps = create_interaction_maps(avg_beta, avg_beta_inter)
            write_results(dataset, f'{term}_interaction_map', 
                          interaction_beta_maps, 'group', None, None, zero_mask, n_vert)
        # If specified, and covariate is specified in a quadratic term, create simple effect maps
        if (quadratic_map is not None) and (term == v1):
            quadratic_beta_maps = create_quadratic_maps(avg_beta, avg_beta_quad)
            write_results(dataset, f'{term}_interaction_map', 
                          quadratic_beta_maps, 'group', None, None, zero_mask, n_vert)






if __name__ == '__main__':
    """Run main analysis"""
    parser = argparse.ArgumentParser(description='Regress functional data on '
                                     'physio time series')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to run analysis on',
                        choices=['chang', 'nki', 'yale', 'hcp', 'lemon'], 
                        required=True,
                        type=str)
    parser.add_argument('-f', '--model_formula',
                        help='model formula string using patsy-style formula. To specify interaction term between two terms, use '
                        'the character ":" between two terms w/ no spaces. Must include main effect terms w/ interaction term. '
                        ' To specify a quadratic term (no higher terms beyond quadratic supported), use the '
                        'following (for example var A): I(A**2)). Must include simple (linear) effect w/ quadratic term.',
                        required=True,
                        type=str)
    parser.add_argument('-i', '--interaction_map',
                        help='interaction effect string from model formula (option -f) that signifies the user would '
                        'like simple effects created. This string should match the interaction effect string used '
                        'in the patsy formula (option -f). This should only be used if an interaction effect was specified in'
                        'the model',
                        default=None,
                        required=False,
                        type=str)
    parser.add_argument('-q', '--quadratic_map',
                        help='quadratic term string from model formula (option -f) that signifies the user would '
                        'like simple effect maps created. This string should match the quadratic term string used '
                        'in the patsy formula (option -f). This should only be used if an quadratic was specified in'
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
             args_dict['quadratic_map'], args_dict['time_lag'], args_dict['physio'], 
             args_dict['convolve'])

