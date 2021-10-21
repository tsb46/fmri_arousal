import argparse
import numpy as np
import pandas as pd
import pickle

from utils.glm_utils import construct_design_matrix, convolve_hrf, \
double_gamma_hrf, lag, linear_regression, create_interaction_maps, \
parse_interaction_string
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
             level, subj_n, scan_n, physio, convolve):
    # Load data
    func_data, physio_sig, physio_labels, zero_mask, n_vert, params = load_data(dataset, level, physio=physio, load_physio=True, 
                                                                                subj_n=subj_n, scan_n=scan_n, group_method='list') 

    # If specified, convolve physio with hemodynamic function
    if convolve:
        hrf = double_gamma_hrf(30, params['tr'])
                
    # Create lagged (if non-zero lag specified) and hrf convolved (if specified) signals 
    physio_sig_proc = []
    for subj_n in range(len(func_data)):
        subj_phys = []
        for p in physio_sig:
            if convolve:
                subj_sig = convolve_hrf(hrf, p[subj_n])
            else:
                subj_sig = p[subj_n]
            physio_sig_lag = lag(subj_sig, time_lag)
            subj_phys.append(physio_sig_lag)
        physio_sig_proc.append(np.stack(subj_phys, axis=1))

    # run linear regression for each subj
    subj_beta_maps = []
    for func, phys in zip(func_data, physio_sig_proc):
        design_mat = construct_design_matrix(model_formula, pd.DataFrame(phys, columns=physio_labels))
        betas = linear_regression(design_mat, func, design_mat.columns)
        subj_beta_maps.append(betas)

    # Parse interaction string (if specified) and identify beta map associated with the interaction term
    if interaction_map is not None:
        if interaction_map not in design_mat.columns:
            raise Exception('Interaction string specified for option -i does not match any string in model formula')
        else:
            i_index = design_mat.columns.tolist().index(interaction_map)
        v1v2 = interaction_map
        v1, v2 = parse_interaction_string(v1v2)
        # Remember, the beta maps are ordered according to the order of the columns in the design_mat dataframe
        avg_beta_inter = np.mean([bmap[i_index] for bmap in subj_beta_maps], axis=0) 

    # Write out individual maps for each covariate in model
    for i, term in enumerate(design_mat.columns):
        avg_beta = np.mean([bmap[i] for bmap in subj_beta_maps], axis=0) 
        write_results(dataset, term, avg_beta[np.newaxis, :], level, subj_n, scan_n, zero_mask, n_vert)
        # If specified, and covariate is part of an interaction, create an interaction map
        if (interaction_map is not None) and ((term == v1) | (term == v2)):
            interaction_beta_maps = create_interaction_maps(avg_beta, avg_beta_inter)
            write_results(dataset, f'{term}_interaction_map', 
                          interaction_beta_maps, level, subj_n, scan_n, zero_mask, n_vert)





if __name__ == '__main__':
    """Run main analysis"""
    parser = argparse.ArgumentParser(description='Regress functional data on '
                                     'physio time series')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to run analysis on',
                        choices=['chang', 'choe', 'gu', 'nki', 'yale', 'hcp'], 
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
    parser.add_argument('-l', '--level',
                        help='subject or group level analysis',
                        default='group',
                        choices=['subject', 'group'],
                        type=str)
    parser.add_argument('-s', '--subject_n',
                        help='subject number for subject level analysis',
                        default=None,
                        type=int)
    parser.add_argument('-scan', '--scan_n',
                        help='scan number for subject level analysis (if multiple runs from same subject',
                        default=None,
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
              args_dict['time_lag'], args_dict['level'], args_dict['subject_n'], args_dict['scan_n'],
              args_dict['physio'], args_dict['convolve'])

