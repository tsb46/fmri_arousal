import argparse
import numpy as np
import pandas as pd
import pickle

from utils.glm_utils import construct_design_matrix, convolve_hrf, \
get_hrf, lag_and_convolve_physio, linear_regression, onsets_to_block, \
create_interaction_maps, get_interaction_map
from scipy.stats import zscore
from utils.load_write import load_data, write_nifti, load_nki_event_file


def average_scans(func_data):
    scan_avg = np.mean(func_data, axis=0)
    return np.squeeze(scan_avg)


def construct_task_blocks(scan_len, event_df, tr):
    event_df = group_task_blocks(event_df)
    group_func = {'onset': 'first', 'duration': 'sum'}
    # Construct 'deep breath' blocks
    breath_df = event_df.loc[event_df.breath]
    breath_df_grp = breath_df.groupby(['breath_indx']).agg(group_func)
    breath_blocks = onsets_to_block(breath_df_grp, scan_len, tr)
    # Construct 'breath hold' blocks
    hold_df = event_df.loc[event_df.hold]
    hold_df_grp = hold_df.groupby(['hold_indx']).agg(group_func)
    hold_blocks = onsets_to_block(hold_df_grp, scan_len, tr)
    return breath_blocks, hold_blocks


def construct_task_regressors(event_df, scan_len, n_scans, c_hrf, r_hrf, tr):
    breath_blocks, hold_blocks = construct_task_blocks(scan_len, event_df, tr)
    breath_blocks_c = convolve_hrf(r_hrf, breath_blocks)
    hold_blocks_c = convolve_hrf(c_hrf, hold_blocks)
    df_blocks = pd.DataFrame({'breath': breath_blocks_c, 
                                'hold': hold_blocks_c})
    return df_blocks, breath_blocks, hold_blocks


def group_task_blocks(event_df):
    event_df['hold'] = event_df.trial_type.str.startswith('H')
    event_df['hold_indx'] = (event_df.hold != event_df.hold.shift(1)).cumsum()
    event_df['breath'] = event_df.trial_type.isin(['Deep','In','Out'])
    event_df['breath_indx'] = (event_df.breath != event_df.breath.shift(1)).cumsum()
    return event_df


def write_results(dataset, term, beta_map, level, subj_n, scan, zero_mask, n_vert):
    if level == 'group':
        analysis_str = f'{dataset}_taskbh_group_{term}'
    elif level == 'subject':
        analysis_str = f'{dataset}_taskbh_s{subj_n}'
        if scan is not None:
            analysis_str += f'_{scan}_{term}'

    write_nifti(beta_map, analysis_str, zero_mask, n_vert)


def run_main(dataset, model_formula, time_lag, interaction_map, 
             level, subj_n, scan_n, convolve, scan_average):
    # Load data
    func_data, physio_sig, physio_labels, zero_mask, n_vert, params = load_data(dataset, level, physio=None, load_physio=True, 
                                                                                subj_n=subj_n, scan_n=scan_n, group_method='list') 
    if level == 'subject':
        func_data = [func_data]
        for i, p in enumerate(physio_sig):
            physio_sig[i] = [p]

                
    # Create lagged physio variables (if non-zero lag)
    physio_sig_proc = lag_and_convolve_physio(physio_sig, physio_labels, 
                                              len(func_data), time_lag, 
                                              convolve, params['tr'])

    # We are assuming that all timing is the same, it looks that this is the case from 
    # manual inspection of several files. NEED TO DOUBLE CHECK!!!!
    # compute across subject task average
    if scan_average:
        scan_avg = average_scans(func_data)
        write_results(dataset, 'scan_avg', scan_avg, level, subj_n, scan_n, zero_mask, n_vert)

    # load NKI event file (assuming timing is the same across subjects, see note above)
    df_events = load_nki_event_file()
    # Get canonical hrf
    c_hrf = get_hrf(30, params['tr'], 'canonical')
    r_hrf = get_hrf(30, params['tr'], 'rvt')
    # Create breathhold task regressors
    df_model, b, h = construct_task_regressors(df_events, func_data[0].shape[0], 
                                         params['nscans'], c_hrf, r_hrf, 
                                         params['tr'])

    # run linear regression for each subj
    subj_beta_maps = []
    for func, phys in zip(func_data, physio_sig_proc):
        df_physio = pd.DataFrame(phys, columns=physio_labels)
        df_full = pd.concat([df_model, df_physio], axis=1)
        design_mat = construct_design_matrix(model_formula, df_full)
        betas = linear_regression(design_mat, func, design_mat.columns)
        subj_beta_maps.append(betas)

    # Parse interaction string (if specified) and identify beta map associated with the interaction term
    if interaction_map is not None:
        avg_beta_inter, v1_i, v2_i = get_interaction_map(interaction_map, design_mat.columns, 
                                                         subj_beta_maps)

    for i, term in enumerate(design_mat.columns):
        avg_beta = np.mean([bmap[i] for bmap in subj_beta_maps], axis=0) 
        write_results(dataset, term, avg_beta[np.newaxis, :], level, subj_n, scan_n, zero_mask, n_vert)
        # If specified, and covariate is part of an interaction, create an interaction map
        if (interaction_map is not None) and ((term == v1_i) | (term == v2_i)):
            interaction_beta_maps = create_interaction_maps(avg_beta, avg_beta_inter)
            write_results(dataset, f'{term}_interaction_map', 
                          interaction_beta_maps, level, subj_n, scan_n, zero_mask, n_vert)


if __name__ == '__main__':
    """Run main analysis"""
    parser = argparse.ArgumentParser(description='Regress functional data on '
                                     'physio time series')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to run analysis on - only option is NKI',
                        choices=['nki'], 
                        default='nki',
                        required=False,
                        type=str)
    parser.add_argument('-f', '--model_formula',
                        help='model formula string using patsy-style formula - task blocks '
                        'regressors are "breath" and "hold"',
                        required=True,
                        type=str)
    parser.add_argument('-lag', '--time_lag', 
                        help='choice of lag (positive - shift to the right, negative - shift to the left) for physio time series',
                        default=0,
                        type=int)
    parser.add_argument('-i', '--interaction_map',
                        help='interaction effect string from model formula (option -f) that signifies the user would '
                        'like an interaction map created. This string should match the interaction effect string used '
                        'in the patsy formula (option -f) This should only be used if an interaction effect was specified'
                        'the model',
                        default=None,
                        required=False,
                        type=str)
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
    parser.add_argument('-c', '--convolve',
                        help='whether to convolve physio time series with canonical '
                        'hemodynamic response function',
                        default=0,
                        choices=[0,1],
                        type=int)
    parser.add_argument('-a', '--scan_average',
                        help='whether to compute a scan average by averaging time series across subjects',
                        default=0,
                        choices=[0,1],
                        type=int)


    args_dict = vars(parser.parse_args())
    run_main(args_dict['dataset'], args_dict['model_formula'], args_dict['time_lag'],
             args_dict['interaction_map'], args_dict['level'], args_dict['subject_n'], args_dict['scan_n'],
             args_dict['convolve'], args_dict['scan_average'])

