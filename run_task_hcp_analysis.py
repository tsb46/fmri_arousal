import argparse
import numpy as np
import pandas as pd
import pickle


from patsy import dmatrix
from scipy.stats import zscore
from run_physio_glm import construct_lag_splines
from utils.glm_utils import convolve_hrf, get_hrf, lag as lag_shift, \
linear_regression, onsets_to_block, construct_lag_splines
from utils.load_write import load_data, write_nifti, load_hcp_task_event_file
from utils.load_utils import fp_hcp_rel, fp_hcp_wm, subject_list_hcp_rel, subject_list_hcp_wm


def conditional_average_block(func_data, physio_data, blocks, block_lag):
    physio_bins = [-10,-1.5,0,1.5,10]
    func_blocks = {n: [] for n in range(1,len(physio_bins))}
    func_blocks_indx = {n: [] for n in range(1,len(physio_bins)+1)}
    for i, (func_subj, physio_subj, block_subj) in enumerate(zip(func_data, physio_data, blocks)):
        blocks_pval = [physio_subj[block[0]] for block in block_subj] 
        blocks_bin_indx = np.digitize(blocks_pval, physio_bins)
        for bin_i in np.unique(blocks_bin_indx):
            blocks_i = np.where(blocks_bin_indx == bin_i)[0]
            func_blocks_indx[bin_i].append((i, blocks_i, [blocks_pval[b] for b in blocks_i]))
            for b in blocks_i: 
                func_tmp = func_subj[block_subj[b][0]:(block_subj[b][-1] + block_lag), :] 
                func_blocks[bin_i].append(func_tmp)
    
    for n in func_blocks:
        func_blocks[n] = np.mean(func_blocks[n], axis=0)

    return func_blocks, func_blocks_indx


def construct_interaction_terms(physio, task):
    n_task = task.shape[1]
    n_physio = physio.shape[1]
    design_mat = np.vstack([physio.iloc[:,p].values*task[:,t] for p in range(n_physio) for t in range(n_task)]).T
    design_mat = np.hstack([task, physio, design_mat])
    return design_mat


def construct_hcp_task_regressors(task, block_len, subject_list, event_loader, tr, scan_len, output_type):
    if output_type == 'block':
        lag_spline = None

    task_design = []
    for subj, scan in zip(subject_list.subject, subject_list.lr):
        event_dir = event_loader('events', subj, scan)
        event_df = load_hcp_task_event_file(event_dir, task)
        event_blocks = onsets_to_block(event_df, scan_len, tr)
        if output_type == 'glm':
            block_ts = np.zeros(scan_len)
            for block_i in event_blocks:
                block_ts[block_i[0]] = 1
            event_blocks = block_ts
        task_design.append(event_blocks)

    if output_type == 'glm':
        task_design = np.hstack(task_design)
        task_design, lag_spline = construct_lag_splines(pd.Series(task_design), block_len, 0, 4)
        task_design = np.nan_to_num(task_design)

    return task_design, lag_spline


def construct_physio_spline(physio, p_type, tr, nknots):
    if (p_type == 'rv_amp') | (p_type == 'rv_rate'):
        hrf = get_hrf(30, tr, 'rvt')
    elif p_type == 'hr':
        hrf = get_hrf(30, tr, 'hr')
    physio_conv = []
    for physio_subj in physio:
        physio_subj_conv = convolve_hrf(hrf, np.squeeze(physio_subj))
        physio_conv.append(physio_subj_conv)
    physio_conv = np.hstack(physio_conv)
    spline_basis = dmatrix("cr(x, df=nknots) - 1", {"x": physio_conv}, return_type='dataframe')
    return physio_conv, spline_basis


def evaluate_model(physio_eval, block_len, model, task_basis, physio_basis, task_eval=1):
    # Create task basis for model evaluation using previously defined design matrix
    block_lag = np.arange(block_len + 1)
    lag_pred = dmatrix(task_basis.design_info, {'x': block_lag}, return_type='dataframe')
    pred_list = [task_eval * lag_pred.iloc[:, l].values for l in range(lag_pred.shape[1])]
    task_pred = np.vstack(pred_list).T
    pred_maps = []
    n_physio = physio_basis.shape[1]
    n_task = task_basis.shape[1]
    for p in physio_eval:
        physio_pred = dmatrix(physio_basis.design_info, {'x': p},return_type='matrix')
        physio_pred = np.repeat(physio_pred, len(block_lag), axis=0)
        pred_design_mat = np.vstack([physio_pred[:,p]*task_pred[:,t] for p in range(n_physio) for t in range(n_task)]).T
        pred_design_mat = np.hstack([task_pred, physio_pred, pred_design_mat])
        # Get predictions from model
        pred_bold = model.predict(pred_design_mat)
        pred_maps.append(pred_bold)
    return pred_maps


def write_results(dataset, analysis, maps, zero_mask, n_vert, params, pred_vec=None, regressors=None, block_indx=None):
    if analysis == 'glm':
        analysis_str = f'{dataset}_task_glm'
        # if time-lag maps specified, get lag of maximum/minimum cross-correlation of each voxel.
        pickle.dump([pred_vec, regressors], open(f'{analysis_str}_results.pkl', 'wb'))
        for i, maps in enumerate(maps):
            write_nifti(maps, f'{analysis_str}_{i}', zero_mask, n_vert, params['mask'])
    elif analysis == 'block':
        analysis_str = f'{dataset}_task_block_avg'
        pickle.dump(block_indx, open(f'{analysis_str}_results.pkl', 'wb'))
        for map_n in maps:
            write_nifti(maps[map_n], f'{analysis_str}_{map_n}', zero_mask, n_vert, params['mask'])



def run_main(dataset, physio, physio_lag, analysis, block_lag, nknots):
    if analysis == 'glm':
        func_group_method = 'stack'
    elif analysis == 'block':
        func_group_method = 'list'

    func_data, physio_sig, physio_labels, zero_mask, n_vert, params = \
    load_data(dataset, 'group', physio=[physio], load_physio=True, group_method=func_group_method,
              physio_group_method='list') 

    # Create task regressors 
    if dataset == 'hcp_rel':
        subject_list = pd.read_csv(subject_list_hcp_rel)
        event_loader = fp_hcp_rel
        block_len = 22
    elif dataset == 'hcp_wm':
        subject_list = pd.read_csv(subject_list_hcp_wm)
        event_loader = fp_hcp_wm
        block_len = 35

    scan_len = len(physio_sig[0][0]) # not ideal
    task_reg, lag_spline = construct_hcp_task_regressors(dataset, block_len, subject_list, event_loader, 
                                                         params['tr'], scan_len, analysis) 

    if analysis == 'glm':
        # Create physio regressors
        physio_conv, physio_reg = construct_physio_spline(physio_sig[0], physio, params['tr'], nknots)

        # Construct Design matrix using patsy style formula
        print('construct spline matrix')
        design_mat = construct_interaction_terms(physio_reg, task_reg)

        # Run regression
        print('run regression')
        lin_reg = linear_regression(design_mat, func_data, return_model=True, 
                                    intercept=True, norm=False)

        # Get predicted maps at all time lags at equally spaced percentiles
        print('get predicted maps at different levels of task and physio')
        physio_eval = np.percentile(physio_conv, [5, 10, 50, 90, 95])
        pred_maps = evaluate_model(physio_eval, block_len, lin_reg, lag_spline, physio_reg)
        write_results(dataset, analysis, pred_maps, zero_mask, n_vert, params,
                      pred_vec = physio_eval, regressors = {'task': task_reg, 'physio': physio_conv})

    elif analysis == 'block':
        block_avg, block_indx = conditional_average_block(func_data, physio_sig[0], task_reg, block_lag)
        write_results(dataset, analysis, block_avg, zero_mask, n_vert, params, block_indx=block_indx)


if __name__ == '__main__':
    """Run main analysis"""
    parser = argparse.ArgumentParser(description='Interaction analysis of physio and task-evoked activity')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to run analysis on',
                        choices=['hcp_rel', 'hcp_wm'], 
                        required=True,
                        type=str)
    parser.add_argument('-p', '--physio',
                        help='select physio - can only provide one',
                        required=True,
                        choices=['hr', 'rv_amp', 'rv_rate'],
                        type=str)
    parser.add_argument('-l', '--lag',
                        help='Number of TRs to shift the physio signal in the positive (forward) direction',
                        required=False,
                        default=0,
                        type=int) 
    parser.add_argument('-a', '--analysis',
                        help='type of analysis to run: glm or block. The GLM anlaysis computes a glm-based interaction analysis of '
                        'physio and task blocks; The block-averaging performs a conditional block averaging of time courses based '
                        'on the level of the physio variable',
                        default='glm',
                        choices=['glm', 'block'],
                        type=str)
    parser.add_argument('-b', '--block_lag',
                        help='How many time points after the end of the block do you want included in the block-based averaging'
                        'analysis only for analysis: "block"',
                        default=10, 
                        required=False,
                        type=int)
    parser.add_argument('-n', '--nknots',
                        help='Number of knots in spline basis'
                        'Knots are placed at equally spaced quantiles based on N knots',
                        default=5, 
                        required=False,
                        type=int)


    args_dict = vars(parser.parse_args())
    run_main(args_dict['dataset'], args_dict['physio'], args_dict['lag'],
             args_dict['analysis'], args_dict['block_lag'], 
             args_dict['nknots'])

