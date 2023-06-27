import argparse
import json
import numpy as np
import pandas as pd
import pickle

from scipy.stats import zscore
from utils.load_write import (
    load_data, write_nifti, 
    load_chang_bh_event_file, load_nki_event_file,
    load_chang_cue_event_file, load_subject_list
)
from utils.glm_utils import onsets_to_block


# Global parameters
chang_block_len = 15 # length of chang task blocks (30s post auditory tone)
chang_trim = 14.7 # remove first 14.7s to account for first 7 volumes removed

# load params file (used to get TR and subject lists)
params = json.load(open('analysis_params.json', 'rb'))

# get nki and chang_bh TR
nki_tr = params['nki']['tr']
chang_tr = params['chang_bh']['tr'] # chang_bh and chang_cue have same TR


def construct_task_blocks(dataset, events, expand):
    # get TR indices of task blocks
    if dataset in ['chang_bh', 'chang_cue']:
        events = events - chang_trim # 14.7 secs
        # Convert event secs to TRs
        events_tr = events/chang_tr
        events_tr = np.round(events_tr).astype(int)
        # Construct blocks
        event_blocks = []
        for event in events_tr:
            event_blocks.append(np.arange(event, event+chang_block_len+expand))
        if dataset == 'chang_bh':
            event_blocks = [('breath', event_blocks)]
        elif dataset == 'chang_cue':
            event_blocks = filter_chang_cue_blocks(event_blocks)
            event_blocks = [('cue', event_blocks)]
    elif dataset == 'nki':
        event_df = group_nki_blocks(events)
        group_func = {'onset': 'first', 'duration': 'sum'}
        # Construct 'deep breath' blocks
        breath_df = event_df.loc[event_df.breath]
        breath_df_grp = breath_df.groupby(['breath_indx']).agg(group_func)
        breath_blocks = onsets_to_block(breath_df_grp, nki_tr, expand)
        # Construct 'breath hold' blocks
        hold_df = event_df.loc[event_df.hold]
        hold_df_grp = hold_df.groupby(['hold_indx']).agg(group_func)
        hold_blocks = onsets_to_block(hold_df_grp, nki_tr, expand)
        event_blocks = [('rest-breath', breath_blocks), ('hold', hold_blocks)]
    return event_blocks


def event_index(func_data, event_blocks, compliance=None):
    if compliance is not None:
        event_blocks = [block for block, c in zip(event_blocks, compliance) if c == 1]
    return func_data[event_blocks, :]


def event_average(func_event_blocks):
    return func_event_blocks.mean(axis=0)


def filter_chang_cue_blocks(blocks):
    """
    some auditory tones in the chang cue dataset are too close together 
    for reliabile separation of overlapping hemodynamic responses
    """
    block_prev = []
    blocks_filt = []
    for block in blocks:
        if (len(set(block) & set(block_prev)) == 0) & (block[-1] <= 692):
            blocks_filt.append(block)
        block_prev = block
    return blocks_filt


def load_chang_compliance():
     # Load trial level compliance file
    compliance = pd.read_csv('data/dataset_chang_bh/compliance.csv')
    trial_cols = [f'trial{n+1}' for n in range(9)]
    # Convert to dict
    compliance_dict = {}
    for i, (subj, scan) in enumerate(zip(compliance.subject, compliance.scan)): 
        if scan < 10:
            scan_str = f'000{scan}'
        else:
            scan_str = f'00{scan}'
        compliance_dict[f'{subj}_{scan_str}'] = compliance.iloc[i, :][trial_cols].values
    return compliance_dict


def group_nki_blocks(event_df, merge_rest_breath=True):
    # group consecutive trials of nki breathhold task
    event_df['hold'] = event_df.trial_type.str.startswith('H')
    event_df['hold_indx'] = (event_df.hold != event_df.hold.shift(1)).cumsum()
    if merge_rest_breath:
        rest_breath_block = ['R', 'G', 'Deep', 'In', 'Out']
        event_df['breath'] = event_df.trial_type.isin(rest_breath_block)
        event_df['breath_indx'] = (event_df.breath != event_df.breath.shift(1)).cumsum()
    else:
        breath_block = ['Deep', 'In', 'Out']
        event_df['breath'] = event_df.trial_type.isin(breath_block)
        event_df['breath_indx'] = (event_df.breath != event_df.breath.shift(1)).cumsum()
        rest_block = ['R', 'G']
        event_df['rest'] = event_df.trial_type.isin(rest_block)
        event_df['rest_indx'] = (event_df.rest != event_df.rest.shift(1)).cumsum()
    return event_df


def write_results(dataset, label, func_avg, zero_mask, n_vert, out_dir):
    if out_dir is not None:
        analysis_str = f'{out_dir}/{dataset}_task_avg_{label}'
    else:
        analysis_str = f'{dataset}_task_avg_{label}'

    pickle.dump(func_avg, open(f'{analysis_str}_results.pkl', 'wb'))
    write_nifti(func_avg, analysis_str, zero_mask, n_vert)


def run_task_avg(dataset, expand, out_dir=None):
    # Load data
    func_data, _, zero_mask, n_vert = load_data(dataset, None, group_method='list')
    # load event file (assuming event timings are the same across subjects)
    if dataset == 'nki':
        events = load_nki_event_file()
    elif dataset == 'chang_cue':
        events = load_chang_cue_event_file()
    elif dataset == 'chang_bh':
        events = load_chang_bh_event_file()
        compliance = load_chang_compliance()

    event_blocks = construct_task_blocks(dataset, events, expand)
    # load subject list for mapping functional to compliance (only used for chang_bh)
    subject, scan = load_subject_list(dataset, params[dataset]['subject_list'])

    # loop through blocks types and subject functional scans and do:
    # segment into task blocks, stack across subjects and average
    # if chang_bh, use compliance text to remove task blocks where there was
    # no compliance (i.e. no 'deep breath')
    for block in event_blocks:
        func_blocks = []
        for i, (func, subj, run) in enumerate(zip(func_data, subject, scan)):
            if dataset == 'chang_bh':
                compliance_subj = compliance[f'{subj}_{run}']
            else:
                compliance_subj = None
            func_segment = event_index(func, block[1], compliance_subj)
            func_blocks.append(func_segment)
        # stack all blocks across subjects
        func_blocks_array = np.vstack(func_blocks)
        # average across all blocks
        func_blocks_avg = func_blocks_array.mean(axis=0)
        # write results
        write_results(dataset, block[0], func_blocks_avg, zero_mask, n_vert, out_dir)
        
        
if __name__ == '__main__':
    """Run main analysis"""
    parser = argparse.ArgumentParser(description='Event averaging of task fMRI datasets')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to run analysis on',
                        choices=['chang_bh', 'chang_cue', 'nki'], 
                        required=True,
                        type=str)
    parser.add_argument('-e', '--expand',
                        help='how many TRs to expand the end of the task block '
                        'to account for lag of HRF',
                        required=False,
                        default=0,
                        type=int)


    args_dict = vars(parser.parse_args())
    run_task_avg(args_dict['dataset'], args_dict['expand'])

