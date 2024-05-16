import argparse
import json
import nibabel as nb
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
block_len = 15 # length of chang task blocks (30s post auditory tone)
chang_trim = 14.7 # remove first 14.7s to account for first 7 volumes removed

# load params file (used to get TR and subject lists)
data_params = json.load(open('analysis_params.json', 'rb'))


def construct_task_blocks(dataset, subject, scan, params_d):
    # loop through subjects and construct task blocks for averaging
    event_blocks = []
    for subj, sc in zip(subject, scan):
        if dataset == 'chang_bh':
            events = load_chang_bh_event_file()
        elif dataset == 'chang_cue': 
            events, rt = load_chang_cue_event_file(subj, sc)
            # filter out events with no reaction time response (convert ms to s)
            events = np.array([e/1000 for e, rt in zip(events, rt) if ~np.isnan(rt)])

        # get TR indices of task blocks
        events = events - chang_trim # 14.7 secs
        # Convert event secs to TRs
        events_tr = events/params_d['tr']
        events_tr = np.round(events_tr).astype(int)
        # Construct blocks
        blocks = []
        for event in events_tr:
            # include TR before onset
            blocks.append(np.arange(event-1, event+block_len))
        # filter out overlapping and 
        if dataset == 'chang_cue':
            blocks = filter_chang_cue_blocks(blocks, subj, sc, params_d)
        event_blocks.append(blocks)
    return event_blocks


# function to detect overlapping trials
def detect_overlap(indx, intervals):
    intervals_c = intervals.copy()
    k_range = set(intervals[indx])
    del intervals_c[indx]
    if any([len(k_range.intersection(k))>2 for k in intervals_c]):
        return True
    else:
        return False


def event_index(func_data, event_blocks, compliance=None):
    if compliance is not None:
        event_blocks = [block for block, c in zip(event_blocks, compliance) if c == 1]
    return func_data[event_blocks, :]


def event_average(func_event_blocks):
    return func_event_blocks.mean(axis=0)


def filter_chang_cue_blocks(blocks, subj, scan, params):
    """
    some auditory tones in the chang cue dataset are too close together 
    for reliabile separation of overlapping hemodynamic responses
    """
    # load number of volumes in functional scan
    func_len = nb.load(params['func'].format(subj, scan)).shape[-1]
    # filter out trials occuring before (trimmed) func start, after end of func
    blocks_filt = [b for i, b in enumerate(blocks)
                   if (b[-1] < func_len) and (b[0] >= 0) 
                   and (not detect_overlap(i, blocks))]
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


def write_results(dataset, func_avg, params, out_dir):
    if out_dir is not None:
        analysis_str = f'{out_dir}/{dataset}_task_avg'
    else:
        analysis_str = f'{dataset}_task_avg'

    write_nifti(func_avg, analysis_str, params)


def run_task_avg(dataset, m_param, out_dir=None):
    # get dataset parameters
    params_d = data_params[dataset] 
    # load subject list 
    subject, scan = load_subject_list(dataset, params_d['subject_list'])
    # Load data
    func_data, _, params = load_data(dataset, None, group_method='list', 
                                     multiecho=m_param)
    # create task block indices
    event_blocks = construct_task_blocks(dataset, subject, scan, params_d)
    # load compliance for chang_bh
    if dataset == 'chang_bh':
        compliance = load_chang_compliance()
    # loop through blocks types and subject functional scans and do:
    # segment into task blocks, stack across subjects and average
    # if chang_bh, use compliance text to remove task blocks where there was
    # no compliance (i.e. no 'deep breath')
    func_blocks = []
    block_iter = zip(func_data, subject, scan, event_blocks)
    for i, (func, subj, sc, block_s) in enumerate(block_iter):
        if dataset == 'chang_bh':
            compliance_s = compliance[f'{subj}_{sc}']
        else:
            compliance_s = np.ones(len(block_s))
        for block, c in zip(block_s, compliance_s):
            if c == 1:
                func_blocks.append(
                    func[block, :]
                )
    # stack all blocks across subjects
    func_blocks_array = np.stack(func_blocks)
    # average across all blocks
    func_blocks_avg = func_blocks_array.mean(axis=0)
    # write results
    write_results(dataset, func_blocks_avg, params, out_dir)
        
        
if __name__ == '__main__':
    """Run main analysis"""
    parser = argparse.ArgumentParser(description='Event averaging of task fMRI datasets')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to run analysis on',
                        choices=['chang_bh', 'chang_cue'], 
                        required=True,
                        type=str)
    parser.add_argument('-m', '--m_param',
                        help='For multiecho data, specify multiecho parameter - kappa or rho',
                        choices=['kappa', 'rho'], 
                        required=False,
                        default=None,
                        type=str)

    args_dict = vars(parser.parse_args())
    run_task_avg(args_dict['dataset'], args_dict['m_param'])

