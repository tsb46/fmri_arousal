import argparse
import numpy as np
import pandas as pd
import pickle

from scipy.stats import zscore
from utils.load_write import load_data, write_nifti, load_chang_bh_event_file



def construct_deep_breath_blocks(events, tr, block_len=15, trim=True):
    # remove first 14.7s to account for first 7 volumes removed
    if trim:
        events = events - 14.7 #secs
    # Convert event secs to TRs
    events_tr = events/tr
    events_tr = np.round(events_tr).astype(int)
    # Construct blocks
    event_blocks = []
    for event in events_tr:
        event_blocks.append(np.arange(event, event+block_len))

    return event_blocks

def event_index(func_data, event_blocks):
    return func_data[event_blocks, :]

def event_average(func_event_blocks):
    return func_event_blocks.mean(axis=0)


def write_results(dataset, func_avg, level, subj_n, scan, zero_mask, n_vert):
    if level == 'subject':
        analysis_str = f'{dataset}_taskbh_{subj_n}_{scan}'
    else:
        analysis_str = f'{dataset}_taskbh_group'
        pickle.dump(func_avg, open(f'{analysis_str}_results.pkl', 'wb'))

    write_nifti(func_avg, analysis_str, zero_mask, n_vert)



def run_main(dataset, level):
    # Load data
    func_data, _, _, zero_mask, n_vert, params = load_data(dataset, level, physio=None, load_physio=False, 
                                                           group_method='list') 
    # load Chang BH event file (assuming timing is the same across subjects
    bh_events = load_chang_bh_event_file()
    bh_event_blocks = construct_deep_breath_blocks(bh_events, params['tr'])
    # load subject list for mapping subj and scan # to output names
    subject_list = pd.read_csv('data/dataset_chang_bh/subject_list_chang_bh.csv')

    func_block_all = []
    for i, func in enumerate(func_data):
        func_blocks = event_index(func, bh_event_blocks)
        if level == 'subject':
            func_avg = event_average(func_blocks)
            subj_n, scan_n = subject_list.iloc[i,:][['subject', 'scan']]
            write_results(dataset, func_avg, level, subj_n, scan_n, zero_mask, n_vert)
        else:
            func_block_all.append(func_blocks)

    if level == 'group':
        func_block_all_array = np.vstack(func_block_all)
        func_all_avg = func_block_all_array.mean(axis=0)
        write_results(dataset, func_all_avg, level, None, None, zero_mask, n_vert)



if __name__ == '__main__':
    """Run main analysis"""
    parser = argparse.ArgumentParser(description='Event averaging of Chang deep breath events')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to run analysis on - only option is Chang_BH',
                        choices=['chang_bh'], 
                        default='chang_bh',
                        required=False,
                        type=str)
    parser.add_argument('-l', '--level',
                        help='subject or group level analysis',
                        default='group',
                        choices=['subject', 'group'],
                        type=str)


    args_dict = vars(parser.parse_args())
    run_main(args_dict['dataset'], args_dict['level'])

