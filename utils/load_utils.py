import json
import numpy as np
import os
import pandas as pd

subject_list_chang = 'data/dataset_chang/subject_list_chang.csv'
subject_list_choe = 'data/dataset_choe/run_list_choe.csv'
subject_list_gu = 'data/dataset_gu/subject_list_gu.csv'
subject_list_yale = 'data/dataset_yale/subject_list_yale.csv'
subject_list_nki = 'data/dataset_nki/subject_list_nki.csv'

physio_chang = ['eeg', 'hr', 'rv']
physio_choe = ['egg']



def find_fps(data, group, physio, subj_n=None, scan=None):
    subj_list = load_subject_list(data)
    if group:
        if data == 'chang':
            search_terms = subj_list[['subject', 'scan']].values.tolist()
        else: 
            search_terms = subj_list.subject.values.tolist()
    else:
        search_subj(data, subj_list, subj, scan)
        if scan is None:
            search_terms = [subj_n]
        else:
            search_terms = [subj_n, scan]
    if len(physio) > 1:
        search_physio(data, physio)
        physio.insert(0,'func')
    else:
        physio = ['func']

    if data == 'chang':
        fps = {d_type: [fp_chang(d_type, subj_scan[0],subj_scan[1]) for subj_scan in search_terms] 
               for d_type in physio}
    if data == 'choe':
        fps = {d_type: [fp_choe(d_type, run) for run in search_terms] 
               for d_type in physio}
    return fps


def fp_chang(data_type, subj, scan):
    if scan < 10:
        scan_str = f'000{scan}'
    else:
        scan_str = f'00{scan}'

    if data_type == 'func':
        f_str = f'data/dataset_chang/func/proc3_filter_norm/sub_00{subj}-mr_{scan_str}-ecr_echo1_w_dspk_blur3mm.nii'
    elif data_type == 'eeg':
        f_str = f'data/dataset_chang/eeg/sub_00{subj}-mr_{scan_str}-ecr_echo1_eeg_at.mat'
    elif data_type == 'hr':
        f_str = f'data/dataset_chang/physio_hr/sub_00{subj}-mr_{scan_str}-ecr_echo1_hr.mat'
    elif data_type == 'rv':
        f_str = f'data/dataset_chang/physio_rv/sub_00{subj}-mr_{scan_str}-ecr_echo1_rv.mat'
    return f_str


def fp_choe(data_type, run):
    if data_type == 'func':
        f_str = f'data/dataset_choe/func/proc3_filter/pb04.20190{run}jp.r01.blur+tlrc.nii'
    elif data_type == 'egg':
        f_str = f'data/dataset_choe/egg/proc2_resample/0{run}.txt'
    return f_str


def load_subject_list(data):
    if data == 'chang':
        subj_list = pd.read_csv(subject_list_chang)
    elif data == 'choe':
        subj_list = pd.read_csv(subject_list_choe)
    return subj_list


def search_physio(data, physio):
    import pdb; pdb.set_trace()
    if data == 'chang':
        found = all([p in physio_chang for p in physio])
    elif data == 'choe':
        found = all([p in physio_chang for p in physio])
    if ~found:
        raise Exception("""
                        Dataset "{0}" does not contain one or more supplied physio terms. See
                        top of utils/load_utils.py for acceptable physio labels for each dataset.
                        """.format(data)
                        )


def search_subj(data, subj_list, subj, scan=None):
    if scan is not None:
        found = subj in subj_list.subject.values
    if scan is not None:
        found = (subj in subj_list.subject.values) & (scan in subj_list.scan.values)
    # Ensure there is only one file associated with the supplied subject number (and scan number - if supplied)
    if ~found:
        raise Exception(f'No file found in dataset "{data}" associated with input subject label {subj} and/or scan number')



    
    