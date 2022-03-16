import json
import numpy as np
import os
import pandas as pd

subject_list_chang = 'data/dataset_chang/subject_list_chang.csv'
subject_list_yale = 'data/dataset_yale/subject_list_yale_subset.csv'
subject_list_nki = 'data/dataset_nki/subject_list_nki.csv'
subject_list_hcp = 'data/dataset_hcp/subject_list_hcp_subset_test.csv'
subject_list_lemon = 'data/dataset_lemon/subject_list_lemon.csv'


def find_fps(data, level, physio, params, events=False, subj_n=None, scan=None):
    subj_list = load_subject_list(data)
    physio_fp = physio.copy()
    if level == 'group':
        if data == 'chang':
            search_terms = subj_list[['subject', 'scan']].values.tolist()
        elif data == 'hcp':
            search_terms = subj_list[['subject', 'lr']].values.tolist()
        else:
            search_terms = subj_list.subject.values.tolist()
    else:
        search_subj(data, subj_list, subj_n, scan)
        if scan is None:
            if data == 'chang':
                scan_chang = subj_list.loc[subj_list.subject == subj_n, 'scan'].values[0]
                search_terms = [[subj_n, scan_chang]] 
            elif data == 'hcp':
                scan_hcp = subj_list.loc[subj_list.subject == subj_n, 'lr'].values[0]
                search_terms = [[subj_n, scan_hcp]]
            else:
                search_terms = [subj_n]
        else:
            search_terms = [[subj_n, scan]]
    if len(physio) > 0:
        search_physio(data, physio_fp, params['physio'])
        physio_fp.insert(0,'func')
    else:
        physio_fp = ['func']

    if events:
        physio_fp.insert(0, 'events')

    if data == 'chang':
        fps = {d_type: [fp_chang(d_type, subj_scan[0],subj_scan[1]) for subj_scan in search_terms] 
               for d_type in physio_fp}
    elif data == 'nki':
        fps = {d_type: [fp_nki(d_type, subj) for subj in search_terms] 
               for d_type in physio_fp}
    elif data == 'yale':
        fps = {d_type: [fp_yale(d_type, subj) for subj in search_terms] 
               for d_type in physio_fp}
    elif data == 'hcp':
        fps = {d_type: [fp_hcp(d_type, subj_scan[0],subj_scan[1]) for subj_scan in search_terms] 
               for d_type in physio_fp}
    elif data == 'lemon':
        fps = {d_type: [fp_lemon(d_type, subj) for subj in search_terms] 
               for d_type in physio_fp}
    return fps


def fp_chang(data_type, subj, scan):
    if scan < 10:
        scan_str = f'000{scan}'
    else:
        scan_str = f'00{scan}'

    if data_type == 'func':
        f_str = f'data/dataset_chang/func/proc4_bandpass/sub_00{subj}-mr_{scan_str}-ecr_echo1_w_dspk_blur3mm.nii.gz' 
    elif data_type == 'alpha':
        f_str = f'data/dataset_chang/eeg/proc1_fbands/sub_00{subj}_mr_{scan_str}_fbands_Alpha.txt'
    elif data_type == 'delta':
        f_str = f'data/dataset_chang/eeg/proc1_fbands/sub_00{subj}_mr_{scan_str}_fbands_Delta.txt'
    elif data_type == 'infraslow':
        f_str = f'data/dataset_chang/eeg/proc1_fbands/sub_00{subj}_mr_{scan_str}_fbands_Infraslow.txt'
    elif data_type == 'hr':
        f_str = f'data/dataset_chang/physio/proc1_physio/sub_00{subj}_mr_{scan_str}_physio_PPG_RATE_NK.txt'
    elif data_type == 'rv':
        f_str = f'data/dataset_chang/physio/proc1_physio/sub_00{subj}_mr_{scan_str}_physio_RESP_AMP_HILBERT.txt'
    elif data_type == 'csf':
        f_str = f'data/dataset_chang/physio/raw_csf/sub_00{subj}_mr_{scan_str}.txt'
    elif data_type == 'vigilance':
        f_str = f'data/dataset_chang/eeg/proc1_fbands/sub_00{subj}_mr_{scan_str}_fbands_vigilance.txt'
    elif data_type == 'pc1':
        f_str = f'data/dataset_chang/physio/proc2_physio_pc/sub_00{subj}-mr_{scan_str}_PC1_angle.txt'
    elif data_type == 'pc2':
        f_str = f'data/dataset_chang/physio/proc2_physio_pc/sub_00{subj}-mr_{scan_str}_PC2_angle.txt'
    elif data_type == 'pc3':
        f_str = f'data/dataset_chang/physio/proc2_physio_pc/sub_00{subj}-mr_{scan_str}_PC3_angle.txt'
    return f_str


def fp_hcp(data_type, subj, scan):
    if data_type == 'func':
        f_str = f'data/dataset_hcp/func/proc3_filter_norm/{subj}_{scan}1_rest.nii.gz'
    elif data_type == 'rv':
        f_str = f'data/dataset_hcp/physio/proc1_physio/{subj}_physio_RESP_VAR_W.txt'
    elif data_type == 'hr':
        f_str = f'data/dataset_hcp/physio/proc1_physio/{subj}_physio_PPG_HR_W.txt'
    return f_str


def fp_lemon(data_type, subj):
    if data_type == 'func':
        f_str = f'data/dataset_lemon/func/proc6_trim/{subj}_task_rest.nii.gz'
    elif data_type == 'rv':
        f_str = f'data/dataset_lemon/physio/proc1_physio/{subj}_physio_RESP_VAR_W.txt'
    elif data_type == 'hr':
        f_str = f'data/dataset_lemon/physio/proc1_physio/{subj}_physio_ECG_HR_NK.txt'
    elif data_type == 'map':
        f_str = f'data/dataset_lemon/physio/proc1_physio/{subj}_physio_MAP_W.txt'
    elif data_type == 'csf':
        f_str =f'data/dataset_lemon/physio/proc1_physio/{subj}_physio_csf.txt'
    return f_str


def fp_nki(data_type, subj):
    if data_type == 'func':
        f_str = f'data/dataset_nki/func/proc5_filter_norm/{subj}_task_breathhold.nii.gz'
    elif data_type == 'hr':
        f_str = f'data/dataset_nki/physio/proc1_physio/{subj}_task_breathhold_physio_hr.txt'
    elif data_type == 'rv':
        f_str = f'data/dataset_nki/physio/proc1_physio/{subj}_task_breathhold_physio_rv.txt'
    elif data_type == 'csf':
        f_str = f'data/dataset_nki/physio/proc1_physio/{subj}_task_breathhold_physio_csf.txt'
    return f_str


def fp_yale(data_type, subj):
    if data_type == 'func':
        f_str = f'data/dataset_yale/func/proc6_trim/{subj}_task-rest_run-01_bold.nii.gz'
    if data_type == 'pupil':
        f_str = f'data/dataset_yale/pupillometry/{subj}_task-rest_run-01_et.tsv'
    return f_str


def load_subject_list(data):
    if data == 'chang':
        subj_list = pd.read_csv(subject_list_chang)
    elif data == 'nki':
        subj_list = pd.read_csv(subject_list_nki)
    elif data == 'yale':
        subj_list = pd.read_csv(subject_list_yale)
    elif data == 'hcp':
        subj_list = pd.read_csv(subject_list_hcp)
    elif data == 'lemon':
        subj_list = pd.read_csv(subject_list_lemon)
    return subj_list


def print_filter_info(params, load_physio):
    # print filter parameters of functional data
    if params['data']['func']['filter_params']['filter_choice'] == 'raw':
        print(f'no filtering applied to functional data \n')
    else:
        filter_choice = params['data']['func']['filter_params']['filter_choice'] 
        if filter_choice == 'bandpass':
            low = params['data']['func']['filter_params'][filter_choice]['low'] 
            high = params['data']['func']['filter_params'][filter_choice]['high']
        elif filter_choice == 'lowpass':
            low = None
            high = params['data']['func']['filter_params'][filter_choice]['high']
        elif filter_choice == 'highpass':
            low = params['data']['func']['filter_params'][filter_choice]['low']
            high = None
        print(f'{filter_choice} filtering applied to functional data: {low} - {high} Hz \n')
    if load_physio:
        # print filter parameters for physio signals
        for p in params['physio']:
            if params['data']['physio']['filter_params']['filter_choice'][p] == 'raw':
                print(f'no filtering applied to {p} data \n')
            else:
                filter_choice = params['data']['physio']['filter_params']['filter_choice'][p]
                if filter_choice == 'bandpass':
                    low = params['data']['physio']['filter_params'][filter_choice]['low'] 
                    high = params['data']['physio']['filter_params'][filter_choice]['high']
                elif filter_choice == 'lowpass':
                    low = None
                    high = params['data']['physio']['filter_params'][filter_choice]['high']
                elif filter_choice == 'highpass':
                    low = params['data']['physio']['filter_params'][filter_choice]['low']
                    high = None
                print(f'{filter_choice} filtering applied to {p} data: {low} - {high} Hz \n')


def search_physio(data, physio, physio_standard):
    found = all([p in physio_standard for p in physio])
    if not found:
        raise Exception("""
                        Dataset "{0}" does not contain one or more supplied physio terms. See
                        top of utils/load_utils.py for acceptable physio labels for each dataset.
                        """.format(data)
                        )


def search_subj(data, subj_list, subj, scan=None):
    if scan is None:
        found = subj in subj_list.subject.values
        if (subj_list.values==subj).sum() > 1:
            raise Exception(f'multiples files found in dataset "{data}" associated with input subject label '
                            '{subj} - please provide scan number')
    else:
        found = (subj in subj_list.subject.values) & (scan in subj_list.scan.values)
    # Ensure there is only one file associated with the supplied subject number (and scan number - if supplied)
    if not found:
        raise Exception(f'No file found in dataset "{data}" associated with input subject label {subj} and/or scan number')



    
    