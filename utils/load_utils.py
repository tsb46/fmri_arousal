import json
import numpy as np
import os
import pandas as pd

subject_list_chang = 'data/dataset_chang/subject_list_chang.csv'
subject_list_chang_bh = 'data/dataset_chang_bh/subject_list_chang_bh.csv'
subject_list_nki = 'data/dataset_nki/subject_list_nki.csv'
subject_list_hcp = 'data/dataset_hcp/subject_list_hcp_rest.csv'
subject_list_hcp_rel = 'data/dataset_hcp_task/subject_list_hcp_relational.csv'
subject_list_hcp_wm = 'data/dataset_hcp_task/subject_list_hcp_wm.csv'
subject_list_monash = 'data/dataset_monash/subject_list_monash_subset.csv'
subject_list_spreng = 'data/dataset_spreng/subject_list_spreng.csv'


def find_fps(data, level, physio, params, subj_n=None, scan=None):
    subj_list = load_subject_list(data)
    physio_fp = physio.copy()
    if level == 'group':
        if (data == 'chang') | (data == 'chang_bh'):
            search_terms = subj_list[['subject', 'scan']].values.tolist()
        elif (data == 'hcp') | (data == 'hcp_fix') | (data == 'hcp_rel') | (data == 'hcp_wm'):
            search_terms = subj_list[['subject', 'lr']].values.tolist()
        else:
            search_terms = subj_list.subject.values.tolist()
    else:
        search_subj(data, subj_list, subj_n, scan)
        if scan is None:
            if (data == 'chang') | (data == 'chang_bh'):
                scan_chang = subj_list.loc[subj_list.subject == subj_n, 'scan'].values[0]
                search_terms = [[subj_n, scan_chang]] 
            elif (data == 'hcp') | (data == 'hcp_fix') | (data == 'hcp_rel') | (data == 'hcp_wm'):
                scan_hcp = subj_list.loc[subj_list.subject == subj_n, 'lr'].values[0]
                search_terms = [[subj_n, scan_hcp]]
            else:
                search_terms = [subj_n]
        else:
            search_terms = [[subj_n, scan]]

    if len(physio) > 0:
        search_physio(data, physio_fp, params['physio'])
        physio_fp.insert(0, 'func')
    else:
        physio_fp = ['func']


    if data == 'chang':
        fps = {d_type: [fp_chang(d_type, subj_scan[0],subj_scan[1]) for subj_scan in search_terms] 
               for d_type in physio_fp}
    elif data == 'chang_bh':
        fps = {d_type: [fp_chang_bh(d_type, subj_scan[0],subj_scan[1]) for subj_scan in search_terms] 
               for d_type in physio_fp}
    elif data == 'nki':
        fps = {d_type: [fp_nki(d_type, subj) for subj in search_terms] for d_type in physio_fp}
    elif (data == 'hcp') | (data == 'hcp_fix'):
        if data == 'hcp_fix':
            fix = True
        else:
            fix = False
        fps = {d_type: [fp_hcp(d_type, subj_scan[0],subj_scan[1], fix) for subj_scan in search_terms] 
               for d_type in physio_fp}
    elif data == 'hcp_rel':
        fps = {d_type: [fp_hcp_rel(d_type, subj_scan[0],subj_scan[1]) for subj_scan in search_terms] 
               for d_type in physio_fp}
    elif data == 'hcp_wm':
        fps = {d_type: [fp_hcp_wm(d_type, subj_scan[0],subj_scan[1]) for subj_scan in search_terms] 
               for d_type in physio_fp}
    elif (data == 'monash') | (data == 'monash_pet'):
        if data ==  'monash_pet':
            modality='pet'
        else:
            modality='func'
        fps = {d_type: [fp_monash(d_type, subj, modality) for subj in search_terms] 
               for d_type in physio_fp}
    elif data == 'spreng':
        fps = {d_type: [fp_spreng(d_type, subj) for subj in search_terms] for d_type in physio_fp}
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
        # f_str = f'data/dataset_chang/physio/proc1_physio/sub_00{subj}_mr_{scan_str}_physio_RESP_AMP_HILBERT.txt'
        f_str = f'data/dataset_chang/physio/proc1_physio/sub_00{subj}_mr_{scan_str}_physio_RESP_RVT_NK.txt'
    elif data_type == 'csf':
        f_str = f'data/dataset_chang/physio/raw_csf/sub_00{subj}_mr_{scan_str}.txt'
    elif data_type == 'vigilance':
        f_str = f'data/dataset_chang/eeg/proc1_fbands/sub_00{subj}_mr_{scan_str}_fbands_vigilance.txt'
    elif data_type == 'ppg_low':
        f_str = f'data/dataset_chang/physio/proc1_physio/sub_00{subj}_mr_{scan_str}_physio_PPG_LOW_NK.txt'
    elif data_type == 'precuneus':
        f_str = f'data/dataset_chang/physio/proc1_physio/sub_00{subj}_mr_{scan_str}_precuneus.txt'
    elif data_type == 'superior_parietal':
        f_str = f'data/dataset_chang/physio/proc1_physio/sub_00{subj}_mr_{scan_str}_superior_parietal.txt'
    elif data_type == 'global_sig':
        f_str = f'data/dataset_chang/physio/proc1_physio/sub_00{subj}_mr_{scan_str}_global_sig.txt'

    return f_str


def fp_chang_bh(data_type, subj, scan):
    if scan < 10:
        scan_str = f'000{scan}'
    else:
        scan_str = f'00{scan}'

    if data_type == 'func':
        f_str = f'data/dataset_chang_bh/func/proc4_bandpass/sub_00{subj}-mr_{scan_str}-adb_echo1_w_dspk_blur3mm.nii.gz' 
    elif data_type == 'alpha':
        f_str = f'data/dataset_chang_bh/eeg/proc1_fbands/sub_00{subj}_mr_{scan_str}_fbands_Alpha.txt'
    elif data_type == 'delta':
        f_str = f'data/dataset_chang_bh/eeg/proc1_fbands/sub_00{subj}_mr_{scan_str}_fbands_Delta.txt'
    elif data_type == 'infraslow':
        f_str = f'data/dataset_chang_bh/eeg/proc1_fbands/sub_00{subj}_mr_{scan_str}_fbands_Infraslow.txt'
    elif data_type == 'hr':
        f_str = f'data/dataset_chang_bh/physio/proc1_physio/sub_00{subj}_mr_{scan_str}_physio_PPG_RATE_NK.txt'
    elif data_type == 'rv':
        f_str = f'data/dataset_chang_bh/physio/proc1_physio/sub_00{subj}_mr_{scan_str}_physio_RESP_AMP_HILBERT.txt'
    elif data_type == 'csf':
        f_str = f'data/dataset_chang_bh/physio/raw_csf/sub_00{subj}_mr_{scan_str}.txt'
    elif data_type == 'vigilance':
        f_str = f'data/dataset_chang_bh/eeg/proc1_fbands/sub_00{subj}_mr_{scan_str}_fbands_vigilance_at.txt'
    elif data_type == 'ppg_low':
        f_str = f'data/dataset_chang_bh/physio/proc1_physio/sub_00{subj}_mr_{scan_str}_physio_PPG_LOW_NK.txt'
    elif data_type == 'global_sig':
        f_str = f'data/dataset_chang_bh/physio/proc1_physio/sub_00{subj}_mr_{scan_str}_global_sig.txt'

    return f_str


def fp_hcp(data_type, subj, scan, fix):
    if data_type == 'func':
        if fix:
            f_str = f'data/dataset_hcp/func_fix/proc4_bandpass/{subj}_{scan}1_rest.nii.gz'
        else:
            f_str = f'data/dataset_hcp/func/proc4_bandpass/{subj}_{scan}1_rest.nii.gz'
    elif (data_type == 'rv') | (data_type == 'rv_amp'):
        f_str = f'data/dataset_hcp/physio/proc1_physio/{subj}_physio_RESP_AMP_HILBERT.txt'
    elif data_type == 'rv_rate':
        f_str = f'data/dataset_hcp/physio/proc1_physio/{subj}_physio_RESP_RATE_NK.txt'
    elif data_type == 'hr':
        f_str = f'data/dataset_hcp/physio/proc1_physio/{subj}_physio_PPG_HR_NK.txt'
    elif data_type == 'ppg_low':
        f_str = f'data/dataset_hcp/physio/proc1_physio/{subj}_physio_PPG_LOW_NK.txt'
    elif data_type == 'precuneus':
        f_str = f'data/dataset_hcp/physio/proc1_physio/{subj}_precuneus.txt'
    elif data_type == 'superior_parietal':
        f_str = f'data/dataset_hcp/physio/proc1_physio/{subj}_superior_parietal.txt'
    elif data_type == 'global_sig':
        f_str = f'data/dataset_hcp/physio/proc1_physio/{subj}_global_sig.txt'
    return f_str


def fp_hcp_rel(data_type, subj, scan):
    if data_type == 'func':
        f_str = f'data/dataset_hcp_task/func_rel/proc4_bandpass/{subj}_{scan}_relational.nii.gz'
    elif (data_type == 'rv') | (data_type == 'rv_amp'):
        f_str = f'data/dataset_hcp_task/physio_rel/proc1_physio/{subj}_physio_RESP_RVT_NK.txt'
    elif data_type == 'rv_rate':
        f_str = f'data/dataset_hcp_task/physio_rel/proc1_physio/{subj}_physio_RESP_RATE_NK.txt'
    elif data_type == 'hr':
        f_str = f'data/dataset_hcp_task/physio_rel/proc1_physio/{subj}_physio_PPG_HR_NK.txt'
    elif data_type == 'ppg_low':
        f_str = f'data/dataset_hcp_task/physio_rel/proc1_physio/{subj}_physio_PPG_LOW_NK.txt'
    elif data_type == 'events':
        f_str = f'data/dataset_hcp_task/events_rel/{subj}_{scan}_EV'
    elif data_type in ('precuneus', 'superior_parietal', 'global_sig'):
        f_str = None
    return f_str


def fp_hcp_wm(data_type, subj, scan):
    if data_type == 'func':
        f_str = f'data/dataset_hcp_task/func_wm/proc4_bandpass/{subj}_{scan}_wm.nii.gz'
    elif (data_type == 'rv') | (data_type == 'rv_amp'):
        f_str = f'data/dataset_hcp_task/physio_wm/proc1_physio/{subj}_physio_RESP_RVT_NK.txt'
    elif data_type == 'rv_rate':
        f_str = f'data/dataset_hcp_task/physio_wm/proc1_physio/{subj}_physio_RESP_RATE_NK.txt'
    elif data_type == 'hr':
        f_str = f'data/dataset_hcp_task/physio_wm/proc1_physio/{subj}_physio_PPG_HR_NK.txt'
    elif data_type == 'ppg_low':
        f_str = f'data/dataset_hcp_task/physio_wm/proc1_physio/{subj}_physio_PPG_LOW_NK.txt'
    elif data_type == 'events':
        f_str = f'data/dataset_hcp_task/events_wm/{subj}_{scan}_EV'
    elif data_type in ('precuneus', 'superior_parietal', 'global_sig'):
        f_str = None
    return f_str


def fp_monash(data_type, subj, modality):
    if data_type == 'func':
        if modality == 'pet':
            f_str = f'data/dataset_monash/pet/proc5_detrend/{subj}_pet.nii.gz'
        else:
            f_str = f'data/dataset_monash/func/proc7_bandpass/{subj}.nii.gz'
    return f_str


def fp_nki(data_type, subj):
    if data_type == 'func':
        f_str = f'data/dataset_nki/func/proc5_filter_norm/{subj}_task_breathhold.nii.gz'
    elif data_type == 'hr':
        f_str = f'data/dataset_nki/physio/proc1_physio/{subj}_task_breathhold_physio_PPG_HR_NK.txt'
    elif data_type == 'rv':
        f_str = f'data/dataset_nki/physio/proc1_physio/{subj}_task_breathhold_physio_RESP_AMP_HILBERT.txt'
    elif data_type == 'csf':
        f_str = f'data/dataset_nki/physio/proc1_physio/{subj}_task_breathhold_physio_csf.txt'
    return f_str


def fp_spreng(data_type, subj):
    if data_type == 'func':
        f_str = f'data/dataset_spreng/func/proc8_bandpass/{subj}_ses-1_task-rest.nii.gz'
    elif data_type == 'hr':
        f_str = f'data/dataset_spreng/physio/proc1_physio/{subj}_ses-1_task-rest_physio_PPG_HR_NK.txt'
    elif data_type == 'rv':
        f_str = f'data/dataset_spreng/physio/proc1_physio/{subj}_ses-1_task-rest_physio_RESP_RVT_NK.txt'
    elif data_type == 'ppg_low':
        f_str = f'data/dataset_spreng/physio/proc1_physio/{subj}_ses-1_task-rest_physio_PPG_LOW_NK.txt'
    return f_str


def load_subject_list(data):
    if data == 'chang':
        subj_list = pd.read_csv(subject_list_chang)
    elif data == 'chang_bh':
        subj_list = pd.read_csv(subject_list_chang_bh)
    elif data == 'nki':
        subj_list = pd.read_csv(subject_list_nki)
    elif (data == 'hcp') | (data == 'hcp_fix'):
        subj_list = pd.read_csv(subject_list_hcp)
    elif data == 'hcp_rel':
        subj_list = pd.read_csv(subject_list_hcp_rel)
    elif data == 'hcp_wm':
        subj_list = pd.read_csv(subject_list_hcp_wm)
    elif (data == 'monash') | (data == 'monash_pet'):
        subj_list = pd.read_csv(subject_list_monash)
    elif data == 'spreng':
        subj_list = pd.read_csv(subject_list_spreng)
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



    
    