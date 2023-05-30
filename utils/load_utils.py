import json
import numpy as np
import os
import pandas as pd

physio_dict = {
    'rv': 'RESP_RVT_NK',
    'hr': 'PPG_RATE_NK',
    'rv_amp': 'RESP_RVT_NK',
    'rv_rate': 'RESP_RATE_NK',
    'resp': 'RESP_RAW',
    'ppg_low': 'PPG_LOW_NK',
    'ppg_amp': 'RESP_AMP_NK',
    'vigilance': 'vigilance_at',
    'alpha': 'fbands_Alpha',
    'delta': 'fbands_Delta',
    'infraslow': 'Infraslow',
    'global_sig': 'global_sig' 
}

physio_type = {
    'func': 'func',
    'rv': 'physio',
    'rv_amp': 'physio',
    'rv_rate': 'physio',
    'resp': 'physio',
    'hr': 'physio',
    'ppg_low': 'physio',
    'ppg_amp': 'physio',
    'vigilance': 'eeg',
    'alpha': 'eeg',
    'delta': 'eeg',
    'infraslow': 'eeg',
    'global_sig': 'global_sig',
    'pupil': 'physio'
}

def find_fps(data, physio, params, subj_n=None, scan=None):
    subj_list = load_subject_list(params['subject_list'])
    physio_fp = physio.copy()
    if data in ['chang', 'chang_bh', 'chang_cue', 'yale', 'spreng', 'natview']: 
        search_terms = subj_list[['subject', 'scan']].values.tolist()
    elif (data == 'hcp'):
        search_terms = subj_list[['subject', 'lr']].values.tolist()
    else:
        search_terms = subj_list.subject.values.tolist()

    if len(physio) > 0:
        search_physio(data, physio_fp, params['physio'])
        physio_fp.insert(0, 'func')
    else:
        physio_fp = ['func']


    if data == 'chang':
        fps = {d_type: [fp_chang(d_type, subj_scan[0],subj_scan[1], params) for subj_scan in search_terms] 
               for d_type in physio_fp}
    elif data == 'chang_bh':
        fps = {d_type: [fp_chang_bh(d_type, subj_scan[0],subj_scan[1], params) for subj_scan in search_terms] 
               for d_type in physio_fp}
    elif data == 'chang_cue':
        fps = {d_type: [fp_chang_cue(d_type, subj_scan[0],subj_scan[1], params) for subj_scan in search_terms] 
               for d_type in physio_fp}
    elif data == 'nki':
        fps = {d_type: [fp_nki(d_type, subj, params) for subj in search_terms] for d_type in physio_fp}
    elif (data == 'hcp'):
        fps = {d_type: [fp_hcp(d_type, subj_scan[0],subj_scan[1], params) for subj_scan in search_terms] 
               for d_type in physio_fp}
    elif data == 'spreng':
        fps = {d_type: [fp_spreng(d_type, subj_scan[0],subj_scan[1], params) for subj_scan in search_terms] 
               for d_type in physio_fp}
    elif data == 'yale':
        fps = {d_type: [fp_yale(d_type, subj_scan[0],subj_scan[1], params) for subj_scan in search_terms] 
               for d_type in physio_fp}
    elif data == 'natview':
        fps = {d_type: [fp_natview(d_type, subj_scan[0],subj_scan[1], params) for subj_scan in search_terms] 
               for d_type in physio_fp}
    return fps


def fp_chang(data_type, subj, scan, params):
    if scan < 10:
        scan_str = f'000{scan}'
    else:
        scan_str = f'00{scan}'

    if physio_type[data_type] == 'func':
        f_str = f'{params["func_dir"]}/sub_00{subj}-mr_{scan_str}-ecr_echo1_w_dspk_blur3mm.nii.gz' 
    elif physio_type[data_type] == 'eeg':
        f_str = f'{params["eeg_dir"]}/sub_00{subj}_mr_{scan_str}_fbands_{physio_dict[data_type]}.txt'
    elif physio_type[data_type] == 'physio':
        f_str = f'{params["physio_dir"]}/sub_00{subj}_mr_{scan_str}_physio_{physio_dict[data_type]}.txt'
    elif physio_type[data_type] == 'global_sig':
        f_str = f'{params["physio_dir"]}/sub_00{subj}_mr_{scan_str}_{physio_dict[data_type]}.txt'

    return f_str


def fp_chang_bh(data_type, subj, scan, params):
    if scan < 10:
        scan_str = f'000{scan}'
    else:
        scan_str = f'00{scan}'

    if physio_type[data_type] == 'func':
        f_str = f'{params["func_dir"]}/sub_00{subj}-mr_{scan_str}-adb_echo1_w_dspk_blur3mm.nii.gz' 
    elif physio_type[data_type] == 'eeg':
        f_str = f'{params["eeg_dir"]}/sub_00{subj}_mr_{scan_str}_fbands_{physio_dict[data_type]}.txt'
    elif physio_type[data_type] == 'physio':
        f_str = f'{params["physio_dir"]}/sub_00{subj}_mr_{scan_str}_physio_{physio_dict[data_type]}.txt'
    elif physio_type[data_type] == 'global_sig':
        f_str = f'{params["physio_dir"]}/sub_00{subj}_mr_{scan_str}_{physio_dict[data_type]}.txt'

    return f_str


def fp_chang_cue(data_type, subj, scan, params):
    if scan < 10:
        scan_str = f'000{scan}'
    else:
        scan_str = f'00{scan}'

    if physio_type[data_type] == 'func':
        f_str = f'{params["func_dir"]}/sub_00{subj}-mr_{scan_str}-ectp_echo1_w_dspk_dtr_blur3mm.nii.gz' 
    elif physio_type[data_type] == 'eeg':
        f_str = f'{params["eeg_dir"]}/sub_00{subj}_mr_{scan_str}_fbands_{physio_dict[data_type]}.txt'
    elif physio_type[data_type] == 'physio':
        f_str = f'{params["physio_dir"]}/sub_00{subj}_mr_{scan_str}_physio_{physio_dict[data_type]}.txt'
    elif physio_type[data_type] == 'global_sig':
        f_str = f'{params["physio_dir"]}/sub_00{subj}_mr_{scan_str}_{physio_dict[data_type]}.txt'

    return f_str


def fp_hcp(data_type, subj, scan, params):
    if physio_type[data_type] == 'func':
        f_str = f'{params["func_dir"]}/{subj}_{scan}1_rest.nii.gz' 
    elif physio_type[data_type] == 'physio':
        f_str = f'{params["physio_dir"]}/{subj}_physio_{physio_dict[data_type]}.txt'
    elif physio_type[data_type] == 'global_sig':
        f_str = f'{params["physio_dir"]}/{subj}_{physio_dict[data_type]}.txt'

    return f_str


def fp_natview(data_type, subj, scan, params):
    if subj < 10:
        subj_str = f'0{subj}'
    else:
        subj_str = f'{subj}'

    if physio_type[data_type] == 'func':
        f_str = f'{params["func_dir"]}/sub-{subj_str}_ses-0{scan}_func_mc.nii.gz'
    elif physio_type[data_type] == 'physio':
        f_str = f'{params["physio_dir"]}/{subj}_task-rest_run-0{scan}_et.txt'

    return f_str


def fp_nki(data_type, subj, params):
    if physio_type[data_type] == 'func':
        f_str = f'{params["func_dir"]}/{subj}_task_breathhold.nii.gz' 
    elif physio_type[data_type] == 'physio':
        f_str = f'{params["physio_dir"]}/{subj}_task_breathhold_physio_{physio_dict[data_type]}.txt'

    return f_str


def fp_spreng(data_type, subj, scan, params):
    if physio_type[data_type] == 'func':
        f_str = f'{params["func_dir"]}/{subj}_task-rest_{scan}_echo-123_bold_medn_afw.nii' 
    elif physio_type[data_type] == 'physio':
        f_str = f'{params["physio_dir"]}/{subj}_task-rest_{scan}_physio_{physio_dict[data_type]}.txt'

    return f_str


def fp_yale(data_type, subj, scan, params):
    if physio_type[data_type] == 'func':
        f_str = f'{params["func_dir"]}/{subj}_task-rest_run-0{scan}_bold.nii.gz' 
    elif physio_type[data_type] == 'physio':
        f_str = f'{params["physio_dir"]}/{subj}_task-rest_run-0{scan}_et.txt'

    return f_str


def load_subject_list(subject_list):
    subj_list = pd.read_csv(subject_list)
    return subj_list


def print_filter_info(params, load_physio, physio_list):
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
        for p in physio_list:
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

    
    