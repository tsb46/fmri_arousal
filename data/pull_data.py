import argparse
import boto3
import botocore
import csv
import openneuro as on
import os
import pandas as pd
import shutil

from botocore import UNSIGNED
from botocore.client import Config


# Datasets for download
dataset = ['hcp', 'natview', 'nki', 'nki_rest', 'spreng', 'yale']


def download_hcp(subjects):
    # Code borrowed and modified from:
    # https://github.com/jokedurnez/HCP_download
    # download HCP data from AWS via boto3

    # Create directories
    os.makedirs('dataset_hcp/func/raw', exist_ok=True)
    os.makedirs('dataset_hcp/physio/raw', exist_ok=True)

    # Set up S3 bucket
    boto3.setup_default_session(profile_name='hcp')
    s3 = boto3.resource('s3')
    bucket = s3.Bucket('hcp-openaccess')

    # Load subject list and iterate through subjects
    for s, lr in zip(subjects.subject, subjects.lr):
        # Set up file path strings
        # define base directories
        s_dir = f'HCP_1200/{s}/MNINonLinear/Results/rfMRI_REST1_{lr}'

        # Pull hcp-fix cleaned cifti file
        func_fp = f'{s_dir}/rfMRI_REST1_{lr}_hp2000_clean.nii.gz'
        func_out = f'dataset_hcp/func/raw/{s}_{lr}_hp2000_clean.nii.gz'
        bucket.download_file(func_fp, func_out)

        # Pull physio .txt file
        phys_fp = f'{s_dir}/rfMRI_REST1_{lr}_Physio_log.txt'
        phys_out = f'dataset_hcp/physio/raw/{s}_{lr}_physio.txt'
        bucket.download_file(phys_fp, phys_out)


def download_nki(subjects, aws_links):
    # Pull NKI Rockland breath-hold sessions via AWS
    # Re-organize folders
    os.makedirs('dataset_nki/anat/raw', exist_ok=True)
    os.makedirs('dataset_nki/func/raw', exist_ok=True)
    os.makedirs('dataset_nki/physio/raw', exist_ok=True)
    os.makedirs('dataset_nki/events', exist_ok=True)

    # Init variables
    s3_bucket_name = 'fcp-indi'
    s3_prefix = 'data/Projects/RocklandSample/RawDataBIDSLatest'

    # Set up S3 bucket
    s3 = boto3.resource('s3')
    s3.meta.client.meta.events.register('choose-signer.s3.*', disable_signing)
    bucket = s3.Bucket(s3_bucket_name)

    # Set up templates
    anat_template = 'sub-{0}/ses-{1}/anat/sub-{0}_ses-{1}_T1w.nii.gz'
    func_template = 'sub-{0}/ses-{1}/func/sub-{0}_ses-{1}_task-BREATHHOLD_acq-1400_bold.nii.gz'
    physio_template = 'sub-{0}/ses-{1}/func/sub-{0}_ses-{1}_task-BREATHHOLD_acq-1400_physio.tsv.gz'
    physio_json_template = 'sub-{0}/ses-{1}/func/sub-{0}_ses-{1}_task-BREATHHOLD_acq-1400_physio.json'
    event_template = 'sub-{0}/ses-{1}/func/sub-{0}_ses-{1}_task-BREATHHOLD_acq-1400_events.tsv'

    # iterate through subjects and download data
    for s, ses in zip(subjects.subject, subjects.session):
        # Pull anatomical file
        anat_fp = anat_template.format(s, ses)
        anat_in = f"{s3_prefix}/{anat_fp}"
        anat_out = f"dataset_nki/anat/raw/sub-{s}_T1w.nii.gz"
        bucket.download_file(anat_in, anat_out)
        # Pull functional file
        func_fp = func_template.format(s, ses)
        func_in = f"{s3_prefix}/{func_fp}"
        func_out = f"dataset_nki/func/raw/{os.path.basename(func_fp)}"
        bucket.download_file(func_in, func_out)
        # Pull physio file
        physio_fp = physio_template.format(s, ses)
        physio_in = f"{s3_prefix}/{physio_fp}"
        physio_out = f"dataset_nki/physio/raw/{os.path.basename(physio_fp)}"
        bucket.download_file(physio_in, physio_out)
        # Pull physio json file
        physio_json_fp = physio_json_template.format(s, ses)
        physio_json_in = f"{s3_prefix}/{physio_json_fp}"
        physio_json_out = f"dataset_nki/physio/raw/{os.path.basename(physio_json_fp)}"
        bucket.download_file(physio_json_in, physio_json_out)
        # Pull event file
        event_fp = event_template.format(s, ses)
        event_in = f"{s3_prefix}/{event_fp}"
        event_out = f"dataset_nki/events/{os.path.basename(event_fp)}"
        bucket.download_file(event_in, event_out)



def download_nki_rest(subjects):
    # Pull NKI Rockland resting-state sessions via AWS
    # Create directories
    os.makedirs('dataset_nki_rest/anat/raw', exist_ok=True)
    os.makedirs('dataset_nki_rest/func/raw', exist_ok=True)
    os.makedirs('dataset_nki_rest/physio/raw', exist_ok=True)

    # Init variables
    s3_bucket_name = 'fcp-indi'
    s3_prefix = 'data/Projects/RocklandSample/RawDataBIDSLatest'

    # Set up S3 bucket
    s3 = boto3.resource('s3')
    s3.meta.client.meta.events.register('choose-signer.s3.*', disable_signing)
    bucket = s3.Bucket(s3_bucket_name)

    # Set up templates
    anat_template = 'sub-{0}/ses-{1}/anat/sub-{0}_ses-{1}_T1w.nii.gz'
    func_template = 'sub-{0}/ses-{1}/func/sub-{0}_ses-{1}_task-rest_acq-1400_bold.nii.gz'
    physio_template = 'sub-{0}/ses-{1}/func/sub-{0}_ses-{1}_task-rest_acq-1400_physio.tsv.gz'
    physio_json_template = 'sub-{0}/ses-{1}/func/sub-{0}_ses-{1}_task-rest_acq-1400_physio.json'


    # iterate through subjects and download data
    for s, ses in zip(subjects.subject, subjects.session):
        # Pull anatomical file
        anat_fp = anat_template.format(s, ses)
        anat_in = f"{s3_prefix}/{anat_fp}"
        anat_out = f"dataset_nki_rest/anat/raw/sub-{s}_T1w.nii.gz"
        bucket.download_file(anat_in, anat_out)
        # Pull functional file
        func_fp = func_template.format(s, ses)
        func_in = f"{s3_prefix}/{func_fp}"
        func_out = f"dataset_nki_rest/func/raw/{os.path.basename(func_fp)}"
        bucket.download_file(func_in, func_out)
        # Pull physio file
        physio_fp = physio_template.format(s, ses)
        physio_in = f"{s3_prefix}/{physio_fp}"
        physio_out = f"dataset_nki_rest/physio/raw/{os.path.basename(physio_fp)}"
        bucket.download_file(physio_in, physio_out)
        # Pull physio json file
        physio_json_fp = physio_json_template.format(s, ses)
        physio_json_in = f"{s3_prefix}/{physio_json_fp}"
        physio_json_out = f"dataset_nki_rest/physio/raw/{os.path.basename(physio_json_fp)}"
        bucket.download_file(physio_json_in, physio_json_out)


def download_natview(subjects):
    # download natview dataset from AWS via boto3
    # Create directories
    os.makedirs('dataset_natview/anat/raw', exist_ok=True)
    os.makedirs('dataset_natview/func/raw', exist_ok=True)
    os.makedirs('dataset_natview/physio/raw', exist_ok=True)
    os.makedirs('dataset_natview/eeg/raw', exist_ok=True)

    # Init variables
    s3_bucket_name = 'fcp-indi'
    s3_prefix_proc = 'data/Projects/NATVIEW_EEGFMRI/preproc_data'
    s3_prefix_raw = 'data/Projects/NATVIEW_EEGFMRI/raw_data'

    # Set up S3 bucket
    s3 = boto3.resource('s3')
    bucket = s3.Bucket(s3_bucket_name)
    s3_client = boto3.client('s3', config=Config(signature_version=UNSIGNED))

    # Set up templates
    anat_template = 'anat/T1w1_denoise.nii.gz'
    func_template = 'func/sub-{0}_ses-{1}_task-rest_bold.nii.gz'
    resp_json_template = 'func/sub-{0}_ses-{1}_task-rest_recording-respiratory_physio.json'
    resp_template = 'func/sub-{0}_ses-{1}_task-rest_recording-respiratory_physio.tsv.gz'
    eye_template = 'eeg/sub-{0}_ses-{1}_task-rest_recording-eyetracking_physio.tsv.gz'
    eye_json_template = 'eeg/sub-{0}_ses-{1}_task-rest_recording-eyetracking_physio.json'
    eeg_template = 'eeg/sub-{0}_ses-{1}_task-rest_eeg.set'
    eeg_channel_template = 'eeg/sub-{0}_ses-{1}_task-rest_channels.tsv'
    eeg_json_template = 'eeg/sub-{0}_ses-{1}_task-rest_eeg.json'


    anat_templates = [anat_template]
    func_templates = [func_template]
    resp_templates = [resp_json_template, resp_template]
    eye_templates = [eye_template, eye_json_template]
    eeg_templates = [eeg_template, eeg_channel_template, eeg_json_template]

    # Load subject list and iterate through subjects
    for s, ses in zip(subjects.subject, subjects.scan):
        # subject session integers to strings
        if s < 10:
            subj = f'0{s}'
        else:
            subj = f'{s}'
        ses_str = f'0{ses}'

        # Pull anatomical files (and warp)
        for template in anat_templates:
            temp_base = os.path.basename(template)
            anat_fp = f'{s3_prefix_proc}/sub-{subj}/ses-{ses_str}/{template}'
            anat_out = f'dataset_natview/anat/raw/sub-{subj}_{temp_base}'
            with open(anat_out, 'wb') as f:
                s3_client.download_fileobj(s3_bucket_name, anat_fp, f)

        # Pull functional files (and coregistration affine)
        for template in func_templates:
            temp_base = os.path.basename(template.format(subj, ses_str))
            func_fp = f'{s3_prefix_raw}/sub-{subj}/ses-{ses_str}/{template.format(subj, ses_str)}'
            func_out = f'dataset_natview/func/raw/{temp_base}'
            with open(func_out, 'wb') as f:
                s3_client.download_fileobj(s3_bucket_name, func_fp, f)

        # Pull respiration files
        for template in resp_templates:
            temp_base = os.path.basename(template.format(subj, ses_str))
            resp_fp = f'{s3_prefix_raw}/sub-{subj}/ses-{ses_str}/{template.format(subj, ses_str)}'
            resp_out = f'dataset_natview/physio/raw/{temp_base}'
            with open(resp_out, 'wb') as f:
                s3_client.download_fileobj(s3_bucket_name, resp_fp, f)

        # # Pull eye-tracking files
        for template in eye_templates:
            temp_base = os.path.basename(template.format(subj, ses_str))
            eye_fp = f'{s3_prefix_proc}/sub-{subj}/ses-{ses_str}/{template.format(subj, ses_str)}'
            eye_out = f'dataset_natview/physio/raw/{temp_base}'
            with open(eye_out, 'wb') as f:
                s3_client.download_fileobj(s3_bucket_name, eye_fp, f)

        # Pull eeg files
        for template in eeg_templates:
            temp_base = os.path.basename(template.format(subj, ses_str))
            eeg_fp = f'{s3_prefix_raw}/sub-{subj}/ses-{ses_str}/{template.format(subj, ses_str)}'
            eeg_out = f'dataset_natview/eeg/raw/{temp_base}'
            with open(eeg_out, 'wb') as f:
                s3_client.download_fileobj(s3_bucket_name, eeg_fp, f)


def download_spreng(subjects):
    # download spreng neurocognitive aging dataset via openneuro python package
    # specify ds number and version #
    ds_num = 'ds003592'
    tag_num = '1.0.13'

    on_anat = lambda x: f'{x}/ses-1/anat/{x}_ses-1_T1w.nii.gz'
    on_func = lambda x, y, z: f'{x}/{y}/func/{x}_{y}_task-rest_echo-{z}_bold.nii.gz'
    on_physio = lambda x, y: f'{x}/{y}/func/{x}_{y}_task-rest_physio.tsv.gz'
    on_physio_json = lambda x, y: f'{x}/{y}/func/{x}_{y}_task-rest_physio.json'
    # get subj and scan ids
    subj_list = subjects.subject
    subj_list_unq = list(dict.fromkeys(subj_list))
    scan_list = subjects.scan
    # get file paths
    anat_list = [on_anat(subj) for subj in subj_list_unq]
    func_list = [on_func(subj, scan, echo+1) 
                 for subj, scan in zip(subj_list, scan_list) 
                 for echo in range(3)]
    physio_list = [on_physio(subj, scan) 
                   for subj, scan in zip(subj_list, scan_list)]
    physio_json_list = [on_physio_json(subj, scan) 
                        for subj, scan in zip(subj_list, scan_list)]

    full_list = anat_list + func_list + physio_list + physio_json_list
    # download with openneuro-py
    on.download(dataset=ds_num, tag=tag_num, target_dir='dataset_spreng/aws_dir', 
                include=full_list)

    # Re-organize folders
    os.makedirs('dataset_spreng/anat/raw', exist_ok=True)
    os.makedirs('dataset_spreng/func/raw', exist_ok=True)
    os.makedirs('dataset_spreng/physio/raw', exist_ok=True)

    # Move all subject files
    for subj, scan in zip(subj_list, scan_list):
        if os.path.exists(f'dataset_spreng/aws_dir/{on_anat(subj)}'):
            os.rename(f'dataset_spreng/aws_dir/{on_anat(subj)}', 
                      f'dataset_spreng/anat/raw/{subj}_T1w.nii.gz')
        for echo in range(3):
            if os.path.exists(f'dataset_spreng/aws_dir/{on_func(subj, scan, echo+1)}'):
                os.rename(f'dataset_spreng/aws_dir/{on_func(subj, scan, echo+1)}', 
                          f'dataset_spreng/func/raw/{subj}_{scan}_task-rest_echo-{echo+1}_bold.nii.gz')
        if os.path.exists(f'dataset_spreng/aws_dir/{on_physio(subj, scan)}'):
            os.rename(f'dataset_spreng/aws_dir/{on_physio(subj, scan)}', 
                      f'dataset_spreng/physio/raw/{subj}_{scan}_task-rest_physio.tsv.gz')
        if os.path.exists(f'dataset_spreng/aws_dir/{on_physio_json(subj, scan)}'):
            os.rename(f'dataset_spreng/aws_dir/{on_physio_json(subj, scan)}', 
                      f'dataset_spreng/physio/raw/{subj}_{scan}_task-rest_physio.json')
    # Delete temporary aws folder
    shutil.rmtree('dataset_spreng/aws_dir')


def download_yale(subjects):
    # download yale pupillometry dataset via openneuro python package
    # specify ds number and version #
    ds_num = 'ds003673'
    tag_num = '2.0.1'

    on_anat = lambda x: f'{x}/anat/{x}_T1w.nii.gz'
    on_func = lambda x, y: f'{x}/func/{x}_task-rest_run-0{y}_bold.nii.gz'
    on_eye = lambda x, y: f'derivatives/{x}/{x}_task-rest_run-0{y}_et.tsv'

    subj_list = subjects.subject
    scan_list = subjects.scan

    eye_list = [on_eye(subj, scan) for subj, scan in zip(subj_list, scan_list)]
    anat_list = [on_anat(subj) for subj in subj_list.unique()]
    func_list = [on_func(subj, scan) for subj, scan in zip(subj_list, scan_list)]

    full_list = anat_list + func_list + eye_list

    on.download(dataset=ds_num, tag=tag_num, target_dir='dataset_yale/aws_dir', 
                include=full_list)

    # Re-organize folders
    os.makedirs('dataset_yale/anat/raw', exist_ok=True)
    os.makedirs('dataset_yale/func/raw', exist_ok=True)
    os.makedirs('dataset_yale/physio/raw', exist_ok=True)

    # Move all subject files
    for subj, scan in zip(subj_list, scan_list):
        if scan == 1:
            if os.path.exists(f'dataset_yale/aws_dir/{on_anat(subj)}'):
                os.rename(f'dataset_yale/aws_dir/{on_anat(subj)}', 
                          f'dataset_yale/anat/raw/{subj}_T1w.nii.gz')
        if os.path.exists(f'dataset_yale/aws_dir/{on_func(subj, scan)}'):
            os.rename(f'dataset_yale/aws_dir/{on_func(subj, scan)}', 
                      f'dataset_yale/func/raw/{subj}_task-rest_run-0{scan}_bold.nii.gz')
        if os.path.exists(f'dataset_yale/aws_dir/{on_eye(subj, scan)}'):
            os.rename(f'dataset_yale/aws_dir/{on_eye(subj, scan)}', 
                      f'dataset_yale/physio/raw/{subj}_task-rest_run-0{scan}_et.tsv')
    # Delete temporary aws folder
    shutil.rmtree('dataset_yale/aws_dir')


def pull_data(dataset):
    print(f'downloading {dataset} dataset')
    if dataset == 'hcp':
        subjects = pd.read_csv('dataset_hcp/subject_list_hcp.csv')
        download_hcp(subjects)
    elif dataset == 'natview':
        subjects = pd.read_csv('dataset_natview/subject_list_natview.csv')
        download_natview(subjects)
    elif dataset == 'nki':
        subjects = pd.read_csv('dataset_nki/subject_list_nki.csv')
        aws_links = 'dataset_nki/aws_links_sample.csv'
        download_nki(subjects, aws_links)
    elif dataset == 'nki_rest':
        subjects = pd.read_csv('dataset_nki_rest/subject_list_nki_rest.csv')
        download_nki_rest(subjects)
    elif dataset == 'spreng':
        subjects = pd.read_csv('dataset_spreng/subject_list_spreng.csv')
        download_spreng(subjects)
    elif dataset == 'yale': 
        subjects = pd.read_csv('dataset_yale/subject_list_yale.csv')
        download_yale(subjects)


if __name__ == '__main__':
    """Run main analysis"""
    parser = argparse.ArgumentParser(description='Pull datasets from AWS')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to pull from AWS - '
                        'to run all datasets use the arg "all"',
                        choices=['all', 'hcp', 'natview', 'nki',
                                 'nki_rest', 'spreng', 'yale'], 
                        required=True,
                        type=str)
    args_dict = vars(parser.parse_args())
    if args_dict['dataset'] == 'all':
        for d in dataset:
            pull_data(d)
    else:
        pull_data(args_dict['dataset'])
