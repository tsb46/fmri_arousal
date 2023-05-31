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
# import downloader script provided NKI
from dataset_nki.download_rockland_raw_bids_ver2 import collect_and_download as \
	collect_and_download_nki


# Datasets for download
dataset = ['nki', 'hcp', 'yale', 'natview']


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
	    print(s)
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
	# Download subjects from NKI dataset via AWS
	# Use downloader script provided by NKI
	collect_and_download_nki('dataset_nki/aws_dir', aws_links)

	# Re-organize folders
	os.makedirs('dataset_nki/anat/raw', exist_ok=True)
	os.makedirs('dataset_nki/func/raw', exist_ok=True)
	os.makedirs('dataset_nki/physio/raw', exist_ok=True)
	os.makedirs('dataset_nki/events', exist_ok=True)

	anat = lambda x, y: f'dataset_nki/aws_dir/sub-{x}/ses-{y}/anat/sub-{x}_ses-{y}_T1w'
	func = lambda x, y: f'dataset_nki/aws_dir/sub-{x}/ses-{y}/func/sub-{x}_ses-{y}_task-BREATHHOLD_acq-1400_bold'
	events = lambda x, y: f'dataset_nki/aws_dir/sub-{x}/ses-{y}/func/sub-{x}_ses-{y}_task-BREATHHOLD_acq-1400_events'
	physio = lambda x, y: f'dataset_nki/aws_dir/sub-{x}/ses-{y}/func/sub-{x}_ses-{y}_task-BREATHHOLD_acq-1400_physio'

	# loop through subjects and move to appopriate location
	for subj in subjects.subject:
	    # structural
	    if os.path.exists(f'{anat(subj,"BAS1")}.json'):
	        os.rename(f'{anat(subj,"BAS1")}.nii.gz', 
	                  f'dataset_nki/anat/raw/{subj}_T1w.nii.gz')
	        os.rename(f'{anat(subj,"BAS1")}.json', 
	                  f'dataset_nki/anat/raw/{subj}_T1w.json')
	    else:
	        if os.path.exists(f'{anat(subj,"FLU2")}.json'):
	            os.rename(f'{anat(subj,"FLU2")}.nii.gz', 
	                      f'dataset_nki/anat/raw/{subj}_T1w.nii.gz')
	            os.rename(f'{anat(subj,"FLU2")}.json', 
	                      f'dataset_nki/anat/raw/{subj}_T1w.json')
	    # functional
	    if os.path.exists(f'{func(subj,"BAS1")}.json'):
	        os.rename(f'{func(subj,"BAS1")}.nii.gz', 
	                  f'dataset_nki/func/raw/{subj}_task_breathhold.nii.gz')
	        os.rename(f'{func(subj,"BAS1")}.json', 
	                  f'dataset_nki/func/raw/{subj}_task_breathold.json')
	    else:
	        if os.path.exists(f'{func(subj,"FLU2")}.json'):
	            os.rename(f'{func(subj,"FLU2")}.nii.gz', 
	                      f'dataset_nki/func/raw/{subj}_task_breathhold.nii.gz')
	            os.rename(f'{func(subj,"FLU2")}.json', 
	                      f'dataset_nki/func/raw/{subj}_task_breathhold.json')
	    # events
	    if os.path.exists(f'{events(subj,"BAS1")}.json'):
	        os.rename(f'{events(subj,"BAS1")}.tsv', 
	                  f'dataset_nki/events/{subj}_task_breathhold_events.tsv')
	        os.rename(f'{events(subj,"BAS1")}.json', 
	                  f'dataset_nki/events/{subj}_task_breathhold_events.json')
	    else:
	        if os.path.exists(f'{events(subj,"FLU2")}.json'):
	            os.rename(f'{events(subj,"FLU2")}.tsv', 
	                      f'dataset_nki/events/{subj}_task_breathhold_events.tsv')
	            os.rename(f'{events(subj,"FLU2")}.json', 
	                      f'dataset_nki/events/{subj}_task_breathhold_events.json')
	     # physio
	    if os.path.exists(f'{physio(subj,"BAS1")}.json'):
	        os.rename(f'{physio(subj,"BAS1")}.tsv.gz', 
	                  f'dataset_nki/physio/raw/{subj}_task_breathhold_physio.tsv.gz')
	        os.rename(f'{physio(subj,"BAS1")}.json', 
	                  f'dataset_nki/physio/raw/{subj}_task_breathhold_physio.json')
	    else:
	        if os.path.exists(f'{physio(subj,"FLU2")}.json'):
	            os.rename(f'{physio(subj,"FLU2")}.tsv.gz', 
	                      f'dataset_nki/physio/raw/{subj}_task_breathhold_physio.tsv.gz')
	            os.rename(f'{physio(subj,"FLU2")}.json', 
	                      f'dataset_nki/physio/raw/{subj}_task_breathhold_physio.json')
	# Delete temporary aws folder
	shutil.rmtree('dataset_nki/aws_dir')


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
	anat2standard_template='anat/reg/highres2standard.mat'
	anatwarp_template='anat/reg/highres2standard_warp.nii.gz'
	func2highres_template='func/sub-{0}_ses-{1}_task-rest_bold/func_reg/example_func2highres.mat'
	func_template = 'func/sub-{0}_ses-{1}_task-rest_bold/func_minimal/func_mc.nii.gz'
	resp_json_template = 'func/sub-{0}_ses-{1}_task-rest_recording-respiratory_physio.json'
	resp_template = 'func/sub-{0}_ses-{1}_task-rest_recording-respiratory_physio.tsv.gz'
	eye_template = 'eeg/sub-{0}_ses-{1}_task-rest_recording-eyetracking_physio.tsv.gz'
	eye_json_template = 'eeg/sub-{0}_ses-{1}_task-rest_recording-eyetracking_physio.json'
	eeg_template = 'eeg/sub-{0}_ses-{1}_task-rest_eeg.set'
	eeg_channel_template = 'eeg/sub-{0}_ses-{1}_task-rest_channels.tsv'
	eeg_json_template = 'eeg/sub-{0}_ses-{1}_task-rest_eeg.json'


	anat_templates = [anat_template, anat2standard_template, anatwarp_template]
	func_templates = [func2highres_template, func_template]
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
	        anat_out = f'dataset_natview/anat/raw/sub-{subj}_ses-{ses_str}_{temp_base}'
	        with open(anat_out, 'wb') as f:
	            s3_client.download_fileobj(s3_bucket_name, anat_fp, f)

	    # Pull functional files (and coregistration affine)
	    for template in func_templates:
	        temp_base = os.path.basename(template.format(subj, ses_str))
	        func_fp = f'{s3_prefix_proc}/sub-{subj}/ses-{ses_str}/{template.format(subj, ses_str)}'
	        func_out = f'dataset_natview/func/raw/sub-{subj}_ses-{ses_str}_{temp_base}'
	        with open(func_out, 'wb') as f:
	            s3_client.download_fileobj(s3_bucket_name, func_fp, f)

	    # Pull functional files (and coregistration affine)
	    for template in func_templates:
	        temp_base = os.path.basename(template.format(subj, ses_str))
	        func_fp = f'{s3_prefix_proc}/sub-{subj}/ses-{ses_str}/{template.format(subj, ses_str)}'
	        func_out = f'dataset_natview/func/raw/sub-{subj}_ses-{ses_str}_{temp_base}'
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
	        eeg_fp = f'{s3_prefix_proc}/sub-{subj}/ses-{ses_str}/{template.format(subj, ses_str)}'
	        eeg_out = f'dataset_natview/eeg/raw/{temp_base}'
	        with open(eeg_out, 'wb') as f:
	            s3_client.download_fileobj(s3_bucket_name, eeg_fp, f)


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
	    if scan == '1':
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
	elif dataset == 'nki':
		subjects = pd.read_csv('dataset_nki/subject_list_nki.csv')
		aws_links = 'dataset_nki/aws_links_sample.csv'
		download_nki(subjects, aws_links)
	elif dataset == 'yale': 
		subjects = pd.read_csv('dataset_yale/subject_list_yale.csv')
		download_yale(subjects)
	elif dataset == 'natview':
		subjects = pd.read_csv('dataset_natview/subject_list_natview.csv')
		download_natview(subjects)





if __name__ == '__main__':
    """Run main analysis"""
    parser = argparse.ArgumentParser(description='Pull datasets from AWS')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to pull from AWS - '
                        'to run all datasets use the arg "all"',
                        choices=['all', 'nki', 'hcp', 'yale', 'natview'], 
                        required=True,
                        type=str)
    args_dict = vars(parser.parse_args())
    if args_dict['dataset'] == 'all':
        for d in dataset:
            pull_data(d)
    else:
        pull_data(args_dict['dataset'])
