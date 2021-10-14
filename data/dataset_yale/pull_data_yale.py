import os
import csv
import openneuro as on

ds_num = 'ds003673'
tag_num = '1.0.1'

on_anat = lambda x: f'{x}/anat/{x}_T1w.nii.gz'
on_func = lambda x: f'{x}/func/{x}_task-rest_run-01_bold.nii.gz'
on_eye = lambda x: f'derivatives/{x}/{x}_task-rest_run-01_et.tsv'

subj_list = []
with open('subject_list_yale.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count > 0:
            subj_list.append(row[0])
        line_count += 1

# Download subject files
for subj in subj_list:
    print(subj)
    on.download(dataset=ds_num, tag=tag_num, target_dir='raw',
                include=[on_anat(subj), on_func(subj), on_eye(subj)])

# Re-organize folders
os.makedirs('anat/raw', exist_ok=True)
os.makedirs('func/raw', exist_ok=True)
os.makedirs('pupillometry', exist_ok=True)

# Move all subject files
for subj in subj_list:
    if os.path.exists(f'raw/{on_anat(subj)}'):
        os.rename(f'raw/{on_anat(subj)}', f'anat/raw/{subj}_T1w.nii.gz')
    if os.path.exists(f'raw/{on_func(subj)}'):
        os.rename(f'raw/{on_func(subj)}', f'func/raw/{subj}_task-rest_run-01_bold.nii.gz')
    if os.path.exists(f'raw/{on_eye(subj)}'):
        os.rename(f'raw/{on_eye(subj)}', f'pupillometry/{subj}_task-rest_run-01_et.tsv')







    

