import os
import csv
import openneuro as on

ds_num = 'ds003768'
tag_num = '1.0.2'

on_anat = lambda x: f'{x}/anat/{x}_T1w.nii.gz'
on_func = lambda x: f'{x}/func/{x}_task-sleep_run-1_bold.nii.gz'
on_eeg = lambda x: f'{x}/eeg/{x}_task-sleep_run-1_eeg.*'
on_sleep_stage = lambda x: f'sourcedata/{x}-sleep-stage.tsv'

subj_list = []
with open('subject_list_gu.csv') as csv_file:
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
                include=[on_anat(subj), on_func(subj), on_eeg(subj), on_sleep_stage(subj)])

# Re-organize folders
os.makedirs('anat/raw', exist_ok=True)
os.makedirs('func/raw', exist_ok=True)
os.makedirs('eeg/raw', exist_ok=True)
os.makedirs('sleep_stage', exist_ok=True)

# Move all subject files
for subj in subj_list:
    if os.path.exists(f'raw/{on_anat(subj)}'):
        os.rename(f'raw/{on_anat(subj)}', f'anat/raw/{subj}_T1w.nii.gz')
    if os.path.exists(f'raw/{on_func(subj)}'):
        os.rename(f'raw/{on_func(subj)}', f'func/raw/{subj}_task-sleep_run-1_bold.nii.gz')
    if os.path.exists(f'raw/{subj}/eeg/{subj}_task-sleep_run-1_eeg.eeg'):
        os.rename(f'raw/{subj}/eeg/{subj}_task-sleep_run-1_eeg.eeg', 
                  f'eeg/raw/{subj}_task-sleep_run-1_eeg.eeg')
    if os.path.exists(f'raw/{subj}/eeg/{subj}_task-sleep_run-1_eeg.vhdr'):
        os.rename(f'raw/{subj}/eeg/{subj}_task-sleep_run-1_eeg.vhdr', 
                  f'eeg/raw/{subj}_task-sleep_run-1_eeg.vhdr')
    if os.path.exists(f'raw/{subj}/eeg/{subj}_task-sleep_run-1_eeg.vmrk'):
        os.rename(f'raw/{subj}/eeg/{subj}_task-sleep_run-1_eeg.vmrk', 
                  f'eeg/raw/{subj}_task-sleep_run-1_eeg.vmrk')
    if os.path.exists(f'raw/{on_sleep_stage(subj)}'):
        os.rename(f'raw/{on_sleep_stage(subj)}', f'sleep_stage/{subj}-sleep-stage.tsv')







    

