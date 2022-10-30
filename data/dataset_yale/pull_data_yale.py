import os
import csv
import numpy as np
import openneuro as on
import pandas as pd 

ds_num = 'ds003673'
tag_num = '2.0.1'

on_anat = lambda x: f'{x}/anat/{x}_T1w.nii.gz'
on_func = lambda x, y: f'{x}/func/{x}_task-rest_run-0{y}_bold.nii.gz'
on_eye = lambda x, y: f'derivatives/{x}/{x}_task-rest_run-0{y}_et.tsv'

subj_list = []
with open('subject_list_yale.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count > 0:
            subj_list.append((row[0], row[1]))
        line_count += 1

eye_list = [on_eye(subj, scan) for subj, scan in subj_list]
anat_list = [on_anat(subj) for subj in set(s[0] for s in subj_list)]
func_list = [on_func(subj, scan) for subj, scan in subj_list]

full_list = anat_list + func_list + eye_list

on.download(dataset=ds_num, tag=tag_num, target_dir='raw', 
            include=full_list)

# Re-organize folders
os.makedirs('anat/raw', exist_ok=True)
os.makedirs('func/raw', exist_ok=True)
os.makedirs('physio/raw', exist_ok=True)

# Move all subject files
for subj, scan in subj_list:
    if scan == '1':
        if os.path.exists(f'raw/{on_anat(subj)}'):
            os.rename(f'raw/{on_anat(subj)}', f'anat/raw/{subj}_T1w.nii.gz')
    if os.path.exists(f'raw/{on_func(subj, scan)}'):
        os.rename(f'raw/{on_func(subj, scan)}', f'func/raw/{subj}_task-rest_run-0{scan}_bold.nii.gz')
    if os.path.exists(f'raw/{on_eye(subj, scan)}'):
        os.rename(f'raw/{on_eye(subj, scan)}', f'physio/raw/{subj}_task-rest_run-0{scan}_et.tsv')
        # Convert .tsv to .txt
        pupil_dil = pd.read_csv(f'physio/raw/{subj}_task-rest_run-0{scan}_et.tsv', 
                        delimiter='\t', header=None)
        np.savetxt(f'physio/raw/{subj}_task-rest_run-0{scan}_et.txt', pupil_dil.values)








    

