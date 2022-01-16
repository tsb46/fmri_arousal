import csv
import os

# Re-organize folders
os.makedirs('anat/raw', exist_ok=True)
os.makedirs('func/raw', exist_ok=True)
os.makedirs('physio/raw', exist_ok=True)


anat = lambda x, y: f'mri_raw/{x}/ses-0{y}/anat/{x}_ses-0{y}_acq-mp2rage'
anat_inv = lambda x, y: f'mri_raw/{x}/ses-0{y}/anat/{x}_ses-0{y}_inv-2_mp2rage'
func = lambda x, y: f'mri_raw/{x}/ses-0{y}/func/{x}_ses-0{y}_task-rest_acq-AP_run-01_bold'
physio = lambda x, y, z: f'physio_raw/{x}/ses-0{y}/func/{x}_ses-0{y}_task-rest_acq-AP_run-01_recording-{z}_physio'


# Move all subject files
subj_list = []
with open('subject_list_lemon.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count > 0:
            subj_list.append(row[0])
        line_count += 1

for subj in subj_list:
    # anat
    if os.path.exists(f'{anat(subj,1)}_T1w.nii.gz'):
        os.rename(f'{anat(subj,1)}_T1w.nii.gz', f'anat/raw/{subj}_T1w.nii.gz')
    # anat - inv
    if os.path.exists(f'{anat_inv(subj,1)}.nii.gz'):
        os.rename(f'{anat_inv(subj,1)}.nii.gz', f'anat/raw/{subj}_inv2.nii.gz')
    # functional
    if os.path.exists(f'{func(subj,1)}.json'):
        os.rename(f'{func(subj,1)}.nii.gz', f'func/raw/{subj}_task_rest.nii.gz')
        os.rename(f'{func(subj,1)}.json', f'func/raw/{subj}_task_rest.json')
     # physio
    for phys_label in ['bpp', 'ppg', 'ecg', 'resp']:
        if os.path.exists(f'{physio(subj,1,phys_label)}.tsv.gz'):
            os.rename(f'{physio(subj,1,phys_label)}.tsv.gz', f'physio/raw/{subj}_task_rest_physio_{phys_label}.tsv.gz')
            os.rename(f'{physio(subj,1,phys_label)}.json', f'physio/raw/{subj}_task_rest_physio_{phys_label}.json')
            if phys_label == 'bpp':
                if os.path.exists(f'{physio(subj,1,"dbp")}.png'):
                    os.rename(f'{physio(subj,1,"dbp")}.png', f'physio/raw/{subj}_task_rest_physio_dbp.png')
                    os.rename(f'{physio(subj,1,"sbp")}.png', f'physio/raw/{subj}_task_rest_physio_sbp.png')
            else:
                if os.path.exists(f'{physio(subj,1,phys_label)}.png'):
                    os.rename(f'{physio(subj,1,phys_label)}.png', f'physio/raw/{subj}_task_rest_physio_{phys_label}.png')




