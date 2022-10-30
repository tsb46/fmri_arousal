import csv
import os

# Re-organize folders
os.makedirs('anat/raw', exist_ok=True)
os.makedirs('func/raw', exist_ok=True)
os.makedirs('physio/raw', exist_ok=True)
os.makedirs('events', exist_ok=True)

anat = lambda x, y: f'raw/sub-{x}/ses-{y}/anat/sub-{x}_ses-{y}_T1w'
func = lambda x, y: f'raw/sub-{x}/ses-{y}/func/sub-{x}_ses-{y}_task-BREATHHOLD_acq-1400_bold'
events = lambda x, y: f'raw/sub-{x}/ses-{y}/func/sub-{x}_ses-{y}_task-BREATHHOLD_acq-1400_events'
physio = lambda x, y: f'raw/sub-{x}/ses-{y}/func/sub-{x}_ses-{y}_task-BREATHHOLD_acq-1400_physio'

# Move all subject files
subj_list = []
with open('aws_links_sample.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count > 0:
            subj_list.append(row[0])
        line_count += 1

for subj in subj_list:
    # structural
    if os.path.exists(f'{anat(subj,"BAS1")}.json'):
        os.rename(f'{anat(subj,"BAS1")}.nii.gz', f'anat/raw/{subj}_T1w.nii.gz')
        os.rename(f'{anat(subj,"BAS1")}.json', f'anat/raw/{subj}_T1w.json')
    else:
        if os.path.exists(f'{anat(subj,"FLU2")}.json'):
            os.rename(f'{anat(subj,"FLU2")}.nii.gz', f'anat/raw/{subj}_T1w.nii.gz')
            os.rename(f'{anat(subj,"FLU2")}.json', f'anat/raw/{subj}_T1w.json')
    # functional
    if os.path.exists(f'{func(subj,"BAS1")}.json'):
        os.rename(f'{func(subj,"BAS1")}.nii.gz', f'func/raw/{subj}_task_breathhold.nii.gz')
        os.rename(f'{func(subj,"BAS1")}.json', f'func/raw/{subj}_task_breathold.json')
    else:
        if os.path.exists(f'{func(subj,"FLU2")}.json'):
            os.rename(f'{func(subj,"FLU2")}.nii.gz', f'func/raw/{subj}_task_breathhold.nii.gz')
            os.rename(f'{func(subj,"FLU2")}.json', f'func/raw/{subj}_task_breathhold.json')
    # events
    if os.path.exists(f'{events(subj,"BAS1")}.json'):
        os.rename(f'{events(subj,"BAS1")}.tsv', f'events/{subj}_task_breathhold_events.tsv')
        os.rename(f'{events(subj,"BAS1")}.json', f'events/{subj}_task_breathhold_events.json')
    else:
        if os.path.exists(f'{events(subj,"FLU2")}.json'):
            os.rename(f'{events(subj,"FLU2")}.tsv', f'events/{subj}_task_breathhold_events.tsv')
            os.rename(f'{events(subj,"FLU2")}.json', f'events/{subj}_task_breathhold_events.json')
     # physio
    if os.path.exists(f'{physio(subj,"BAS1")}.json'):
        os.rename(f'{physio(subj,"BAS1")}.tsv.gz', f'physio/raw/{subj}_task_breathhold_physio.tsv.gz')
        os.rename(f'{physio(subj,"BAS1")}.json', f'physio/raw/{subj}_task_breathhold_physio.json')
    else:
        if os.path.exists(f'{physio(subj,"FLU2")}.json'):
            os.rename(f'{physio(subj,"FLU2")}.tsv.gz', f'physio/raw/{subj}_task_breathhold_physio.tsv.gz')
            os.rename(f'{physio(subj,"FLU2")}.json', f'physio/raw/{subj}_task_breathhold_physio.json')