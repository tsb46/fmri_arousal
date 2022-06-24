import csv
import os

# Re-organize folders
os.makedirs('anat/raw', exist_ok=True)
os.makedirs('func/raw', exist_ok=True)
os.makedirs('pet/raw', exist_ok=True)


n_scans=6

anat = lambda x: f'ds002898/{x}/anat/{x}_T1w'
func = lambda x, y: f'ds002898/{x}/func/{x}_task-rest_run-{y}_bold'
pet = lambda x: f'ds002898/{x}/pet/{x}_task-rest_pet'


# Move all subject files
subj_list = []
with open('subject_list_monash.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count > 0:
            subj_list.append(row[0])
        line_count += 1

for subj in subj_list:
    # structural
    if os.path.exists(f'{anat(subj)}.json'):
        os.rename(f'{anat(subj)}.nii.gz', f'anat/raw/{subj}_T1w.nii.gz')
        os.rename(f'{anat(subj)}.json', f'anat/raw/{subj}_T1w.json')
    # functional
    for scan in range(1,n_scans+1):
        if os.path.exists(f'{func(subj,scan)}.json'):
            os.rename(f'{func(subj,scan)}.nii.gz', f'func/raw/{subj}_run_{scan}.nii.gz')
            os.rename(f'{func(subj,scan)}.json', f'func/raw/{subj}_run_{scan}.json')
    # pet
    if os.path.exists(f'{pet(subj)}.json'):
        os.rename(f'{pet(subj)}.nii.gz', f'pet/raw/{subj}_pet.nii.gz')
        os.rename(f'{pet(subj)}.json', f'pet/raw/{subj}_pet.json')

    