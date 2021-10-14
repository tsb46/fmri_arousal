import os
import csv

from oct2py import octave

anat = 'fmri_allrun/anat_final.20190703jp+tlrc.nii'
on_func = lambda x: f'fmri_allrun/{x}/run1/pb04.2019{x}jp.r01.blur+tlrc.nii'
egg = 'EGGGICA_prepro_ds_all.mat'

run_list = []
with open('run_list_choe.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    line_count = 0
    for row in csv_reader:
        if line_count > 0:
            run_list.append(row[0])
        line_count += 1

# Re-organize folders
os.makedirs('anat/raw', exist_ok=True)
os.makedirs('func/raw', exist_ok=True)
os.makedirs('egg/raw', exist_ok=True)

# Move all subject files
# os.rename(anat, f'anat/raw/anat_final.20190703jp+tlrc.nii')
os.rename(egg, f'egg/raw/{egg}')
for run in run_list:
    if os.path.exists(on_func(run)):
        os.rename(on_func(run), f'func/raw/pb04.2019{run}jp.r01.blur+tlrc.nii')

octave.run('reorganize_egg.m')







    

