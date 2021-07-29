#!/bin/bash
mask=$1

# Make directories
mkdir -p data/bold/proc1_resample
mkdir -p data/bold/proc2_smooth_mask
mkdir -p data/bold/proc3_filter_norm


for file_path in data/bold/orig/*.nii.gz; do
    filename=$(basename $file_path)
    echo "$filename" 

    # Resample to 3mm
    flirt -in $file_path -ref $mask -out data/bold/proc1_resample/$filename -applyisoxfm 3
    # Mask
    fslmaths data/bold/proc1_resample/$filename -mul $mask data/bold/proc2_smooth_mask/$filename
    # Bandpass filter (0.01 - 0.1Hz) and norm
    python utils/norm_filter.py -n data/bold/proc2_smooth_mask/$filename -m $mask -l 0.01 -u 0.1 -o data/bold/proc3_filter_norm/$filename
done


    # fslmaths data/bold/proc1_resample/$filename -kernel gauss 2.1233226 -fmean data/bold/proc2_smooth_mask/$filename
