#!/bin/bash
dataset=$1

mask="masks/MNI152_T1_3mm_brain_mask.nii.gz"

# Chang - EEG, FMRI, physio data
if [ "$dataset" == "chang" ]; then
    # Make directories
    mkdir -p data/dataset_chang/func/proc1_resample
    mkdir -p data/dataset_chang/func/proc2_smooth_mask
    mkdir -p data/dataset_chang/func/proc3_filter_norm

    for file_path in data/dataset_chang/func/raw/*.nii.gz; do
        filename=$(basename $file_path)
        echo "$filename" 

        # Resample to 3mm
        flirt -in $file_path -ref $mask -out data/dataset_chang/func/proc1_resample/$filename -applyisoxfm 3
        # Mask
        fslmaths data/dataset_chang/func/proc1_resample/$filename -mul $mask data/dataset_chang/func/proc2_smooth_mask/$filename
        # Lowpass filter (0.1Hz) and norm
        python utils/norm_filter.py -f data/dataset_chang/func/proc2_smooth_mask/$filename -m $mask \
         -ch 0.1 -t 2.1 -o data/dataset_chang/func/proc3_filter_norm/$filename
    done
fi

# Gu - Sleep, EEG, FMRI
if [ "$dataset" == "gu" ]; then
    # mkdir -p data/dataset_gu/anat/proc1_bet
    # mkdir -p data/dataset_gu/anat/proc2_affine
    # mkdir -p data/dataset_gu/anat/proc3_fnirt

    # echo "Structural preprocessing..."
    # # Structural preprocessing
    # for file_path in data/dataset_gu/anat/raw/*.nii.gz; do
    #     filename=$(basename $file_path)
    #     echo "$filename"
    #     # Brain Extraction
    #     bet $file_path data/dataset_gu/anat/proc1_bet/$filename -o -m
    #     # Affine registration
    #     flirt -in data/dataset_gu/anat/proc1_bet/$filename -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz \
    #     -out data/dataset_gu/anat/proc2_affine/$filename -omat data/dataset_gu/anat/proc2_affine/$filename.mat
    #     # Nonlinear transformation
    #     fnirt --ref=$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz --in=$file_path \
    #     --iout=data/dataset_gu/anat/proc3_fnirt/$filename --cout=data/dataset_gu/anat/proc3_fnirt/$filename.mat \
    #     --aff=data/dataset_gu/anat/proc2_affine/$filename.mat --config=T1_2_MNI152_2mm --warpres=6,6,6
    # done

    mkdir -p data/dataset_gu/func/proc1_mcflirt
    mkdir -p data/dataset_gu/func/proc2_slicetime
    mkdir -p data/dataset_gu/func/proc3A_firstvolume
    mkdir -p data/dataset_gu/func/proc3B_func2struct
    mkdir -p data/dataset_gu/func/proc4_standard
    mkdir -p data/dataset_gu/func/proc5_mask_smooth
    mkdir -p data/dataset_gu/func/proc6_filter_norm

    echo "Functional preprocessing..."
    for file_path in data/dataset_gu/func/raw/*.nii.gz; do
        filename=$(basename $file_path)
        echo "$filename" 
        # get base subject name to specify path to structural scan
        subj_file=$(cut -d'_' -f1 <<< "${filename}")
        filename_base=$(cut -d'.' -f1 <<< "${filename}")
        # # Motion correction (re-alignment)
        # mcflirt -in $file_path -out data/dataset_gu/func/proc1_mcflirt/$filename -plots
        # # Slices were acquired bottom to top (default option in fsl)
        # slicetimer -i data/dataset_gu/func/proc1_mcflirt/$filename -o data/dataset_gu/func/proc2_slicetime/$filename -r 2.1
        # Select first volume from each functional
        # fslroi data/dataset_gu/func/proc2_slicetime/$filename data/dataset_gu/func/proc3A_firstvolume/$filename 0 1
        # # Co-registration with structural
        # epi_reg --epi=data/dataset_gu/func/proc3A_firstvolume/$filename \
        # --t1="data/dataset_gu/anat/raw/${subj_file}_T1w" --t1brain="data/dataset_gu/anat/proc1_bet/${subj_file}_T1w" \
        # --out=data/dataset_gu/func/proc3B_func2struct/$filename
        # Get transform file to send functional to MNI
        # applywarp --ref=masks/MNI152_T1_3mm_brain.nii.gz --in=data/dataset_gu/func/proc2_slicetime/$filename \
        # --out=data/dataset_gu/func/proc4_standard/$filename --warp="data/dataset_gu/anat/proc3_fnirt/${subj_file}_T1w.nii.gz.mat.nii.gz" \
        # --premat="data/dataset_gu/func/proc3B_func2struct/${filename_base}.mat" 
        # Mask
        # fslmaths data/dataset_gu/func/proc4_standard/$filename -mul $mask -kernel gauss 2.123 \
        # -fmean data/dataset_gu/func/proc5_mask_smooth/$filename
        # low pass filter (0.1Hz) and norm
        python utils/norm_filter.py -f data/dataset_gu/func/proc5_mask_smooth/$filename -m $mask \
         -ch 0.1 -t 2.1 -o data/dataset_gu/func/proc6_filter_norm/$filename



    done
fi

# Yale - FMRI, Pupillometry
if [ "$dataset" == "yale" ]; then
    mkdir -p data/dataset_yale/anat/proc1_bet
    mkdir -p data/dataset_yale/anat/proc2_affine
    mkdir -p data/dataset_yale/anat/proc3_fnirt

    # echo "Structural preprocessing..."
    # Structural preprocessing
    # for file_path in data/dataset_yale/anat/raw/*.nii.gz; do
    #     filename=$(basename $file_path)
    #     echo "$filename"
    #     # Brain Extraction
    #     bet $file_path data/dataset_yale/anat/proc1_bet/$filename -o -m
    #     # Affine registration
    #     flirt -in data/dataset_yale/anat/proc1_bet/$filename -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz \
    #     -out data/dataset_yale/anat/proc2_affine/$filename -omat data/dataset_yale/anat/proc2_affine/$filename.mat
    #     # Nonlinear transformation
    #     fnirt --ref=$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz --in=$file_path \
    #     --iout=data/dataset_yale/anat/proc3_fnirt/$filename --cout=data/dataset_yale/anat/proc3_fnirt/$filename.mat \
    #     --aff=data/dataset_yale/anat/proc2_affine/$filename.mat --config=T1_2_MNI152_2mm --warpres=6,6,6
    # done

    mkdir -p data/dataset_yale/func/proc1_mcflirt
    mkdir -p data/dataset_yale/func/proc2A_firstvolume
    mkdir -p data/dataset_yale/func/proc2B_func2struct
    mkdir -p data/dataset_yale/func/proc3_standard
    mkdir -p data/dataset_yale/func/proc4_mask_smooth
    mkdir -p data/dataset_yale/func/proc5_filter_norm
    mkdir -p data/dataset_yale/func/proc6_trim

    echo "Functional preprocessing..."
    for file_path in data/dataset_yale/func/raw/*.nii.gz; do
        filename=$(basename $file_path)
        echo "$filename" 
        # get base subject name to specify path to structural scan
        subj_file=$(cut -d'_' -f1 <<< "${filename}")
        filename_base=$(cut -d'.' -f1 <<< "${filename}")
        # # Motion correction (re-alignment)
        # mcflirt -in $file_path -out data/dataset_yale/func/proc1_mcflirt/$filename -plots
        # Select first volume from each functional
        # fslroi data/dataset_yale/func/proc1_mcflirt/$filename data/dataset_yale/func/proc2A_firstvolume/$filename 0 1
        # # Co-registration with structural
        # epi_reg --epi=data/dataset_yale/func/proc2A_firstvolume/$filename \
        # --t1="data/dataset_yale/anat/raw/${subj_file}_T1w" --t1brain="data/dataset_yale/anat/proc1_bet/${subj_file}_T1w" \
        # --out=data/dataset_yale/func/proc2B_func2struct/$filename
        # # Get transform file to send functional to MNI
        # applywarp --ref=masks/MNI152_T1_3mm_brain.nii.gz --in=data/dataset_yale/func/proc1_mcflirt/$filename \
        # --out=data/dataset_yale/func/proc3_standard/$filename --warp="data/dataset_yale/anat/proc3_fnirt/${subj_file}_T1w.nii.gz.mat.nii.gz" \
        # --premat="data/dataset_yale/func/proc2B_func2struct/${filename_base}.mat" 
        # Mask
        # fslmaths data/dataset_yale/func/proc3_standard/$filename -mul $mask -kernel gauss 2.123 \
        # -fmean data/dataset_yale/func/proc4_mask_smooth/$filename
        # low pass filter (0.1Hz) and norm
        python utils/norm_filter.py -f data/dataset_yale/func/proc4_mask_smooth/$filename -m $mask \
         -ch 0.1 -t 1 -o data/dataset_yale/func/proc5_filter_norm/$filename
        # Trim first ten volumes
        python utils/trim.py -f data/dataset_yale/func/proc5_filter_norm/$filename -o data/dataset_yale/func/proc6_trim/$filename -n 10


    done
fi

# Cheo - FMRI, EGG
if [ "$dataset" == "choe" ]; then

    mkdir -p data/dataset_choe/func/proc1_resample
    mkdir -p data/dataset_choe/func/proc2_mask
    mkdir -p data/dataset_choe/func/proc3_filter

    echo "Functional preprocessing..."
    for file_path in data/dataset_choe/func/raw/*.nii; do
        filename=$(basename $file_path)
        echo "$filename" 
        # # Resample to 3mm
        # flirt -in $file_path -ref $mask -out data/dataset_choe/func/proc1_resample/$filename -applyisoxfm 3
        # # # Mask
        # fslmaths data/dataset_choe/func/proc1_resample/$filename -mul $mask data/dataset_choe/func/proc2_mask/$filename
        # Lowpass filter (0.1Hz) and norm
        python utils/norm_filter.py -f data/dataset_choe/func/proc2_mask/$filename.gz -m $mask -t 2 \
         -ch 0.1 -o data/dataset_choe/func/proc3_filter/$filename
    done

    # mkdir -p data/dataset_choe/egg/proc1_filt
    # mkdir -p data/dataset_choe/egg/proc2_resample
    echo "EGG preprocessing..."
    sed 1d data/dataset_choe/run_list_choe.csv | while read subj || [ -n "$subj" ];
    do
        subj=$(echo $subj|tr -d '\r')
        echo "$subj"
        # Lowpass filter (0.1 Hz) and norm - signal is 10 Hz (or 0.1 TR)
        python utils/norm_filter.py -f "data/dataset_choe/egg/raw/0${subj}_run1_EGG.txt" -n 0 -t 0.1 \
         -ch 0.1 -o "data/dataset_choe/egg/proc1_filt/0${subj}_run1_EGG.txt" 
        # Resample signal
        python utils/resample_signal.py -f "data/dataset_choe/egg/proc1_filt/0${subj}_run1_EGG.txt"  \
         -s "data/dataset_choe/func/proc3_filter/pb04.20190${subj}jp.r01.blur+tlrc.nii" \
         -o "data/dataset_choe/egg/proc2_resample/0${subj}_run1_EGG.txt" 
    done 
fi


# NKI - FMRI Breath-hold task
if [ "$dataset" == "nki" ]; then
    mkdir -p data/dataset_nki/anat/proc1_bet
    mkdir -p data/dataset_nki/anat/proc2_affine
    mkdir -p data/dataset_nki/anat/proc3_fnirt

    # echo "Structural preprocessing..."
    # # Structural preprocessing
    # for file_path in data/dataset_nki/anat/raw/*.nii.gz; do
    #     filename=$(basename $file_path)
    #     echo "$filename"
    #     # Brain Extraction
    #     bet $file_path data/dataset_nki/anat/proc1_bet/$filename -o -m
    #     # Affine registration
    #     flirt -in data/dataset_nki/anat/proc1_bet/$filename -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz \
    #     -out data/dataset_nki/anat/proc2_affine/$filename -omat data/dataset_nki/anat/proc2_affine/$filename.mat
    #     # Nonlinear transformation
    #     fnirt --ref=$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz --in=$file_path \
    #     --iout=data/dataset_nki/anat/proc3_fnirt/$filename --cout=data/dataset_nki/anat/proc3_fnirt/$filename.mat \
    #     --aff=data/dataset_nki/anat/proc2_affine/$filename.mat --config=T1_2_MNI152_2mm --warpres=6,6,6
    # done

    mkdir -p data/dataset_nki/func/proc1_mcflirt
    mkdir -p data/dataset_nki/func/proc2A_firstvolume
    mkdir -p data/dataset_nki/func/proc2B_func2struct
    mkdir -p data/dataset_nki/func/proc3_standard
    mkdir -p data/dataset_nki/func/proc4_mask_smooth
    mkdir -p data/dataset_nki/func/proc5_filter_norm

    echo "Functional preprocessing..."
    for file_path in data/dataset_nki/func/raw/*.nii.gz; do
        filename=$(basename $file_path)
        echo "$filename" 
        # get base subject name to specify path to structural scan
        subj_file=$(cut -d'_' -f1 <<< "${filename}")
        filename_base=$(cut -d'.' -f1 <<< "${filename}")
        # # Motion correction (re-alignment)
        # mcflirt -in $file_path -out data/dataset_nki/func/proc1_mcflirt/$filename -plots
        # # Select first volume from each functional
        # fslroi data/dataset_nki/func/proc1_mcflirt/$filename data/dataset_nki/func/proc2A_firstvolume/$filename 0 1
        # # # Co-registration with structural
        # epi_reg --epi=data/dataset_nki/func/proc2A_firstvolume/$filename \
        # --t1="data/dataset_nki/anat/raw/${subj_file}_T1w" --t1brain="data/dataset_nki/anat/proc1_bet/${subj_file}_T1w" \
        # --out=data/dataset_nki/func/proc2B_func2struct/$filename
        # # Get transform file to send functional to MNI
        # applywarp --ref=masks/MNI152_T1_3mm_brain.nii.gz --in=data/dataset_nki/func/proc1_mcflirt/$filename \
        # --out=data/dataset_nki/func/proc3_standard/$filename --warp="data/dataset_nki/anat/proc3_fnirt/${subj_file}_T1w.nii.gz.mat.nii.gz" \
        # --premat="data/dataset_nki/func/proc2B_func2struct/${filename_base}.mat" 
        # # Mask
        # fslmaths data/dataset_nki/func/proc3_standard/$filename -mul $mask -kernel gauss 2.123 \
        # -fmean data/dataset_nki/func/proc4_mask_smooth/$filename
        # low filter (0.1Hz) and norm
        python utils/norm_filter.py -f data/dataset_nki/func/proc4_mask_smooth/$filename -m $mask \
         -ch 0.1 -t 1.4 -o data/dataset_nki/func/proc5_filter_norm/$filename

    done
fi

# HCP - FMRI, hr, rv
if [ "$dataset" == "hcp" ]; then
    # Make directories
    mkdir -p data/dataset_hcp/func/proc1_resample
    mkdir -p data/dataset_hcp/func/proc2_mask_smooth
    mkdir -p data/dataset_hcp/func/proc3_filter_norm

    for file_path in data/dataset_hcp/func/raw/*.nii.gz; do
        filename=$(basename $file_path)
        echo "$filename" 

        # Resample to 3mm
        flirt -in $file_path -ref $mask -out data/dataset_hcp/func/proc1_resample/$filename -applyisoxfm 3
        # Mask and smooth
        fslmaths data/dataset_hcp/func/proc1_resample/$filename -mul $mask -kernel gauss 2.123 \
        -fmean data/dataset_hcp/func/proc2_mask_smooth/$filename
        # Lowpass filter (0.1Hz) and norm
        python utils/norm_filter.py -f data/dataset_hcp/func/proc2_mask_smooth/$filename -m $mask \
         -ch 0.1 -t 2.1 -o data/dataset_hcp/func/proc3_filter_norm/$filename
    done
fi

    # fslmaths data/bold/proc1_resample/$filename -kernel gauss 1.27399 -fmean data/bold/proc2_smooth_mask/$filename
