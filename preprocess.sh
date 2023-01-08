#!/bin/bash
dataset=$1

# Dilated mask that includes sinuses slightly outside gray matter tissue
mask="masks/MNI152_T1_3mm_brain_mask_dilated.nii.gz"

# Chang - EEG, FMRI, physio data
if [ "$dataset" == "chang" ]; then

    # mkdir -p data/dataset_chang/anat/proc1_reorient
    # mkdir -p data/dataset_chang/anat/proc2_bet
    # mkdir -p data/dataset_chang/anat/proc3_affine
    # mkdir -p data/dataset_chang/anat/proc4_fnirt
    # mkdir -p data/dataset_chang/anat/proc5_csfmask
    # echo "Structural preprocessing..."
    # for file_path in data/dataset_chang/anat/raw/*.nii; do
    #     filename=$(basename $file_path)
    #     echo "$filename"
    #     filename_base=$(cut -d'.' -f1 <<< "${filename}")
    #     subj_name=$(cut -d'-' -f1 <<< "${filename_base}")
        # # Reorient to standard orientation
        # fslreorient2std $file_path data/dataset_chang/anat/proc1_reorient/$filename
        # # Brain Extraction
        # bet data/dataset_chang/anat/proc1_reorient/$filename data/dataset_chang/anat/proc2_bet/$filename -o -m
        # # FAST segmentation
        # fast data/dataset_chang/anat/proc2_bet/$filename 
        # # Affine registration
        # flirt -in data/dataset_chang/anat/proc2_bet/$filename -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz \
        # -out data/dataset_chang/anat/proc3_affine/$filename -omat data/dataset_chang/anat/proc3_affine/$filename.mat
        # # Nonlinear transformation
        # fnirt --ref=$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz --in=data/dataset_chang/anat/proc1_reorient/$filename \
        # --iout=data/dataset_chang/anat/proc4_fnirt/$filename --cout=data/dataset_chang/anat/proc4_fnirt/$filename.mat \
        # --aff=data/dataset_chang/anat/proc3_affine/$filename.mat --config=T1_2_MNI152_2mm --warpres=6,6,6
        # # Send CSF mask to MNI space
        # applywarp --ref=$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz --in=data/dataset_chang/anat/proc2_bet/${subj_name}-anat_pve_0  \
        # --out=data/dataset_chang/anat/proc5_csfmask/$filename_base --warp="data/dataset_chang/anat/proc4_fnirt/${subj_name}-anat.nii.mat.nii.gz" 
        # Binarize CSF mask
    #     fslmaths data/dataset_chang/anat/proc5_csfmask/$filename_base -thr 0.9 -bin \
    #     data/dataset_chang/anat/proc5_csfmask/$filename_base
    # done

    # # Create group CSF mask
    # # Add together subject CSF masks in MNI space
    # fslmaths data/dataset_chang/anat/proc5_csfmask/sub_0010-anat.nii.gz -add data/dataset_chang/anat/proc5_csfmask/sub_0016-anat.nii.gz \
    # -add data/dataset_chang/anat/proc5_csfmask/sub_0028-anat.nii.gz -add data/dataset_chang/anat/proc5_csfmask/sub_0031-anat.nii.gz \
    # -add data/dataset_chang/anat/proc5_csfmask/sub_0012-anat.nii.gz -add data/dataset_chang/anat/proc5_csfmask/sub_0021-anat.nii.gz \
    # -add data/dataset_chang/anat/proc5_csfmask/sub_0029-anat.nii.gz -add data/dataset_chang/anat/proc5_csfmask/sub_0032-anat.nii.gz \
    # -add data/dataset_chang/anat/proc5_csfmask/sub_0014-anat.nii.gz -add data/dataset_chang/anat/proc5_csfmask/sub_0027-anat.nii.gz \
    # -add data/dataset_chang/anat/proc5_csfmask/sub_0030-anat.nii.gz data/dataset_chang/anat/proc5_csfmask/group_csf_mask
    # # Threshold to >= 10 subjects overlap in CSF masks
    # fslmaths data/dataset_chang/anat/proc5_csfmask/group_csf_mask -thr 10 -bin data/dataset_chang/anat/proc5_csfmask/group_csf_mask
    # # Mask by MNI prior probability CSF mask (thresholded)
    # fslmaths data/dataset_chang/anat/proc5_csfmask/group_csf_mask -mul masks/MNI152_T1_2mm_csf_mask \
    # data/dataset_chang/anat/proc5_csfmask/group_csf_mask

    # mkdir -p data/dataset_chang/func/proc1_resample
    # mkdir -p data/dataset_chang/func/proc2_smooth_mask
    # mkdir -p data/dataset_chang/func/proc3_filter_norm
    # mkdir -p data/dataset_chang/func/proc4_bandpass

    # echo "Functional preprocessing - post ME-ICA..."
    # for file_path in data/dataset_chang/func/raw/*.nii.gz; do
    #     filename=$(basename $file_path)
    #     echo "$filename" 
    #     # Resample to 3mm
    #     flirt -in $file_path -ref $mask -out data/dataset_chang/func/proc1_resample/$filename -applyisoxfm 3
    #     # Mask
    #     fslmaths data/dataset_chang/func/proc1_resample/$filename -mul $mask data/dataset_chang/func/proc2_smooth_mask/$filename
    #     # Lowpass filter (0.1Hz) and norm
    #     python -m utils.signal.norm_filter -f data/dataset_chang/func/proc2_smooth_mask/$filename -m $mask \
    #      -ch 0.1 -t 2.1 -o data/dataset_chang/func/proc3_filter_norm/$filename
    #     # band pass filter (0.01-0.1Hz) and norm 
    #     python -m utils.signal.norm_filter -f data/dataset_chang/func/proc2_smooth_mask/${filename} -m $mask \
    #      -ch 0.1 -cl 0.01 -t 2.1 -o data/dataset_chang/func/proc4_bandpass/$filename
    # done


    mkdir -p data/dataset_chang/physio/proc1_physio
    mkdir -p data/dataset_chang/eeg/proc1_fbands
    echo "Physio preprocessing..."
    for file_path in data/dataset_chang/func/raw/*.nii.gz; do
        filename=$(basename $file_path)
        echo "$filename" 
        # get base subject name to specify path to mask
        subj_file=$(cut -d'-' -f1 <<< "${filename}")
        sess_n=$(cut -d'-' -f2 <<< "${filename}")
        filename_base=$(cut -d'.' -f1 <<< "${filename}")
        subj_out=${subj_file}_${sess_n}
        # Preprocess EEG and Physio Data
        python -m utils.dataset.preprocess_chang -e data/dataset_chang/eeg/raw/${subj_file}-${sess_n}_eeg_pp.mat \
        -p data/dataset_chang/physio/raw/${subj_file}-${sess_n}-ecr_echo1_physOUT.mat \
         -f 693 -om data/dataset_chang/eeg/raw/${subj_out} \
         -oe data/dataset_chang/eeg/proc1_fbands/${subj_out}_fbands \
         -op data/dataset_chang/physio/proc1_physio/${subj_out}_physio 
        # # Extract global BOLD signal from smoothed functional data
        # fslmeants -i data/dataset_chang/func/proc2_smooth_mask/${filename}\
        # -o data/dataset_chang/physio/proc1_physio/${subj_out}_global_sig.txt \
        # -m masks/MNI152_T1_3mm_gray_mask.nii.gz
        # # Extract CSF signal from raw functional data (pre-ME-ICA)
        # fslmeants -i data/dataset_chang/func/raw/$filename \
        # -o data/dataset_chang/physio/proc1_physio/${subj_out}_csf.txt \
        # -m data/dataset_chang/anat/proc5_csfmask/group_csf_mask
        # # Extract precuneus BOLD signal from preprocessed band-pass functional data
        # fslmeants -i data/dataset_chang/func/proc4_bandpass/${filename}\
        # -o data/dataset_chang/physio/proc1_physio/${subj_out}_precuneus.txt \
        # -m masks/precuneus_sphere_6mm.nii.gz
        # # Extract superior parietal BOLD signal from preprocessed band-pass functional data
        # fslmeants -i data/dataset_chang/func/proc4_bandpass/${filename}\
        # -o data/dataset_chang/physio/proc1_physio/${subj_out}_superior_parietal.txt \
        # -m masks/superior_parietal_sphere_6mm.nii.gz
    done
fi


# Chang Breath Hold Task - EEG, FMRI, physio data for breath hold task
if [ "$dataset" == "chang_bh" ]; then
    # mkdir -p data/dataset_chang_bh/anat/proc1_reorient
    # mkdir -p data/dataset_chang_bh/anat/proc2_bet
    # mkdir -p data/dataset_chang_bh/anat/proc3_affine
    # mkdir -p data/dataset_chang_bh/anat/proc4_fnirt
    # mkdir -p data/dataset_chang_bh/anat/proc5_csfmask
    # echo "Structural preprocessing..."
    # for file_path in data/dataset_chang_bh/anat/raw/*.nii; do
    #     filename=$(basename $file_path)
    #     echo "$filename"
    #     filename_base=$(cut -d'.' -f1 <<< "${filename}")
    #     subj_name=$(cut -d'-' -f1 <<< "${filename_base}")
        # # Reorient to standard orientation
        # fslreorient2std $file_path data/dataset_chang_bh/anat/proc1_reorient/$filename
        # # Brain Extraction
        # bet data/dataset_chang_bh/anat/proc1_reorient/$filename data/dataset_chang_bh/anat/proc2_bet/$filename -o -m
        # # FAST segmentation
        # fast data/dataset_chang_bh/anat/proc2_bet/$filename 
        # # Affine registration
        # flirt -in data/dataset_chang_bh/anat/proc2_bet/$filename -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz \
        # -out data/dataset_chang_bh/anat/proc3_affine/$filename -omat data/dataset_chang_bh/anat/proc3_affine/$filename.mat
        # # Nonlinear transformation
        # fnirt --ref=$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz --in=data/dataset_chang_bh/anat/proc1_reorient/$filename \
        # --iout=data/dataset_chang_bh/anat/proc4_fnirt/$filename --cout=data/dataset_chang_bh/anat/proc4_fnirt/$filename.mat \
        # --aff=data/dataset_chang_bh/anat/proc3_affine/$filename.mat --config=T1_2_MNI152_2mm --warpres=6,6,6
        # # Send CSF mask to MNI space
        # applywarp --ref=$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz --in=data/dataset_chang_bh/anat/proc2_bet/${subj_name}-anat_pve_0  \
        # --out=data/dataset_chang_bh/anat/proc5_csfmask/$filename_base --warp="data/dataset_chang_bh/anat/proc4_fnirt/${subj_name}-anat.nii.mat.nii.gz" 
        # # Binarize CSF mask
        # fslmaths data/dataset_chang_bh/anat/proc5_csfmask/$filename_base -thr 0.9 -bin \
        # data/dataset_chang_bh/anat/proc5_csfmask/$filename_base
    # done

    # # Create group CSF mask
    # # Add together subject CSF masks in MNI space
    # fslmaths data/dataset_chang_bh/anat/proc5_csfmask/sub_0015-anat.nii.gz -add data/dataset_chang_bh/anat/proc5_csfmask/sub_0017-anat.nii.gz \
    # -add data/dataset_chang_bh/anat/proc5_csfmask/sub_0019-anat.nii.gz -add data/dataset_chang_bh/anat/proc5_csfmask/sub_0020-anat.nii.gz \
    # -add data/dataset_chang_bh/anat/proc5_csfmask/sub_0021-anat.nii.gz data/dataset_chang_bh/anat/proc5_csfmask/group_csf_mask
    # # Threshold to >4 subjects overlap in CSF masks
    # fslmaths data/dataset_chang_bh/anat/proc5_csfmask/group_csf_mask -thr 4 -bin data/dataset_chang_bh/anat/proc5_csfmask/group_csf_mask
    # # Mask by MNI prior probability CSF mask (thresholded)
    # fslmaths data/dataset_chang_bh/anat/proc5_csfmask/group_csf_mask -mul masks/MNI152_T1_2mm_csf_mask \
    # data/dataset_chang_bh/anat/proc5_csfmask/group_csf_mask

    # mkdir -p data/dataset_chang_bh/func/proc1_resample
    # mkdir -p data/dataset_chang_bh/func/proc2_smooth_mask
    # mkdir -p data/dataset_chang_bh/func/proc3_filter_norm
    # mkdir -p data/dataset_chang_bh/func/proc4_bandpass

    # echo "Functional preprocessing - post ME-ICA..."
    # for file_path in data/dataset_chang_bh/func/raw/*.nii.gz; do
    #     filename=$(basename $file_path)
    #     echo "$filename" 
    #     # Resample to 3mm
    #     flirt -in $file_path -ref $mask -out data/dataset_chang_bh/func/proc1_resample/$filename -applyisoxfm 3
    #     # Mask
    #     fslmaths data/dataset_chang_bh/func/proc1_resample/$filename -mul $mask data/dataset_chang_bh/func/proc2_smooth_mask/$filename
    #     # Lowpass filter (0.1Hz) and norm
    #     python -m utils.signal.norm_filter -f data/dataset_chang_bh/func/proc2_smooth_mask/$filename -m $mask \
    #      -ch 0.1 -t 2.1 -o data/dataset_chang_bh/func/proc3_filter_norm/$filename
    #     # band pass filter (0.01-0.1Hz) and norm 
    #     python -m utils.signal.norm_filter -f data/dataset_chang_bh/func/proc2_smooth_mask/${filename} -m $mask \
    #      -ch 0.1 -cl 0.01 -t 2.1 -o data/dataset_chang_bh/func/proc4_bandpass/$filename
    # done


    mkdir -p data/dataset_chang_bh/physio/proc1_physio
    mkdir -p data/dataset_chang_bh/eeg/proc1_fbands
    echo "Physio preprocessing..."
    sed -n '2,$p' data/dataset_chang_bh/subject_list_chang_bh.csv | while IFS=, read -r subj scan nframes;
    do 
        if [[ $scan -gt 10 ]]
        then
          sess_n="mr_00${scan}"
        else
          sess_n="mr_000${scan}"
        fi
        subj_file="sub_00${subj}"
        subj_out=${subj_file}_${sess_n}
        echo "$subj_file $sess_n"
        # Preprocess EEG and Physio Data
        python -m utils.dataset.preprocess_chang \
        -e data/dataset_chang_bh/eeg/raw/${subj_file}-${sess_n}-adb_echo1_EEG_pp.mat \
        -p data/dataset_chang_bh/physio/raw/${subj_file}-${sess_n}-adb_echo1_physOUT.mat \
         -f $nframes -om data/dataset_chang_bh/eeg/raw/${subj_out} \
         -oe data/dataset_chang_bh/eeg/proc1_fbands/${subj_out}_fbands \
         -op data/dataset_chang_bh/physio/proc1_physio/${subj_out}_physio  
        # # Extract global BOLD signal from preprocessed low-pass functional data
        # fslmeants -i data/dataset_chang_bh/func/proc3_filter_norm/${subj_file}-${sess_n}-adb_echo1_w_dspk_blur3mm \
        # -o data/dataset_chang_bh/physio/proc1_physio/${subj_out}_global_sig.txt \
        # -m masks/MNI152_T1_3mm_gray_mask.nii.gz
        # # Extract CSF signal from raw functional data (post ME-ICA)
        # fslmeants -i data/dataset_chang_bh/func/raw/${subj_file}-${sess_n}-adb_echo1_w_dspk_blur3mm  \
        # -o data/dataset_chang_bh/physio/proc1_physio/${subj_out}_csf.txt \
        # -m data/dataset_chang_bh/anat/proc5_csfmask/group_csf_mask
    done
fi

# NKI - FMRI Breath-hold task
if [ "$dataset" == "nki" ]; then
    # mkdir -p data/dataset_nki/anat/proc1_bet
    # mkdir -p data/dataset_nki/anat/proc2_affine
    # mkdir -p data/dataset_nki/anat/proc3_fnirt
    # mkdir -p data/dataset_nki/anat/proc4_csfmask

    # echo "Structural preprocessing..."
    # # Structural preprocessing
    # for file_path in data/dataset_nki/anat/raw/*.nii.gz; do
    #     filename=$(basename $file_path)
    #     filename_base=$(cut -d'.' -f1 <<< "${filename}")
    #     echo "$filename"
        # # Brain Extraction
        # bet $file_path data/dataset_nki/anat/proc1_bet/$filename -o -m
        # # FAST segmentation
        # fast data/dataset_nki/anat/proc1_bet/$filename_base
        # # Affine registration
        # flirt -in data/dataset_nki/anat/proc1_bet/$filename -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz \
        # -out data/dataset_nki/anat/proc2_affine/$filename -omat data/dataset_nki/anat/proc2_affine/$filename.mat
        # # Nonlinear transformation
        # fnirt --ref=$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz --in=$file_path \
        # --iout=data/dataset_nki/anat/proc3_fnirt/$filename --cout=data/dataset_nki/anat/proc3_fnirt/$filename.mat \
        # --aff=data/dataset_nki/anat/proc2_affine/$filename.mat --config=T1_2_MNI152_2mm --warpres=6,6,6
        # # Send CSF mask to MNI space
        # applywarp --ref=masks/MNI152_T1_3mm_brain.nii.gz --in=data/dataset_nki/anat/proc1_bet/${filename_base}_pve_0  \
        # --out=data/dataset_nki/anat/proc4_csfmask/$filename_base --warp="data/dataset_nki/anat/proc3_fnirt/${filename}.mat.nii.gz" 
        # # Binarize CSF mask
        # fslmaths data/dataset_nki/anat/proc4_csfmask/$filename_base -thr 0.9 -bin \
        # data/dataset_nki/anat/proc4_csfmask/$filename_base

    # done

    # # Create group CSF mask
    # ## remove mask if exists
    # rm -f data/dataset_nki/anat/proc4_csfmask/group_csfmask.nii.gz
    # ## Create empty mask
    # tmp_files=(data/dataset_nki/anat/proc4_csfmask/*)    
    # fslmaths ${tmp_files[0]} -thr 2 data/dataset_nki/anat/proc4_csfmask/group_csfmask
    # ## Loop through subjects and add masks
    # for file_path in data/dataset_nki/anat/raw/*.nii.gz; do
    #     filename=$(basename $file_path)
    #     filename_base=$(cut -d'.' -f1 <<< "${filename}")
    #     fslmaths data/dataset_nki/anat/proc4_csfmask/group_csfmask \
    #     -add data/dataset_nki/anat/proc4_csfmask/${filename_base} \
    #     data/dataset_nki/anat/proc4_csfmask/group_csfmask
    # done

    # # Threshold to >40 subjects overlap in CSF masks
    # fslmaths data/dataset_nki/anat/proc4_csfmask/group_csfmask -thr 40 -bin \
    # data/dataset_nki/anat/proc4_csfmask/group_csfmask

    # # # Mask by MNI prior probability CSF mask (thresholded)
    # fslmaths data/dataset_nki/anat/proc4_csfmask/group_csfmask -mul masks/MNI152_T1_3mm_csf_mask \
    # data/dataset_nki/anat/proc4_csfmask/group_csfmask

    # mkdir -p data/dataset_nki/func/proc1_mcflirt
    # mkdir -p data/dataset_nki/func/proc2A_firstvolume
    # mkdir -p data/dataset_nki/func/proc2B_func2struct
    # mkdir -p data/dataset_nki/func/proc3_standard
    # mkdir -p data/dataset_nki/func/proc4_mask_smooth
    # mkdir -p data/dataset_nki/func/proc5_filter_norm

    # echo "Functional preprocessing..."
    # for file_path in data/dataset_nki/func/raw/*.nii.gz; do
    #     filename=$(basename $file_path)
    #     echo "$filename" 
    #     # get base subject name to specify path to structural scan
    #     subj_file=$(cut -d'_' -f1 <<< "${filename}")
    #     filename_base=$(cut -d'.' -f1 <<< "${filename}")
    #     # Motion correction (re-alignment)
    #     mcflirt -in $file_path -out data/dataset_nki/func/proc1_mcflirt/$filename -plots
    #     # Select first volume from each functional
    #     fslroi data/dataset_nki/func/proc1_mcflirt/$filename data/dataset_nki/func/proc2A_firstvolume/$filename 0 1
    #     # # Co-registration with structural
    #     epi_reg --epi=data/dataset_nki/func/proc2A_firstvolume/$filename \
    #     --t1="data/dataset_nki/anat/raw/${subj_file}_T1w" --t1brain="data/dataset_nki/anat/proc1_bet/${subj_file}_T1w" \
    #     --out=data/dataset_nki/func/proc2B_func2struct/$filename
    #     # Get transform file to send functional to MNI
    #     applywarp --ref=masks/MNI152_T1_3mm_brain.nii.gz --in=data/dataset_nki/func/proc1_mcflirt/$filename \
    #     --out=data/dataset_nki/func/proc3_standard/$filename --warp="data/dataset_nki/anat/proc3_fnirt/${subj_file}_T1w.nii.gz.mat.nii.gz" \
    #     --premat="data/dataset_nki/func/proc2B_func2struct/${filename_base}.mat" 
    #     # Mask
    #     fslmaths data/dataset_nki/func/proc3_standard/$filename -mul $mask -kernel gauss 2.123 \
    #     -fmean data/dataset_nki/func/proc4_mask_smooth/$filename
    #     # low filter (0.1Hz) and norm
    #     python -m utils.signal.norm_filter -f data/dataset_nki/func/proc4_mask_smooth/$filename -m $mask \
    #      -ch 0.1 -t 1.4 -o data/dataset_nki/func/proc5_filter_norm/$filename
    # done

    mkdir -p data/dataset_nki/physio/proc1_physio
    echo "Physio preprocessing..."
    for file_path in data/dataset_nki/physio/raw/*.tsv.gz; do
        filename=$(basename $file_path)
        echo "$filename" 
        # get base subject name to specify path to structural scan
        subj_file=$(cut -d'_' -f1 <<< "${filename}")
        filename_base=$(cut -d'.' -f1 <<< "${filename}")
        # Physio extraction
        python -m utils.dataset.preprocess_nki -s $file_path -o data/dataset_nki/physio/proc1_physio/$filename_base
        # # CSF extraction
        # fslmeants -i data/dataset_nki/func/proc3_standard/${subj_file}_task_breathhold \
        # -o data/dataset_nki/physio/proc1_physio/${subj_file}_task_breathhold_physio_csf.txt \
        # -m data/dataset_nki/anat/proc4_csfmask/group_csfmask
    done

fi

# HCP Resting - FMRI, hr, rv
if [ "$dataset" == "hcp" ]; then

    # mkdir -p data/dataset_hcp/anat/proc1_csfmask
    # echo "Structural preprocessing..."
    # # Structural preprocessing
    # for file_path in data/dataset_hcp/anat/raw_fs_seg/*.nii.gz; do
    #     filename=$(basename $file_path)
    #     echo "$filename" 
    #     filename_base=$(cut -d'.' -f1 <<< "${filename}")
    #     # Extract left/right lateral ventricle and 4th ventricle
    #     ## Left lateral ventricle = 4
    #     fslmaths $file_path -thr 4 -uthr 4 data/dataset_hcp/anat/proc1_csfmask/${filename_base}_l_lv
    #     ## rigth lateral ventricle = 43
    #     fslmaths $file_path -thr 43 -uthr 43 data/dataset_hcp/anat/proc1_csfmask/${filename_base}_r_lv
    #     ## fourth ventricle = 15
    #     fslmaths $file_path -thr 15 -uthr 15 data/dataset_hcp/anat/proc1_csfmask/${filename_base}_4v
    #     ## Create csf mask
    #     fslmaths data/dataset_hcp/anat/proc1_csfmask/${filename_base}_l_lv \
    #     -add data/dataset_hcp/anat/proc1_csfmask/${filename_base}_r_lv \
    #     -add data/dataset_hcp/anat/proc1_csfmask/${filename_base}_4v -bin \
    #     data/dataset_hcp/anat/proc1_csfmask/${filename_base}_csfmask
    #     ## Send csf mask to 2mm MNI space
    #     flirt -in data/dataset_hcp/anat/proc1_csfmask/${filename_base}_csfmask \
    #     -ref $FSLDIR/data/standard/MNI152_T1_2mm.nii.gz \
    #     -out data/dataset_hcp/anat/proc1_csfmask/${filename_base}_csfmask -applyxfm -interp nearestneighbour
    #     # Clean up csf masks
    #     rm data/dataset_hcp/anat/proc1_csfmask/${filename_base}_l_lv.nii.gz
    #     rm data/dataset_hcp/anat/proc1_csfmask/${filename_base}_r_lv.nii.gz
    #     rm data/dataset_hcp/anat/proc1_csfmask/${filename_base}_4v.nii.gz
    # done

    # # Create group CSF mask
    # ## remove mask if exists
    # rm -f data/dataset_hcp/anat/proc1_csfmask/group_csfmask.nii.gz
    # ## Create empty mask
    # tmp_files=(data/dataset_hcp/anat/proc1_csfmask/*)    
    # fslmaths ${tmp_files[0]} -thr 100 data/dataset_hcp/anat/proc1_csfmask/group_csfmask
    # ## Loop through subjects and add masks
    # for file_path in data/dataset_hcp/anat/raw_fs_seg/*.nii.gz; do
    #     filename=$(basename $file_path)
    #     filename_base=$(cut -d'.' -f1 <<< "${filename}")
    #     fslmaths data/dataset_hcp/anat/proc1_csfmask/group_csfmask \
    #     -add data/dataset_hcp/anat/proc1_csfmask/${filename_base}_csfmask \
    #     data/dataset_hcp/anat/proc1_csfmask/group_csfmask
    # done

    # # Threshold to >28 subjects overlap in CSF masks
    # fslmaths data/dataset_hcp/anat/proc1_csfmask/group_csfmask -thr 28 -bin \
    # data/dataset_hcp/anat/proc1_csfmask/group_csfmask


    # mkdir -p data/dataset_hcp/func/proc1_resample
    # mkdir -p data/dataset_hcp/func/proc2_mask_smooth
    # mkdir -p data/dataset_hcp/func/proc3_filter_norm
    # mkdir -p data/dataset_hcp/func/proc4_bandpass

    # echo "Functional preprocessing for minimally preprocessed data..."
    # for file_path in data/dataset_hcp/func/raw/*.nii.gz; do
    #     filename=$(basename $file_path)
    #     echo "$filename" 
    #     # Resample to 3mm
    #     flirt -in $file_path -ref $mask -out data/dataset_hcp/func/proc1_resample/$filename -applyisoxfm 3
    #     # Mask and smooth
    #     fslmaths data/dataset_hcp/func/proc1_resample/$filename -mul $mask -kernel gauss 2.123 \
    #     -fmean data/dataset_hcp/func/proc2_mask_smooth/$filename
    #     # Lowpass filter (0.1Hz) and norm
    #     python -m utils.signal.norm_filter -f data/dataset_hcp/func/proc2_mask_smooth/$filename -m $mask \
    #      -ch 0.1 -t 0.72 -o data/dataset_hcp/func/proc3_filter_norm/$filename
    #     # band pass filter (0.01-0.1Hz)
    #     python -m utils.signal.norm_filter -f data/dataset_hcp/func/proc2_mask_smooth/$filename -m $mask \
    #      -ch 0.1 -cl 0.01 -t 0.72 -o data/dataset_hcp/func/proc4_bandpass/$filename
    # done

    # mkdir -p data/dataset_hcp/func_fix/proc1_resample
    # mkdir -p data/dataset_hcp/func_fix/proc2_mask_smooth
    # mkdir -p data/dataset_hcp/func_fix/proc3_filter_norm
    # mkdir -p data/dataset_hcp/func_fix/proc4_bandpass

    # echo "Functional preprocessing for ICA-FIX data..."
    # for file_path in data/dataset_hcp/func_fix/raw/*.nii.gz; do
    #     filename=$(basename $file_path)
    #     echo "$filename" 
    #     # Resample to 3mm
    #     flirt -in $file_path -ref $mask -out data/dataset_hcp/func_fix/proc1_resample/$filename -applyisoxfm 3
    #     # Mask and smooth
    #     fslmaths data/dataset_hcp/func_fix/proc1_resample/$filename -mul $mask -kernel gauss 2.123 \
    #     -fmean data/dataset_hcp/func_fix/proc2_mask_smooth/$filename
    #     # Lowpass filter (0.1Hz) and norm
    #     python -m utils.signal.norm_filter -f data/dataset_hcp/func_fix/proc2_mask_smooth/$filename -m $mask \
    #      -ch 0.1 -t 0.72 -o data/dataset_hcp/func_fix/proc3_filter_norm/$filename
    #     # band pass filter (0.01-0.1Hz)
    #     python -m utils.signal.norm_filter -f data/dataset_hcp/func_fix/proc2_mask_smooth/$filename -m $mask \
    #      -ch 0.1 -cl 0.01 -t 0.72 -o data/dataset_hcp/func_fix/proc4_bandpass/$filename
    # done

    mkdir -p data/dataset_hcp/physio/proc1_physio
    echo "Physio preprocessing..."
    for file_path in data/dataset_hcp/physio/raw/*.txt; do
        filename=$(basename $file_path)
        echo "$filename" 
        # get base subject name to specify path to structural scan
        subj_file=$(cut -d'_' -f1 <<< "${filename}")
        # subj_func=$(cut -d'_physio' -f1 <<< "${filename}")
        subj_func=$(echo $filename | awk 'BEGIN {FS="_physio.txt" } ; { print $1 }')

        # Physio extraction
        python -m utils.dataset.preprocess_hcp -s $file_path -o data/dataset_hcp/physio/proc1_physio/${subj_file}_physio -d rest
        # # Extract precuneus BOLD signal from preprocessed band-pass functional data
        # fslmeants -i data/dataset_hcp/func_fix/proc4_bandpass/${subj_func}_rest.nii.gz\
        # -o data/dataset_hcp/physio/proc1_physio/${subj_file}_precuneus.txt \
        # -m masks/precuneus_sphere_6mm.nii.gz
        # # Extract global BOLD signal from preprocessed low-pass functional data
        # fslmeants -i data/dataset_hcp/func_fix/proc3_filter_norm/${subj_func}_rest.nii.gz \
        # -o data/dataset_hcp/physio/proc1_physio/${subj_file}_global_sig.txt \
        # -m masks/MNI152_T1_3mm_gray_mask.nii.gz
        # # Extract superior parietal BOLD signal from preprocessed band-pass functional data
        # fslmeants -i data/dataset_hcp/func_fix/proc4_bandpass/${subj_func}_rest.nii.gz\
        # -o data/dataset_hcp/physio/proc1_physio/${subj_file}_superior_parietal.txt \
        # -m masks/superior_parietal_sphere_6mm.nii.gz
        
        # # Extract CSF signal from HCP CSF mask
        # fslmeants -i data/dataset_hcp/func_fix/raw/${subj_func}_rest.nii.gz\
        # -o data/dataset_hcp/physio/proc1_physio/${subj_file}_csf.txt \
        # -m data/dataset_hcp/anat/proc1_csfmask/group_csfmask

    done
fi

# Spreng - FMRI, Physio
if [ "$dataset" == "spreng" ]; then
    # mkdir -p data/dataset_spreng/anat/raw_reorient
    # mkdir -p data/dataset_spreng/anat/proc1_bet
    # mkdir -p data/dataset_spreng/anat/proc2_affine
    # mkdir -p data/dataset_spreng/anat/proc3_fnirt
    # mkdir -p data/dataset_spreng/anat/proc4_csfmask
    # echo "Structural preprocessing..."
    # # Structural preprocessing
    # for file_path in data/dataset_spreng/anat/raw/*.nii.gz; do
    #     filename=$(basename $file_path)
    #     echo "$filename"
    #     filename_base=$(cut -d'.' -f1 <<< "${filename}")
        # # Reorient to standard orientation
        # fslreorient2std $file_path data/dataset_spreng/anat/raw_reorient/$filename
        # # Brain Extraction
        # bet data/dataset_spreng/anat/raw_reorient/$filename data/dataset_spreng/anat/proc1_bet/$filename -o -m
        # FAST segmentation
        # fast data/dataset_spreng/anat/proc1_bet/$filename_base
        # # Affine registration
        # flirt -in data/dataset_spreng/anat/proc1_bet/$filename -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz \
        # -out data/dataset_spreng/anat/proc2_affine/$filename -omat data/dataset_spreng/anat/proc2_affine/$filename.mat
        # # Nonlinear transformation
        # fnirt --ref=$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz --in=data/dataset_spreng/anat/raw_reorient/$filename \
        # --iout=data/dataset_spreng/anat/proc3_fnirt/$filename --cout=data/dataset_spreng/anat/proc3_fnirt/$filename.mat \
        # --aff=data/dataset_spreng/anat/proc2_affine/$filename.mat --config=T1_2_MNI152_2mm --warpres=6,6,6
        # Send CSF mask to MNI space
        # applywarp --ref=masks/MNI152_T1_3mm_brain.nii.gz --in=data/dataset_spreng/anat/proc1_bet/${filename_base}_pve_0  \
        # --out=data/dataset_spreng/anat/proc4_csfmask/$filename_base --warp="data/dataset_spreng/anat/proc3_fnirt/${filename}.mat.nii.gz" 
        # # Binarize CSF mask
        # fslmaths data/dataset_spreng/anat/proc4_csfmask/$filename_base -thr 0.9 -bin \
        # data/dataset_spreng/anat/proc4_csfmask/$filename_base
    # done

    # # # Create group CSF mask
    # ## remove mask if exists
    # rm -f data/dataset_spreng/anat/proc4_csfmask/group_csfmask.nii.gz
    # ## Create empty mask
    # tmp_files=(data/dataset_spreng/anat/proc4_csfmask/*)    
    # fslmaths ${tmp_files[0]} -thr 2 data/dataset_spreng/anat/proc4_csfmask/group_csfmask
    # ## Loop through subjects and add masks
    # sed -n '2,$p' data/dataset_spreng/subject_list_spreng.csv | while IFS=, read -r id subj other_cols; do
    #     subj_fp=${subj}_ses-1_T1w
    #     fslmaths data/dataset_spreng/anat/proc4_csfmask/group_csfmask \
    #     -add data/dataset_spreng/anat/proc4_csfmask/${subj_fp} \
    #     data/dataset_spreng/anat/proc4_csfmask/group_csfmask
    # done

    # # Threshold to >41 subjects overlap in CSF masks
    # fslmaths data/dataset_spreng/anat/proc4_csfmask/group_csfmask -thr 41 -bin \
    # data/dataset_spreng/anat/proc4_csfmask/group_csfmask

    # # # Mask by MNI prior probability CSF mask (thresholded)
    # fslmaths data/dataset_spreng/anat/proc4_csfmask/group_csfmask -mul masks/MNI152_T1_3mm_csf_mask \
    # data/dataset_spreng/anat/proc4_csfmask/group_csfmask

    mkdir -p data/dataset_spreng/func/proc1_resample
    mkdir -p data/dataset_spreng/func/proc2_mask_smooth
    mkdir -p data/dataset_spreng/func/proc3_bandpass
    echo "Functional preprocessing..."
    for file_path in data/dataset_spreng/func/raw/*.nii; do
        filename=$(basename $file_path)
        echo "$filename" 
        # # Resample to 3mm
        # flirt -in $file_path -ref masks/MNI152_T1_3mm_brain.nii.gz -out data/dataset_spreng/func/proc1_resample/$filename -applyisoxfm 3 \
        # -usesqform
        # # Mask and smooth
        # fslmaths data/dataset_spreng/func/proc1_resample/$filename -mul $mask -kernel gauss 2.123 \
        # -fmean data/dataset_spreng/func/proc2_mask_smooth/$filename
        # band pass filter (0.01-0.1Hz)
        python -m utils.signal.norm_filter -f data/dataset_spreng/func/proc2_mask_smooth/${filename}.gz -m $mask \
         -ch 0.1 -cl 0.01 -t 3 -o data/dataset_spreng/func/proc3_bandpass/${filename}
    done

    # mkdir -p data/dataset_spreng/physio/proc1_physio
    # echo "Physio preprocessing..."
    # for file_path in data/dataset_spreng/physio/raw/*.tsv.gz; do
    #     filename=$(basename $file_path)
    #     echo "$filename" 
    #     subj_file=$(cut -d'_' -f1 <<< "${filename}")
    #     filename_base=$(cut -d'.' -f1 <<< "${filename}")
    #     # HR and RV extraction
    #     python -m utils.dataset.preprocess_spreng -s $file_path -o data/dataset_spreng/physio/proc1_physio/$filename_base
    #     # Extract global BOLD signal from smoothed functional data
    #     # fslmeants -i data/dataset_spreng/func/proc6_standard/${subj_file}_ses-1_task-rest \
    #     # -o data/dataset_spreng/physio/proc1_physio/${filename_base}_global_sig.txt \
    #     # -m masks/MNI152_T1_3mm_gray_mask.nii.gz 
    #     # # CSF extraction
    #     # fslmeants -i data/dataset_spreng/func/proc6_standard/${subj_file}_ses-1_task-rest \
    #     # -o data/dataset_spreng/physio/proc1_physio/${filename_base}_csf.txt \
    #     # -m data/dataset_spreng/anat/proc4_csfmask/group_csfmask
    # done

fi 

# Yale - FMRI, Pupillometry
if [ "$dataset" == "yale" ]; then
    # mkdir -p data/dataset_yale/anat/proc1_bet
    # mkdir -p data/dataset_yale/anat/proc2_affine
    # mkdir -p data/dataset_yale/anat/proc3_fnirt
    # mkdir -p data/dataset_yale/anat/proc4_csfmask


    # echo "Structural preprocessing..."
    # # Structural preprocessing
    # for file_path in data/dataset_yale/anat/raw/*.nii.gz; do
    #     filename=$(basename $file_path)
    #     filename_base=$(cut -d'.' -f1 <<< "${filename}")
    #     echo "$filename"
    #     # Brain Extraction
    #     bet $file_path data/dataset_yale/anat/proc1_bet/$filename -o -m
    #     # FAST segmentation
    #     fast data/dataset_yale/anat/proc1_bet/$filename_base
    #     # Affine registration
    #     flirt -in data/dataset_yale/anat/proc1_bet/$filename -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz \
    #     -out data/dataset_yale/anat/proc2_affine/$filename -omat data/dataset_yale/anat/proc2_affine/$filename.mat
    #     # Nonlinear transformation
    #     fnirt --ref=$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz --in=$file_path \
    #     --iout=data/dataset_yale/anat/proc3_fnirt/$filename --cout=data/dataset_yale/anat/proc3_fnirt/$filename.mat \
    #     --aff=data/dataset_yale/anat/proc2_affine/$filename.mat --config=T1_2_MNI152_2mm --warpres=6,6,6
    #     Send CSF mask to MNI space
    #     applywarp --ref=masks/MNI152_T1_3mm_brain.nii.gz --in=data/dataset_yale/anat/proc1_bet/${filename_base}_pve_0  \
    #     --out=data/dataset_yale/anat/proc4_csfmask/$filename_base --warp="data/dataset_yale/anat/proc3_fnirt/${filename}.mat.nii.gz" 
    #     # Binarize CSF mask
    #     fslmaths data/dataset_yale/anat/proc4_csfmask/$filename_base -thr 0.9 -bin \
    #     data/dataset_yale/anat/proc4_csfmask/$filename_base
    # done

    # # # Create group CSF mask
    # ## remove mask if exists
    # rm -f data/dataset_yale/anat/proc4_csfmask/group_csfmask.nii.gz
    # ## Create empty mask
    # tmp_files=(data/dataset_spreng/anat/proc4_csfmask/*)    
    # fslmaths ${tmp_files[0]} -thr 2 data/dataset_yale/anat/proc4_csfmask/group_csfmask
    # ## Loop through subjects and add masks
    # for subj_mask in data/dataset_yale/anat/proc4_csfmask/*_T1w.nii.gz; do
    #     fslmaths data/dataset_yale/anat/proc4_csfmask/group_csfmask \
    #     -add $subj_mask data/dataset_yale/anat/proc4_csfmask/group_csfmask
    # done

    # # Threshold to >25 subjects overlap in CSF masks
    # fslmaths data/dataset_yale/anat/proc4_csfmask/group_csfmask -thr 25 -bin \
    # data/dataset_yale/anat/proc4_csfmask/group_csfmask

    # # Mask by MNI prior probability CSF mask (thresholded)
    # fslmaths data/dataset_yale/anat/proc4_csfmask/group_csfmask -mul masks/MNI152_T1_3mm_csf_mask \
    # data/dataset_yale/anat/proc4_csfmask/group_csfmask

    mkdir -p data/dataset_yale/func/procA_trim
    mkdir -p data/dataset_yale/func/proc1_mcflirt
    mkdir -p data/dataset_yale/func/proc2A_firstvolume
    mkdir -p data/dataset_yale/func/proc2B_func2struct
    mkdir -p data/dataset_yale/func/proc3_standard
    mkdir -p data/dataset_yale/func/proc4_mask_smooth
    mkdir -p data/dataset_yale/func/proc5_bandpass
    
    echo "Functional preprocessing..."
    for file_path in data/dataset_yale/func/raw/*.nii.gz; do
        filename=$(basename $file_path)
        echo "$filename" 
        # get base subject name to specify path to structural scan
        subj_file=$(cut -d'_' -f1 <<< "${filename}")
        filename_base=$(cut -d'.' -f1 <<< "${filename}")

        # Trim first ten volumes (index starts at 0)
        python -m utils.signal.trim -f $file_path -o data/dataset_yale/func/procA_trim/$filename -n 10
        # Motion correction (re-alignment)
        mcflirt -in data/dataset_yale/func/procA_trim/$filename -out data/dataset_yale/func/proc1_mcflirt/$filename -plots
        # Select first volume from each functional
        fslroi data/dataset_yale/func/proc1_mcflirt/$filename data/dataset_yale/func/proc2A_firstvolume/$filename 0 1
        # Co-registration with structural
        epi_reg --epi=data/dataset_yale/func/proc2A_firstvolume/$filename \
        --t1="data/dataset_yale/anat/raw/${subj_file}_T1w" --t1brain="data/dataset_yale/anat/proc1_bet/${subj_file}_T1w" \
        --out=data/dataset_yale/func/proc2B_func2struct/$filename
        # Get transform file to send functional to MNI
        applywarp --ref=masks/MNI152_T1_3mm_brain.nii.gz --in=data/dataset_yale/func/proc1_mcflirt/$filename \
        --out=data/dataset_yale/func/proc3_standard/$filename --warp="data/dataset_yale/anat/proc3_fnirt/${subj_file}_T1w.nii.gz.mat.nii.gz" \
        --premat="data/dataset_yale/func/proc2B_func2struct/${filename_base}.mat" 
        # Mask
        fslmaths data/dataset_yale/func/proc3_standard/$filename -mul $mask -kernel gauss 2.123 \
        -fmean data/dataset_yale/func/proc4_mask_smooth/$filename
        # bandpass filter (0.1Hz) and norm
        python -m utils.signal.norm_filter -f data/dataset_yale/func/proc4_mask_smooth/$filename -m $mask \
         -ch 0.1 -cl 0.01 -t 1 -o data/dataset_yale/func/proc5_bandpass/$filename

    done

    # mkdir -p data/dataset_yale/physio/proc1_physio
    # echo "Physio preprocessing..."
    # for file_path in data/dataset_yale/func/proc3_standard/*.nii.gz; do
    #     filename=$(basename $file_path)
    #     echo "$filename" 
    #     # get base subject name to specify path to structural scan
    #     subj_file=$(cut -d'_' -f1 <<< "${filename}")
    #     filename_base=$(cut -d'.' -f1 <<< "${filename}")
    #     # CSF extraction
    #     fslmeants -i $file_path -o data/dataset_yale/physio/proc1_physio/${filename_base}_csf.txt \
    #     -m data/dataset_yale/anat/proc4_csfmask/group_csfmask
    # done
fi



