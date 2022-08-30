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
    #     # Preprocess EEG and Physio Data
    #     python -m utils.dataset.preprocess_chang \
    #     -e data/dataset_chang_bh/eeg/raw/${subj_file}-${sess_n}-adb_echo1_EEG_pp.mat \
    #     -p data/dataset_chang_bh/physio/raw/${subj_file}-${sess_n}-adb_echo1_physOUT.mat \
    #      -f $nframes -om data/dataset_chang_bh/eeg/raw/${subj_out} \
    #      -oe data/dataset_chang_bh/eeg/proc1_fbands/${subj_out}_fbands \
    #      -op data/dataset_chang_bh/physio/proc1_physio/${subj_out}_physio  
        # Extract global BOLD signal from preprocessed low-pass functional data
        fslmeants -i data/dataset_chang_bh/func/proc3_filter_norm/${subj_file}-${sess_n}-adb_echo1_w_dspk_blur3mm \
        -o data/dataset_chang_bh/physio/proc1_physio/${subj_out}_global_sig.txt \
        -m masks/MNI152_T1_3mm_gray_mask.nii.gz
        # Extract CSF signal from raw functional data (post ME-ICA)
        fslmeants -i data/dataset_chang_bh/func/raw/${subj_file}-${sess_n}-adb_echo1_w_dspk_blur3mm  \
        -o data/dataset_chang_bh/physio/proc1_physio/${subj_out}_csf.txt \
        -m data/dataset_chang_bh/anat/proc5_csfmask/group_csf_mask
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
        # HR and RV extraction
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

    # Threshold to >40 subjects overlap in CSF masks
    # fslmaths data/dataset_hcp/anat/proc1_csfmask/group_csfmask -thr 40 -bin \
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

        # # Physio extraction
        # python -m utils.dataset.preprocess_hcp -s $file_path -o data/dataset_hcp/physio/proc1_physio/${subj_file}_physio -d rest
        # # Extract precuneus BOLD signal from preprocessed band-pass functional data
        # fslmeants -i data/dataset_hcp/func_fix/proc4_bandpass/${subj_func}_rest.nii.gz\
        # -o data/dataset_hcp/physio/proc1_physio/${subj_file}_precuneus.txt \
        # -m masks/precuneus_sphere_6mm.nii.gz
        # Extract global BOLD signal from preprocessed low-pass functional data
        fslmeants -i data/dataset_hcp/func_fix/proc3_filter_norm/${subj_func}_rest.nii.gz \
        -o data/dataset_hcp/physio/proc1_physio/${subj_file}_global_sig.txt \
        -m masks/MNI152_T1_3mm_gray_mask.nii.gz
        # # Extract superior parietal BOLD signal from preprocessed band-pass functional data
        # fslmeants -i data/dataset_hcp/func_fix/proc4_bandpass/${subj_func}_rest.nii.gz\
        # -o data/dataset_hcp/physio/proc1_physio/${subj_file}_superior_parietal.txt \
        # -m masks/superior_parietal_sphere_6mm.nii.gz
        
        # Extract CSF signal from HCP CSF mask
        fslmeants -i data/dataset_hcp/func_fix/raw/${subj_func}_rest.nii.gz\
        -o data/dataset_hcp/physio/proc1_physio/${subj_file}_csf.txt \
        -m data/dataset_hcp/anat/proc1_csfmask/group_csfmask

    done

fi

# HCP Relational Task - FMRI, hr, rv
if [ "$dataset" == "hcp_rel" ]; then

    # mkdir -p data/dataset_hcp_task/func_rel/proc1_resample
    # mkdir -p data/dataset_hcp_task/func_rel/proc2_mask_smooth
    # mkdir -p data/dataset_hcp_task/func_rel/proc3_filter_norm
    # mkdir -p data/dataset_hcp_task/func_rel/proc4_bandpass

    # echo "Functional preprocessing for minimally preprocessed data..."
    # for file_path in data/dataset_hcp_task/func_rel/raw/*.nii.gz; do
    #     filename=$(basename $file_path)
    #     echo "$filename" 
    #     # Resample to 3mm
    #     flirt -in $file_path -ref $mask -out data/dataset_hcp_task/func_rel/proc1_resample/$filename -applyisoxfm 3
    #     # Mask and smooth
    #     fslmaths data/dataset_hcp_task/func_rel/proc1_resample/$filename -mul $mask -kernel gauss 2.123 \
    #     -fmean data/dataset_hcp_task/func_rel/proc2_mask_smooth/$filename
    #     # Lowpass filter (0.1Hz) and norm
    #     python -m utils.signal.norm_filter -f data/dataset_hcp_task/func_rel/proc2_mask_smooth/$filename -m $mask \
    #      -ch 0.1 -t 0.72 -o data/dataset_hcp_task/func_rel/proc3_filter_norm/$filename
    #     # band pass filter (0.01-0.1Hz)
    #     python -m utils.signal.norm_filter -f data/dataset_hcp_task/func_rel/proc2_mask_smooth/$filename -m $mask \
    #      -ch 0.1 -cl 0.01 -t 0.72 -o data/dataset_hcp_task/func_rel/proc4_bandpass/$filename
    # done

    mkdir -p data/dataset_hcp_task/physio_rel/proc1_physio
    echo "Physio preprocessing..."
    for file_path in data/dataset_hcp_task/physio_rel/raw/*.txt; do
        filename=$(basename $file_path)
        echo "$filename" 
        # get base subject name to specify path to structural scan
        subj_file=$(cut -d'_' -f1 <<< "${filename}")
        subj_func=$(echo $filename | awk 'BEGIN {FS="_physio.txt" } ; { print $1 }')
        # Physio extraction
        python -m utils.dataset.preprocess_hcp -s $file_path -o data/dataset_hcp_task/physio_rel/proc1_physio/${subj_file}_physio -d rel
    done
fi

# HCP WM Task - FMRI, hr, rv
if [ "$dataset" == "hcp_wm" ]; then

    # mkdir -p data/dataset_hcp_task/func_wm/proc1_resample
    # mkdir -p data/dataset_hcp_task/func_wm/proc2_mask_smooth
    # mkdir -p data/dataset_hcp_task/func_wm/proc3_filter_norm
    # mkdir -p data/dataset_hcp_task/func_wm/proc4_bandpass

    # echo "Functional preprocessing for minimally preprocessed data..."
    # for file_path in data/dataset_hcp_task/func_wm/raw/*.nii.gz; do
    #     filename=$(basename $file_path)
    #     echo "$filename" 
    #     # Resample to 3mm
    #     flirt -in $file_path -ref $mask -out data/dataset_hcp_task/func_wm/proc1_resample/$filename -applyisoxfm 3
    #     # Mask and smooth
    #     fslmaths data/dataset_hcp_task/func_wm/proc1_resample/$filename -mul $mask -kernel gauss 2.123 \
    #     -fmean data/dataset_hcp_task/func_wm/proc2_mask_smooth/$filename
    #     # Lowpass filter (0.1Hz) and norm
    #     python -m utils.signal.norm_filter -f data/dataset_hcp_task/func_wm/proc2_mask_smooth/$filename -m $mask \
    #      -ch 0.1 -t 0.72 -o data/dataset_hcp_task/func_wm/proc3_filter_norm/$filename
    #     # band pass filter (0.01-0.1Hz)
    #     python -m utils.signal.norm_filter -f data/dataset_hcp_task/func_wm/proc2_mask_smooth/$filename -m $mask \
    #      -ch 0.1 -cl 0.01 -t 0.72 -o data/dataset_hcp_task/func_wm/proc4_bandpass/$filename
    # done

    mkdir -p data/dataset_hcp_task/physio_wm/proc1_physio
    echo "Physio preprocessing..."
    for file_path in data/dataset_hcp_task/physio_wm/raw/*.txt; do
        filename=$(basename $file_path)
        echo "$filename" 
        # get base subject name to specify path to structural scan
        subj_file=$(cut -d'_' -f1 <<< "${filename}")
        subj_func=$(echo $filename | awk 'BEGIN {FS="_physio.txt" } ; { print $1 }')
        # Physio extraction
        python -m utils.dataset.preprocess_hcp -s $file_path -o data/dataset_hcp_task/physio_wm/proc1_physio/${subj_file}_physio -d wm
    done
fi


# Monash - FMRI, PET
if [ "$dataset" == "monash" ]; then
    # mkdir -p data/dataset_monash/anat/proc1_bet
    # mkdir -p data/dataset_monash/anat/proc2_affine
    # mkdir -p data/dataset_monash/anat/proc3_fnirt

    # echo "Structural preprocessing..."
    # # Structural preprocessing
    # for file_path in data/dataset_monash/anat/raw/*.nii.gz; do
    #     filename=$(basename $file_path)
    #     echo "$filename"
    #     # Brain Extraction
    #     bet $file_path data/dataset_monash/anat/proc1_bet/$filename -o -m
    #     # Affine registration
    #     flirt -in data/dataset_monash/anat/proc1_bet/$filename -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz \
    #     -out data/dataset_monash/anat/proc2_affine/$filename -omat data/dataset_monash/anat/proc2_affine/$filename.mat
    #     # Nonlinear transformation
    #     fnirt --ref=$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz --in=$file_path \
    #     --iout=data/dataset_monash/anat/proc3_fnirt/$filename --cout=data/dataset_monash/anat/proc3_fnirt/$filename.mat \
    #     --aff=data/dataset_monash/anat/proc2_affine/$filename.mat --config=T1_2_MNI152_2mm --warpres=6,6,6
    # done

    mkdir -p data/dataset_monash/func/proc1_mcflirt
    mkdir -p data/dataset_monash/func/proc2_slicetime
    mkdir -p data/dataset_monash/func/proc3A_firstvolume
    mkdir -p data/dataset_monash/func/proc3B_aligned
    mkdir -p data/dataset_monash/func/proc4_func2struct
    mkdir -p data/dataset_monash/func/proc5_standard
    mkdir -p data/dataset_monash/func/proc6_mask_smooth
    mkdir -p data/dataset_monash/func/proc7_bandpass

    echo "Functional preprocessing..."
    sed -n '2,$p' data/dataset_monash/subject_list_monash.csv | while IFS=, read -r subj; do
        echo ${subj}
        # # Motion correction for all scans per subject
        # echo 'motion correct and slice time correct all scans'
        # for scan in {1..6}; do 
        #     echo $scan
        #     # Motion correction (re-alignment)
        #     mcflirt -in data/dataset_monash/func/raw/${subj}_run_${scan} -out data/dataset_monash/func/proc1_mcflirt/${subj}_run_${scan} -plots
        #     # Slice time correct with .txt file
        #     slicetimer -i data/dataset_monash/func/proc1_mcflirt/${subj}_run_${scan} \
        #     -o data/dataset_monash/func/proc2_slicetime/${subj}_run_${scan} --tcustom=data/dataset_monash/slicetime_monash.txt

        # done
        # echo 'align other scans to first functional scan'
        # # Align functional scans to the first 
        # echo 1 # scan 1
        # # Select first volume from each functional
        # fslroi data/dataset_monash/func/proc2_slicetime/${subj}_run_1 data/dataset_monash/func/proc3A_firstvolume/${subj}_run_1 0 1
        # for scan in {2..6}; do
        #     echo $scan
        #     # Select first volume
        #     fslroi data/dataset_monash/func/proc2_slicetime/${subj}_run_${scan} data/dataset_monash/func/proc3A_firstvolume/${subj}_run_${scan} 0 1
        #     # Align first volume to first volume of scan 1
        #     flirt -dof 6 -in data/dataset_monash/func/proc3A_firstvolume/${subj}_run_${scan} \
        #     -ref data/dataset_monash/func/proc3A_firstvolume/${subj}_run_1 \
        #     -omat data/dataset_monash/func/proc3A_firstvolume/${subj}_run_${scan}.mat
        #     # Align functional scan n to functional scan 1
        #     flirt -in data/dataset_monash/func/proc2_slicetime/${subj}_run_${scan} -ref data/dataset_monash/func/proc3A_firstvolume/${subj}_run_1 \
        #     -applyxfm -init data/dataset_monash/func/proc3A_firstvolume/${subj}_run_${scan}.mat \
        #     -out data/dataset_monash/func/proc3B_aligned/${subj}_run_${scan}
        # done

        
        # # Co-registration with structural
        # echo 'co-register first functional scan to structural T1'
        # epi_reg --epi=data/dataset_monash/func/proc3A_firstvolume/${subj}_run_1 \
        # --t1="data/dataset_monash/anat/raw/${subj}_T1w" --t1brain="data/dataset_monash/anat/proc1_bet/${subj}_T1w" \
        # --out=data/dataset_monash/func/proc4_func2struct/${subj}

        # echo 'functional to MNI space and smooth'
        # for scan in {1..6}; do 
        #     echo $scan
            # if [ $scan -eq 1 ]; then
            #     # Get transform file to send functional to MNI
            #     applywarp --ref=masks/MNI152_T1_3mm_brain.nii.gz --in=data/dataset_monash/func/proc2_slicetime/${subj}_run_1 \
            #     --out=data/dataset_monash/func/proc5_standard/${subj}_run_1 \
            #     --warp="data/dataset_monash/anat/proc3_fnirt/${subj}_T1w.nii.gz.mat.nii.gz" \
            #     --premat=data/dataset_monash/func/proc4_func2struct/${subj}.mat
            # else
            #     # Get transform file to send functional to MNI
            #     applywarp --ref=masks/MNI152_T1_3mm_brain.nii.gz --in=data/dataset_monash/func/proc3B_aligned/${subj}_run_${scan} \
            #     --out=data/dataset_monash/func/proc5_standard/${subj}_run_${scan} \
            #     --warp="data/dataset_monash/anat/proc3_fnirt/${subj}_T1w.nii.gz.mat.nii.gz" \
            #     --premat=data/dataset_monash/func/proc4_func2struct/${subj}.mat
            # fi

            # Mask
            # fslmaths data/dataset_monash/func/proc5_standard/${subj}_run_${scan} -mul masks/MNI152_T1_3mm_brain_mask \
            # -kernel gauss 2.123 -fmean data/dataset_monash/func/proc6_mask_smooth/${subj}_run_${scan} 
                
        # done

        # Concatenate scans together
        # python -m utils.signal.concatenate_func -f data/dataset_monash/func/proc6_mask_smooth/${subj}_run_1.nii.gz \
        # -f data/dataset_monash/func/proc6_mask_smooth/${subj}_run_2.nii.gz -f data/dataset_monash/func/proc6_mask_smooth/${subj}_run_3.nii.gz \
        # -f data/dataset_monash/func/proc6_mask_smooth/${subj}_run_4.nii.gz -f data/dataset_monash/func/proc6_mask_smooth/${subj}_run_5.nii.gz \
        # -f data/dataset_monash/func/proc6_mask_smooth/${subj}_run_6.nii.gz -o data/dataset_monash/func/proc6_mask_smooth/${subj}.nii.gz

        # # remove intermediate files
        # rm data/dataset_monash/func/proc6_mask_smooth/${subj}_run_1.nii.gz \
        # data/dataset_monash/func/proc2B_concat/${subj}_run_2.nii.gz data/dataset_monash/func/proc2B_concat/${subj}_run_3.nii.gz \
        # data/dataset_monash/func/proc2B_concat/${subj}_run_4.nii.gz data/dataset_monash/func/proc2B_concat/${subj}_run_5.nii.gz \
        # data/dataset_monash/func/proc2B_concat/${subj}_run_6.nii.gz
        
        
        # band pass filter (0.01-0.1Hz)
        python -m utils.signal.norm_filter -f data/dataset_monash/func/proc6_mask_smooth/${subj}.nii.gz -m masks/MNI152_T1_3mm_brain_mask.nii.gz \
         -ch 0.1 -cl 0.01 -t 2.45 -o data/dataset_monash/func/proc7_bandpass/${subj}.nii.gz
    done

    # mkdir -p data/dataset_monash/pet/raw_reduce
    # mkdir -p data/dataset_monash/pet/proc1_mcflirt
    # mkdir -p data/dataset_monash/pet/proc2A_firstvolume
    # mkdir -p data/dataset_monash/pet/proc2B_func2struct
    # mkdir -p data/dataset_monash/pet/proc3_standard
    # mkdir -p data/dataset_monash/pet/proc4_mask_smooth
    # mkdir -p data/dataset_monash/pet/proc5_detrend
    # mkdir -p data/dataset_monash/pet/proc6_convolve


    # echo "PET preprocessing..."
    # for file_path in data/dataset_monash/pet/raw/*.nii.gz; do
    #     filename=$(basename $file_path)
    #     # get base subject name to specify path to structural scan
    #     subj_file=$(cut -d'_' -f1 <<< "${filename}")
    #     filename_base=$(cut -d'.' -f1 <<< "${filename}")
    #     echo $subj_file
        # # Trim first 131 volumes (consistent with Jamadar et al. 2020; Scientific Data) to leave 225 volumes
        # fslroi data/dataset_monash/pet/raw/$filename data/dataset_monash/pet/raw_reduce/$filename  131 225
        # # Resample to 3mm isotropic for easier preprocessing and reduce FOV 
        # ## Reduce FOV
        # fslroi data/dataset_monash/pet/raw_reduce/$filename data/dataset_monash/pet/raw_reduce/$filename 85 260 80 260 -1 -1
        # ## Resample
        # flirt -in data/dataset_monash/pet/raw_reduce/$filename -ref data/dataset_monash/pet/raw_reduce/$filename \
        # -applyisoxfm 3.0 -nosearch -out data/dataset_monash/pet/raw_reduce/$filename
        # # Motion correction (re-alignment)
        # mcflirt -in data/dataset_monash/pet/raw_reduce/$filename -out data/dataset_monash/pet/proc1_mcflirt/$filename -plots
        # # Select first volume from pet images
        # fslroi data/dataset_monash/pet/proc1_mcflirt/$filename data/dataset_monash/pet/proc2A_firstvolume/$filename 0 1
        # # Co-registration with structural
        # epi_reg --epi=data/dataset_monash/pet/proc2A_firstvolume/$filename \
        # --t1="data/dataset_monash/anat/raw/${subj_file}_T1w" --t1brain="data/dataset_monash/anat/proc1_bet/${subj_file}_T1w" \
        # --out=data/dataset_monash/pet/proc2B_func2struct/$filename
        # # Get transform file to send pet to MNI
        # applywarp --ref=masks/MNI152_T1_3mm_brain.nii.gz --in=data/dataset_monash/pet/proc1_mcflirt/$filename \
        # --out=data/dataset_monash/pet/proc3_standard/$filename --warp="data/dataset_monash/anat/proc3_fnirt/${subj_file}_T1w.nii.gz.mat.nii.gz" \
        # --premat="data/dataset_monash/pet/proc2B_func2struct/${filename_base}.mat" 
        # # Mask
        # fslmaths data/dataset_monash/pet/proc3_standard/$filename -mul $mask -kernel gauss 2.123 \
        # -fmean data/dataset_monash/pet/proc4_mask_smooth/$filename
        # Polynomial Detrending (2nd order)
        # python -m utils.signal.poly_detrend -f data/dataset_monash/pet/proc4_mask_smooth/$filename -p 2 \
        # -m $mask -o data/dataset_monash/pet/proc5_detrend/$filename

    # done
fi


# Spreng - FMRI, Physio
if [ "$dataset" == "spreng" ]; then
    mkdir -p data/dataset_spreng/anat/raw_reorient
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
    # fslmaths data/dataset_nki/anat/proc4_csfmask/group_csfmask -mul masks/MNI152_T1_3mm_csf_mask \
    # data/dataset_nki/anat/proc4_csfmask/group_csfmask

    # echo "Functional preprocessing..."
    # mkdir -p data/dataset_spreng/func/proc1_mcflirt
    # mkdir -p data/dataset_spreng/func/proc2_slicetime
    # mkdir -p data/dataset_spreng/func/proc3_combination
    # mkdir -p data/dataset_spreng/func/proc3B_meica_combination
    # mkdir -p data/dataset_spreng/func/proc4_trim
    # mkdir -p data/dataset_spreng/func/proc5A_firstvolume
    # mkdir -p data/dataset_spreng/func/proc5B_func2struct
    # mkdir -p data/dataset_spreng/func/proc6_standard
    # mkdir -p data/dataset_spreng/func/proc7_mask_smooth
    # mkdir -p data/dataset_spreng/func/proc8_bandpass


    # sed -n '2,$p' data/dataset_spreng/subject_list_spreng.csv | while IFS=, read -r id subj other_cols; do
    #     echo ${subj}
    #     subj_fp=${subj}_ses-1_task-rest

        # # Motion correction for first echo
        # mcflirt -in data/dataset_spreng/func/raw/${subj_fp}_echo-1_bold \
        # -out data/dataset_spreng/func/proc1_mcflirt/${subj_fp}_echo-1_bold -mats -meanvol 
        # # Apply mcflirt transform params to second echo
        # applyxfm4D data/dataset_spreng/func/raw/${subj_fp}_echo-2_bold data/dataset_spreng/func/proc1_mcflirt/${subj_fp}_echo-1_bold_mean_reg \
        # data/dataset_spreng/func/proc1_mcflirt/${subj_fp}_echo-2_bold data/dataset_spreng/func/proc1_mcflirt/${subj_fp}_echo-1_bold.mat -fourdigit
        # # Apply mcflirt transform params to third echo
        # applyxfm4D data/dataset_spreng/func/raw/${subj_fp}_echo-3_bold data/dataset_spreng/func/proc1_mcflirt/${subj_fp}_echo-1_bold_mean_reg \
        # data/dataset_spreng/func/proc1_mcflirt/${subj_fp}_echo-3_bold data/dataset_spreng/func/proc1_mcflirt/${subj_fp}_echo-1_bold.mat -fourdigit
        # # Apply slicetiming correction to each echo
        # for echo in {1..3}; do 
        #     # Slice time correct with .txt file
        #     slicetimer -i data/dataset_spreng/func/proc1_mcflirt/${subj_fp}_echo-${echo}_bold \
        #     -o data/dataset_spreng/func/proc2_slicetime/${subj_fp}_echo-${echo}_bold --tcustom=data/dataset_spreng/slicetime_spreng.txt
        # done 

        # # Optimal echo combination through tedana
        # t2smap -d data/dataset_spreng/func/proc2_slicetime/${subj_fp}_echo-1_bold.nii.gz \
        # data/dataset_spreng/func/proc2_slicetime/${subj_fp}_echo-2_bold.nii.gz \
        # data/dataset_spreng/func/proc2_slicetime/${subj_fp}_echo-3_bold.nii.gz -e 13.7 30.0 47.0 \
        # --out-dir=data/dataset_spreng/func/proc3_combination/${subj_fp} --prefix=${subj_fp} \
        # --convention=orig  

        # # Multi-echo ICA and optimal combination through tedana
        # tedana -d data/dataset_spreng/func/proc2_slicetime/${subj_fp}_echo-1_bold.nii.gz \
        # data/dataset_spreng/func/proc2_slicetime/${subj_fp}_echo-2_bold.nii.gz \
        # data/dataset_spreng/func/proc2_slicetime/${subj_fp}_echo-3_bold.nii.gz -e 13.7 30.0 47.0 \
        # --out-dir=data/dataset_spreng/func/proc3B_meica_combination/${subj_fp} --prefix=${subj_fp} \
        # --convention=orig

    #     # Trim first 4 volumes 
    #     python -m utils.signal.trim -f data/dataset_spreng/func/proc3B_meica_combination/${subj_fp}/${subj_fp}_dn_ts_OC.nii.gz \
    #     -o data/dataset_spreng/func/proc4_trim/${subj_fp}.nii.gz -n 4

    #     # # Select first volume from each functional
    #     fslroi data/dataset_spreng/func/proc4_trim/${subj_fp} data/dataset_spreng/func/proc5A_firstvolume/${subj_fp} 0 1
    #     # Co-registration with structural
    #     epi_reg --epi=data/dataset_spreng/func/proc5A_firstvolume/${subj_fp} \
    #     --t1="data/dataset_spreng/anat/raw_reorient/${subj}_ses-1_T1w" --t1brain="data/dataset_spreng/anat/proc1_bet/${subj}_ses-1_T1w" \
    #     --out=data/dataset_spreng/func/proc5B_func2struct/${subj_fp}
    #     # Get transform file to send functional to MNI
    #     applywarp --ref=masks/MNI152_T1_3mm_brain.nii.gz --in=data/dataset_spreng/func/proc4_trim/${subj_fp} \
    #     --out=data/dataset_spreng/func/proc6_standard/${subj_fp} --warp="data/dataset_spreng/anat/proc3_fnirt/${subj}_ses-1_T1w.nii.gz.mat.nii.gz" \
    #     --premat="data/dataset_spreng/func/proc5B_func2struct/${subj_fp}.mat" 
    #     # Mask
    #     fslmaths data/dataset_spreng/func/proc6_standard/${subj_fp} -mul $mask -kernel gauss 2.123 \
    #     -fmean data/dataset_spreng/func/proc7_mask_smooth/${subj_fp}
    #     # bandpass filter (0.1Hz) and norm
    #     python -m utils.signal.norm_filter -f data/dataset_spreng/func/proc7_mask_smooth/${subj_fp}.nii.gz -m $mask \
    #     -ch 0.1 -cl 0.01 -t 3 -o data/dataset_spreng/func/proc8_bandpass/${subj_fp}.nii.gz

    # done  


    mkdir -p data/dataset_spreng/physio/proc1_physio
    echo "Physio preprocessing..."
    for file_path in data/dataset_spreng/physio/raw/*.tsv.gz; do
        filename=$(basename $file_path)
        echo "$filename" 
        subj_file=$(cut -d'_' -f1 <<< "${filename}")
        filename_base=$(cut -d'.' -f1 <<< "${filename}")
        # # HR and RV extraction
        # python -m utils.dataset.preprocess_spreng -s $file_path -o data/dataset_spreng/physio/proc1_physio/$filename_base
        # Extract global BOLD signal from smoothed functional data
        fslmeants -i data/dataset_spreng/func/proc6_standard/${subj_file}_ses-1_task-rest \
        -o data/dataset_spreng/physio/proc1_physio/${filename_base}_global_sig.txt \
        -m masks/MNI152_T1_3mm_gray_mask.nii.gz 
        # # CSF extraction
        # fslmeants -i data/dataset_spreng/func/proc6_standard/${subj_file}_ses-1_task-rest \
        # -o data/dataset_spreng/physio/proc1_physio/${filename_base}_csf.txt \
        # -m data/dataset_spreng/anat/proc4_csfmask/group_csfmask
    done

   

fi 



