#!/bin/bash
dataset=$1

# Dilated mask that includes sinuses slightly outside gray matter tissue
mask="masks/MNI152_T1_3mm_brain_mask_dilated2.nii.gz"

# Chang - EEG, FMRI, physio data
if [ "$dataset" == "chang" ]; then
    # mkdir -p data/dataset_chang/anat/proc1_reorient
    # mkdir -p data/dataset_chang/anat/proc2_bet
    # mkdir -p data/dataset_chang/anat/proc3_affine
    # mkdir -p data/dataset_chang/anat/proc4_fnirt
    # echo "Structural preprocessing..."
    # for file_path in data/dataset_chang/anat/raw/*.nii; do
    #     filename=$(basename $file_path)
    #     echo "$filename"
    #     # Reorient to standard orientation
    #     fslreorient2std $file_path data/dataset_chang/anat/proc1_reorient/$filename
    #     # Brain Extraction
    #     bet data/dataset_chang/anat/proc1_reorient/$filename data/dataset_chang/anat/proc2_bet/$filename -o -m
    #     # Affine registration
    #     flirt -in data/dataset_chang/anat/proc2_bet/$filename -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz \
    #     -out data/dataset_chang/anat/proc3_affine/$filename -omat data/dataset_chang/anat/proc3_affine/$filename.mat
    #     # Nonlinear transformation
    #     fnirt --ref=$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz --in=data/dataset_chang/anat/proc1_reorient/$filename \
    #     --iout=data/dataset_chang/anat/proc4_fnirt/$filename --cout=data/dataset_chang/anat/proc4_fnirt/$filename.mat \
    #     --aff=data/dataset_chang/anat/proc3_affine/$filename.mat --config=T1_2_MNI152_2mm --warpres=6,6,6
    # done


    # mkdir -p data/dataset_chang/func/proc1A_firstvolume
    # mkdir -p data/dataset_chang/func/proc1B_func2struct
    # mkdir -p data/dataset_chang/func/proc2_standard
    mkdir -p data/dataset_chang/func/proc3_mask_smooth
    mkdir -p data/dataset_chang/func/proc4_filter_norm
    mkdir -p data/dataset_chang/func/proc5_bandpass


    echo "Functional preprocessing..."
    for file_path in data/dataset_chang/func/raw/*.nii; do
        filename=$(basename $file_path)
        echo "$filename" 
        # get base subject name to specify path to structural scan
        subj_file=$(cut -d'-' -f1 <<< "${filename}")
        sess_n=$(cut -d'-' -f2 <<< "${filename}")
        filename_base=$(cut -d'.' -f1 <<< "${filename}")
        subj_out=${subj_file}_${sess_n}
        # # Select first volume from each functional
        # fslroi $file_path data/dataset_chang/func/proc1A_firstvolume/$subj_out 0 1
        # Co-registration with structural
        # epi_reg --epi=data/dataset_chang/func/proc1A_firstvolume/$subj_out \
        # --t1="data/dataset_chang/anat/proc1_reorient/${subj_file}-anat" \
        # --t1brain="data/dataset_chang/anat/proc2_bet/${subj_file}-anat" \
        # --out=data/dataset_chang/func/proc1B_func2struct/$subj_out
        # # Get transform file to send functional to MNI
        # applywarp --ref=masks/MNI152_T1_3mm_brain.nii.gz --in=$file_path \
        # --out=data/dataset_chang/func/proc2_standard/${subj_out} \
        # --warp="data/dataset_chang/anat/proc4_fnirt/${subj_file}-anat.nii.mat.nii.gz" \
        # --premat="data/dataset_chang/func/proc1B_func2struct/${subj_out}.mat" 
        # Mask
        fslmaths data/dataset_chang/func/proc2_standard/${subj_out} -mul $mask -kernel gauss 2.123 \
        -fmean data/dataset_chang/func/proc3_mask_smooth/${subj_out}
        # low pass filter (0.1Hz) and norm
        python -m utils.signal.norm_filter -f data/dataset_chang/func/proc3_mask_smooth/${subj_out}.nii.gz -m $mask \
         -ch 0.1 -t 1 -o data/dataset_chang/func/proc4_filter_norm/${subj_out}.nii.gz
        # band pass filter (0.01-0.1Hz) and norm for viz purposes
        python -m utils.signal.norm_filter -f data/dataset_chang/func/proc3_mask_smooth/${subj_out}.nii.gz -m $mask \
         -ch 0.1 -cl 0.01 -t 1 -o data/dataset_chang/func/proc5_bandpass/${subj_out}.nii.gz

    done

    # mkdir -p data/dataset_chang/physio_csf/raw
    # mkdir -p data/dataset_chang/physio_csf/proc1_filt
    # echo "Physio preprocessing..."
    # # PPG and Resp time courses were received preprocessed - only CSF inflow extraction
    # # This code assumes that csf inflow masks have been manually created and placed in a directory "data/dataset_chang/func/inflow_mask"
    # for file_path in data/dataset_chang/func/raw/*.nii; do
    #     filename=$(basename $file_path)
    #     echo "$filename" 
    #     # get base subject name to specify path to mask
    #     subj_file=$(cut -d'-' -f1 <<< "${filename}")
    #     sess_n=$(cut -d'-' -f2 <<< "${filename}")
    #     filename_base=$(cut -d'.' -f1 <<< "${filename}")
    #     subj_out=${subj_file}_${sess_n}
    #     # Extract CSF inflow signals
    #     fslmeants -i $file_path -o data/dataset_chang/physio_csf/raw/${subj_out}.txt \
    #     -m data/dataset_chang/func/inflow_mask/${filename_base}_mask
        
    # done
fi


# Yale - FMRI, Pupillometry
if [ "$dataset" == "yale" ]; then
    mkdir -p data/dataset_yale/anat/proc1_bet
    mkdir -p data/dataset_yale/anat/proc2_affine
    mkdir -p data/dataset_yale/anat/proc3_fnirt

    echo "Structural preprocessing..."
    Structural preprocessing
    for file_path in data/dataset_yale/anat/raw/*.nii.gz; do
        filename=$(basename $file_path)
        echo "$filename"
        # Brain Extraction
        bet $file_path data/dataset_yale/anat/proc1_bet/$filename -o -m
        # Affine registration
        flirt -in data/dataset_yale/anat/proc1_bet/$filename -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz \
        -out data/dataset_yale/anat/proc2_affine/$filename -omat data/dataset_yale/anat/proc2_affine/$filename.mat
        # Nonlinear transformation
        fnirt --ref=$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz --in=$file_path \
        --iout=data/dataset_yale/anat/proc3_fnirt/$filename --cout=data/dataset_yale/anat/proc3_fnirt/$filename.mat \
        --aff=data/dataset_yale/anat/proc2_affine/$filename.mat --config=T1_2_MNI152_2mm --warpres=6,6,6
    done

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
        # Motion correction (re-alignment)
        mcflirt -in $file_path -out data/dataset_yale/func/proc1_mcflirt/$filename -plots
        Select first volume from each functional
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
        # low pass filter (0.1Hz) and norm
        python -m utils.signal.norm_filter -f data/dataset_yale/func/proc4_mask_smooth/$filename -m $mask \
         -ch 0.1 -t 1 -o data/dataset_yale/func/proc5_filter_norm/$filename
        # Trim first ten volumes
        python -m utils.signal.trim -f data/dataset_yale/func/proc5_filter_norm/$filename -o data/dataset_yale/func/proc6_trim/$filename -n 10


    done
fi


# NKI - FMRI Breath-hold task
if [ "$dataset" == "nki" ]; then
    # mkdir -p data/dataset_nki/anat/proc1_bet
    # mkdir -p data/dataset_nki/anat/proc2_affine
    # mkdir -p data/dataset_nki/anat/proc3_fnirt

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
        # python -m utils.preprocess_resp_ppg -d nki -f $file_path -r 2 -p 1 -g 3 -t 1.4 -s 62.5 \
        # -n 186 -o data/dataset_nki/physio/proc1_physio/$filename_base
        # CSF extraction
        fslmeants -i data/dataset_nki/func/raw/${subj_file}_task_breathhold \
        -o data/dataset_nki/physio/proc1_physio/${subj_file}task_breathhold_physio_csf.txt \
        -m data/dataset_nki/func/inflow_mask/${subj_file}_task_breathhold_mask
    done
fi

# HCP - FMRI, hr, rv
if [ "$dataset" == "hcp" ]; then
    # Make directories
    mkdir -p data/dataset_hcp/func/proc1_resample
    mkdir -p data/dataset_hcp/func/proc2_mask_smooth
    mkdir -p data/dataset_hcp/func/proc3_filter_norm

    echo "Functional preprocessing..."
    for file_path in data/dataset_hcp/func/raw/*.nii.gz; do
        filename=$(basename $file_path)
        echo "$filename" 

        # # Resample to 3mm
        # flirt -in $file_path -ref $mask -out data/dataset_hcp/func/proc1_resample/$filename -applyisoxfm 3
        # # Mask and smooth
        fslmaths data/dataset_hcp/func/proc1_resample/$filename -mul $mask -kernel gauss 2.123 \
        -fmean data/dataset_hcp/func/proc2_mask_smooth/$filename
        # Lowpass filter (0.1Hz) and norm
        python -m utils.signal.norm_filter -f data/dataset_hcp/func/proc2_mask_smooth/$filename -m $mask \
         -ch 0.1 -t 2.1 -o data/dataset_hcp/func/proc3_filter_norm/$filename
    done

    # mkdir -p data/dataset_hcp/physio/proc1_hr_rv
    # echo "Physio preprocessing..."
    # for file_path in data/dataset_hcp/physio/raw/*.txt; do
    #     filename=$(basename $file_path)
    #     echo "$filename" 
    #     # get base subject name to specify path to structural scan
    #     subj_file=$(cut -d'_' -f1 <<< "${filename}")
    #     filename_base=$(cut -d'.' -f1 <<< "${filename}")
    #     # HR and RV extraction
    #     python -m utils.preprocess_resp_ppg -d hcp -f $file_path -r 1 -p 2 -g 0 -t 0.72 -s 400 \
    #     -n 1200 -o data/dataset_hcp/physio/proc1_hr_rv/$filename_base
    # done
fi


# DREAM - EEG, physio
if [ "$dataset" == "dream" ]; then
    # Make directories
    mkdir -p data/dataset_dream/eeg/proc1_fbands
    mkdir -p data/dataset_dream/physio

    echo "EEG and physio preprocessing..."

    for file_path in data/dataset_dream/eeg/raw/*.edf; do
        filename=$(basename $file_path)
        echo "$filename" 
        # get base subject name to specify path to structural scan
        filename_base=$(cut -d'.' -f1 <<< "${filename}")
        python -m utils.dataset.preprocess_dream -p "data/dataset_dream/eeg/raw/${filename_base}.edf" \
         -g "data/dataset_dream/hypnogram_rk/HypnogramR&K_${filename_base}.txt" \
        -oe data/dataset_dream/eeg/proc1_fbands/${filename_base}.csv \
        -op data/dataset_dream/physio/${filename_base}.csv
    done 
    
fi

# EKE - TCD, fNIRS, physio
if [ "$dataset" == "eke" ]; then
    # Make directories
    mkdir -p data/dataset_eke/data/preprocessed
    echo "Physio preprocessing..."
    sed 1d data/dataset_eke/subject_list_eke.csv | while IFS=, read -r subject extra
    do
        echo $subject
        # Preprocess EEG and physio
        python -m utils.dataset.preprocess_eke -d "data/dataset_eke/data/orig/subj0${subject}.csv" \
        -c "data/dataset_eke/data/orig/subj0${subject}.marker.csv" \
        -o "data/dataset_eke/data/preprocessed/subj0${subject}.csv"
    done 
fi

# LEMON - fMRI, hr, rv, bpp
if [ "$dataset" == "lemon" ]; then

    # mkdir -p data/dataset_lemon/anat/proc1A_masks
    # mkdir -p data/dataset_lemon/anat/proc1B_masked
    # mkdir -p data/dataset_lemon/anat/proc2_crop
    # mkdir -p data/dataset_lemon/anat/proc3_bet
    # mkdir -p data/dataset_lemon/anat/proc4_affine
    # mkdir -p data/dataset_lemon/anat/proc5_fnirt

    # echo "Structural preprocessing..."
    # # Structural preprocessing
    # for file_path in data/dataset_lemon/anat/raw/*T1w.nii.gz; do
    #     filename=$(basename $file_path)
    #     # get base subject name to specify path to structural scan
    #     subj_file=$(cut -d'_' -f1 <<< "${filename}")
    #     echo "$subj_file"
    #     # Mask inverse mprage image
    #     fslmaths data/dataset_lemon/anat/raw/${subj_file}_inv2 -thrP 5 -bin -fillh data/dataset_lemon/anat/proc1A_masks/${subj_file}_mask
    #     # Apply mask to T1w images
    #     fslmaths data/dataset_lemon/anat/raw/$filename -mul data/dataset_lemon/anat/proc1A_masks/${subj_file}_mask \
    #     data/dataset_lemon/anat/proc1B_masked/$filename
    #     # Robust fov
    #     robustfov -i data/dataset_lemon/anat/proc1B_masked/$filename -r data/dataset_lemon/anat/proc2_crop/$filename 
    #     # Brain Extraction of inverse mprage image
    #     bet data/dataset_lemon/anat/proc2_crop/$filename data/dataset_lemon/anat/proc3_bet/$filename -o -m -f 0.25 -B
    #     # Affine registration
    #     flirt -in data/dataset_lemon/anat/proc3_bet/$filename -ref $FSLDIR/data/standard/MNI152_T1_2mm_brain.nii.gz \
    #     -out data/dataset_lemon/anat/proc4_affine/$filename -omat data/dataset_lemon/anat/proc4_affine/$filename.mat
    #     # Nonlinear transformation
    #     fnirt --ref=$FSLDIR/data/standard/MNI152_T1_2mm.nii.gz --in=data/dataset_lemon/anat/proc2_crop/$filename \
    #     --iout=data/dataset_lemon/anat/proc5_fnirt/$filename --cout=data/dataset_lemon/anat/proc5_fnirt/$filename.mat \
    #     --aff=data/dataset_lemon/anat/proc4_affine/$filename.mat --config=T1_2_MNI152_2mm --warpres=6,6,6
    # done

    # mkdir -p data/dataset_lemon/func/proc1_mcflirt
    # mkdir -p data/dataset_lemon/func/proc2A_firstvolume
    # mkdir -p data/dataset_lemon/func/proc2B_func2struct
    # mkdir -p data/dataset_lemon/func/proc3_standard
    # mkdir -p data/dataset_lemon/func/proc4_mask_smooth
    # mkdir -p data/dataset_lemon/func/proc5_filter_norm
    # mkdir -p data/dataset_lemon/func/proc6_trim

    # echo "Functional preprocessing..."
    # for file_path in data/dataset_lemon/func/raw/*.nii.gz; do
    #     filename=$(basename $file_path)
    #     echo "$filename" 
    #     # get base subject name to specify path to structural scan
    #     subj_file=$(cut -d'_' -f1 <<< "${filename}")
    #     filename_base=$(cut -d'.' -f1 <<< "${filename}")
        # # Motion correction (re-alignment)
        # mcflirt -in $file_path -out data/dataset_lemon/func/proc1_mcflirt/$filename -plots
        # # Select first volume from each functional
        # fslroi data/dataset_lemon/func/proc1_mcflirt/$filename data/dataset_lemon/func/proc2A_firstvolume/$filename 0 1
        # # # Co-registration with structural
        # epi_reg --epi=data/dataset_lemon/func/proc2A_firstvolume/$filename \
        # --t1="data/dataset_lemon/anat/proc2_crop/${subj_file}_T1w" --t1brain="data/dataset_lemon/anat/proc3_bet/${subj_file}_T1w" \
        # --out=data/dataset_lemon/func/proc2B_func2struct/$filename
        # # Get transform file to send functional to MNI
        # applywarp --ref=masks/MNI152_T1_3mm_brain.nii.gz --in=data/dataset_lemon/func/proc1_mcflirt/$filename \
        # --out=data/dataset_lemon/func/proc3_standard/$filename --warp="data/dataset_lemon/anat/proc5_fnirt/${subj_file}_T1w.nii.gz.mat.nii.gz" \
        # --premat="data/dataset_lemon/func/proc2B_func2struct/${filename_base}.mat" 
        # # Mask
        # fslmaths data/dataset_lemon/func/proc3_standard/$filename -mul $mask -kernel gauss 2.123 \
        # -fmean data/dataset_lemon/func/proc4_mask_smooth/$filename
        # # low filter (0.1Hz) and norm
        # python -m utils.signal.norm_filter -f data/dataset_lemon/func/proc4_mask_smooth/$filename -m $mask \
        #  -ch 0.1 -t 1.4 -o data/dataset_lemon/func/proc5_filter_norm/$filename
        #  # trim time series (it seems PPG and BP only extend to 470 TRs)
        #  python -m utils.signal.trim -f data/dataset_lemon/func/proc5_filter_norm/$filename -n 0 -n_end 470 \
        #  -o data/dataset_lemon/func/proc6_trim/$filename
    # done

    mkdir -p data/dataset_lemon/physio/proc1_physio
    echo "Physio preprocessing..."
    for file_path in data/dataset_lemon/anat/raw/*_T1w.nii.gz; do
        filename=$(basename $file_path) 
        # get base subject name to specify path to physio filesd
        subj_file=$(cut -d'_' -f1 <<< "${filename}")
        echo "$subj_file"        
        # # preprocess physio files (func time samples = 657)a
        # python -m utils.dataset.preprocess_lemon -s $subj_file -o data/dataset_lemon/physio/proc1_physio/${subj_file}_physio
        # Extract CSF inflow signals
        fslmeants -i data/dataset_lemon/func/raw/${subj_file}_task_rest -o data/dataset_lemon/physio/proc1_physio/${subj_file}_physio_csf.txt \
        -m data/dataset_lemon/func/inflow_mask/${subj_file}_task_rest_mask
    done
    
fi

