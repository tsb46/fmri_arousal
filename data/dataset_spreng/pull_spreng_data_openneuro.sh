#!/bin/bash

mkdir -p func/raw
mkdir -p anat/raw
mkdir -p physio/raw

sed -n '2,$p' subject_list_spreng.csv | while IFS=, read -r subj other_cols;
do 
	subj_n=$(cut -d'-' -f2 <<< "${subj}")
    if [[ $subj_n -lt 10 ]]
        then
        	subj_new="sub-000${subj_n}"
    elif [[ $subj_n -lt 100 ]]
		then
			subj_new="sub-00${subj_n}"
    else
          subj_new="sub-0${subj_n}"
    fi

    echo $subj_new
    # Pull all session data
    aws s3 sync --no-sign-request s3://openneuro.org/ds003592 tmp/ --exclude "*" \
    --include "*${subj}/ses-1/*"

    # Rename and move functional echos
    mv "tmp/${subj}/ses-1/func/${subj}_ses-1_task-rest_echo-1_bold.nii.gz" "func/raw/${subj_new}_ses-1_task-rest_echo-1_bold.nii.gz"
    mv "tmp/${subj}/ses-1/func/${subj}_ses-1_task-rest_echo-1_bold.json" "func/raw/${subj_new}_ses-1_task-rest_echo-1_bold.json"
    mv "tmp/${subj}/ses-1/func/${subj}_ses-1_task-rest_echo-2_bold.nii.gz" "func/raw/${subj_new}_ses-1_task-rest_echo-2_bold.nii.gz"
    mv "tmp/${subj}/ses-1/func/${subj}_ses-1_task-rest_echo-2_bold.json" "func/raw/${subj_new}_ses-1_task-rest_echo-2_bold.json"
    mv "tmp/${subj}/ses-1/func/${subj}_ses-1_task-rest_echo-3_bold.nii.gz" "func/raw/${subj_new}_ses-1_task-rest_echo-3_bold.nii.gz"
    mv "tmp/${subj}/ses-1/func/${subj}_ses-1_task-rest_echo-3_bold.json" "func/raw/${subj_new}_ses-1_task-rest_echo-3_bold.json"

    # Rename and move physio data
    mv "tmp/${subj}/ses-1/func/${subj}_ses-1_task-rest_physio.tsv.gz" "physio/raw/${subj_new}_ses-1_task-rest_physio.tsv.gz"
    mv "tmp/${subj}/ses-1/func/${subj}_ses-1_task-rest_physio.json" "physio/raw/${subj_new}_ses-1_task-rest_physio.json"

    # Rename and move T1w data
    mv "tmp/${subj}/ses-1/anat/${subj}_ses-1_T1w.nii.gz" "anat/raw/${subj_new}_ses-1_T1w.nii.gz"

    # Remove tmp dir
    rm -r tmp

done

