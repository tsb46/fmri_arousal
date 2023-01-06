#!/bin/bash
# This is a runner script to replicate all main results in the study
# Note: This script must run to completion for code to run in the analysis.ipynb 

# Run PCA analyses across datasets
mkdir -p results/pca 

# ME-REST
python run_pca.py -n 10 -d chang 
mv chang_pca_group.nii results/pca/chang_pca_group.nii
mv chang_pca_group_results.pkl results/pca/chang_pca_group_results.pkl

# HCP-REST
python run_pca.py -n 10 -d hcp 
mv hcp_pca_group.nii results/pca/hcp_pca_group.nii
mv hcp_pca_group_results.pkl results/pca/hcp_pca_group_results.pkl

# ME-REST-SUPP
python run_pca.py -n 10 -d spreng
mv spreng_pca_group.nii results/pca/spreng_pca_group.nii
mv spreng_pca_group_results.pkl results/pca/spreng_pca_group_results.pkl

# ME-TASK
python run_pca.py -n 10 -d chang_bh
mv chang_bh_pca_group.nii results/pca/chang_bh_pca_group.nii
mv chang_bh_pca_group_results.pkl results/pca/chang_bh_pca_group_results.pkl

# NKI-TASK
python run_pca.py -n 10 -d nki
mv nki_pca_group.nii results/pca/nki_pca_group.nii
mv nki_pca_group_results.pkl results/pca/nki_pca_group_results.pkl

# YALE-REST
python run_pca.py -n 10 -d yale
mv yale_pca_group.nii results/pca/yale_pca_group.nii
mv yale_pca_group_results.pkl results/pca/yale_pca_group_results.pkl


# Run Complex PCA analyses across datasets
mkdir -p results/cpca 

# ME-REST
python run_pca.py -n 10 -d chang -t complex
mv chang_pca* results/cpca

# HCP-REST
python run_pca.py -n 10 -d hcp -t complex
mv hcp_pca* results/cpca

# ME-REST-SUPP
python run_pca.py -n 10 -d spreng -t complex
mv spreng_pca* results/cpca

# ME-TASK
python run_pca.py -n 10 -d chang_bh -t complex
mv chang_bh_pca* results/cpca

# NKI-TASK
python run_pca.py -n 10 -d nki -t complex
mv nki_pca* results/cpca

# YALE-REST
python run_pca.py -n 10 -d yale -t complex
mv yale_pca* results/cpca


# Run Distributed-Lag Non-linear Model (DLNM) across datasets
mkdir -p results/physio_dlnm

# ME-REST
python run_physio_dlnm.py -d chang -p rv -pl 15 -nl 5 
mv chang_dlnm* results/physio_dlnm

# HCP-REST
python run_physio_dlnm.py -d hcp -p rv -pl 45 -nl 15 
mv hcp_dlnm* results/physio_dlnm

# ME-REST-SUPP
python run_physio_dlnm.py -d spreng -p rv -pl 10 -nl 4 
mv spreng_dlnm* results/physio_dlnm


# Get Cross-Correlation Maps with Low-Frequency PPG Signals (Figure 2)
mkdir -p results/physio_crosscorr

# ME-REST
python run_crosscorr.py -d chang -p ppg_low -c 20
mv chang_cc* results/physio_crosscorr 

# HCP-REST
python run_crosscorr.py -d hcp -p ppg_low -c 30
mv hcp_cc* results/physio_crosscorr 





