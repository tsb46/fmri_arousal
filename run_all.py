import json
import numpy as np
import os

from run_pca import run_pca
from run_corr import run_crosscorr

# load analysis parameter json
params = json.load(open('analysis_params.json', 'rb'))


# run pca and cpca across all datasets (n_comps = 10)
# specify datasets for PCA
datasets_pca = [
	'chang', 
	'chang_bh', 
	'chang_cue', 
	'natview',
	'nki', 
	'hcp', 
	'spreng', 
	'yale',
	'toronto'
]
os.makedirs('results/pca', exist_ok=True)
for d in datasets_pca:
	print(f'run pca for {d}')
	# run PCA
	run_pca(d, n_comps=10, pca_type='real', rotate=None, 
	        recon=False, regress_global=False, 
	        out_dir='results/pca')
	# run CPCA
	print(f'run cpca for {d}')
	run_pca(
		d, n_comps=10, pca_type='complex', rotate=None, 
		recon=True, regress_global=False, out_dir='results/pca'
	)


# run cross correlation with PPG for ME-REST-SUPP dataset (& multiecho t2* and s0 signals)
os.makedirs('results/crosscorr', exist_ok=True)
print(f'run crosscorr for spreng dataset - t2')
run_crosscorr(
	'spreng', physio='PPG_PEAK_AMP', max_lag=10,
	n_samples=40, m_param='t2', 
	out_dir='results/crosscorr'
)
# rename to avoid naming conflict
os.rename(
	'results/crosscorr/spreng_cc_PPG_PEAK_AMP.nii', 
	'results/crosscorr/spreng_cc_t2_PPG_PEAK_AMP.nii'
)
print(f'run crosscorr for spreng dataset - s0')
run_crosscorr(
	'spreng', physio='PPG_PEAK_AMP', max_lag=10,
	n_samples=40, m_param='s0', 
	out_dir='results/crosscorr'
)
# rename to avoid naming conflict
os.rename(
	'results/crosscorr/spreng_cc_PPG_PEAK_AMP.nii', 
	'results/crosscorr/spreng_cc_s0_PPG_PEAK_AMP.nii'
)
print(f'run crosscorr for spreng dataset')
run_crosscorr(
	'spreng', physio='PPG_PEAK_AMP', max_lag=10,
	n_samples=40, m_param=None, 
	out_dir='results/crosscorr'
)







