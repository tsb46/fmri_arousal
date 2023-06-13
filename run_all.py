import json
import numpy as np
import os

from run_crosscorr import run_crosscorr
from run_pca import run_pca
from run_peak_average import run_peak_average
from run_task_avg import run_task_avg

# load analysis parameter json
params = json.load(open('analysis_params.json', 'rb'))

# specify datasets
datasets = [
	'chang', 
	'chang_bh', 
	'chang_cue', 
	'natview',
	'nki', 
	'hcp', 
	'spreng', 
	'yale'
]

# run pca and cpca across all datasets (n_comps = 10)
os.makedirs('results/pca', exist_ok=True)
for d in datasets:
	print(f'run pca for {d}')
	# run PCA
	run_pca(d, n_comps=10, pca_type='real', rotate=None, 
	        recon=False, regress_global=False, 
	        out_dir='results/pca')
	# run CPCA
	print(f'run cpca for {d}')
	run_pca(d, n_comps=10, pca_type='complex', rotate=None, 
	        recon=True, regress_global=False, out_dir='results/pca')

# run peak averaging across all datasets and physiological signals
# set left and right window around peaks as 15 and 15 secs, respectively
# set peak threshold as 2 z-score units
os.makedirs('results/peak_avg', exist_ok=True)
for d in datasets:
	print(f'run peak average for {d}')
	# get left and right 15s windows (in TRs)
	w_len = int(np.ceil(15/params[d]['tr']))
	# loop through physio signals
	for p in params[d]['physio_labels']:
		run_peak_average(d, p, l_window=w_len, r_window=w_len, min_peak_thres=2, 
		                 peak_distance=w_len, out_dir='results/peak_avg')


# get cross-correlation peak maps across all datasets and physiological signals
# set maximum lag as 30 secs
os.makedirs('results/crosscorr', exist_ok=True)
for d in datasets:
	print(f'run cross correlation maps for {d}')
	# get max lag (30s in TRs)
	max_lag = int(np.ceil(30/params[d]['tr']))
	# loop through physio signals
	for p in params[d]['physio_labels']:
		run_crosscorr(d, p, max_lag=max_lag, out_dir='results/crosscorr')


# # run task averaging across task-fmri datasets 
# the expand argument (int) adds TRs to expand the 
# end of the task block to account for lag of HRF
os.makedirs('results/task_avg', exist_ok=True)
print('run task averaging for chang_bh')
run_task_avg('chang_bh', expand=0, out_dir='results/task_avg')
print('run task averaging for chang_cue')
run_task_avg('chang_cue', expand=0, out_dir='results/task_avg')
print('run task averaging for nki')
run_task_avg('nki', expand=2, out_dir='results/task_avg')






