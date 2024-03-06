import json
import numpy as np
import os

from run_pca import run_pca

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






