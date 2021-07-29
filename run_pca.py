import argparse
import numpy as np 
import fbpca
import pickle

from utils.load_write import load_data, write_nifti
from scipy.signal import hilbert
from scipy.stats import zscore


def hilbert_transform(input_data):
    complex_data = hilbert(input_data, axis=0)
    return complex_data.conj()


def pca(input_data, n_comps, n_iter=10):
	n_samples = input_data.shape[0]
	(U, s, Va) = fbpca.pca(input_data, k=n_comps, n_iter=n_iter)
	explained_variance_ = (s ** 2) / (n_samples - 1)
	total_var = explained_variance_.sum()
	pc_scores = input_data @ Va.T
	loadings =  Va.T @ np.diag(s) 
	loadings /= np.sqrt(input_data.shape[0]-1)
	output_dict = {'U': U,
				   's': s,
				   'Va': Va,
				   'loadings': loadings.T,
				   'exp_var': explained_variance_,
				   'pc_scores': pc_scores
				   }   
	return output_dict


def write_results(level, mask, pca_output, pca_type, comp_weights, 
                  subj_n, zero_mask, n_vert):
	if level == 'group':
		analysis_str = 'pca_group'
	elif level == 'subject':
		analysis_str = f'pca_s{subj_n}'

	# Write nifti
	if pca_type == 'complex': 
		analysis_str = 'c'+ analysis_str
		comp_weights_real = np.real(comp_weights)
		comp_weights_imag = np.imag(comp_weights)
		comp_weights_ang = np.angle(comp_weights)
		write_nifti(comp_weights_real, mask, f'{analysis_str}_real', zero_mask, n_vert)
		write_nifti(comp_weights_imag, mask, f'{analysis_str}_imag', zero_mask, n_vert)
		write_nifti(comp_weights_ang, mask, f'{analysis_str}_ang', zero_mask, n_vert)        
	elif pca_type == 'real':
		write_nifti(comp_weights, mask, f'{analysis_str}', zero_mask, n_vert)
	pickle.dump(pca_output, open(f'{analysis_str}_results.pkl', 'wb'))


def run_main(n_comps, level, subj_n, pca_type, center, physio_params, mask):
	if level == 'subject' and subj_n is None:
		raise Exception('Subject number must be provided for subject-level analysis')

	func_data, eeg, rv, hr, zero_mask, n_vert = load_data(level, mask, physio_params, subj_n)
	# Normalize data
	func_data = zscore(func_data)
	# If specified, center along rows
	if center == 'r':
		func_data -= func_data.mean(axis=1, keepdims=True)

	if pca_type == 'complex':
		func_data = hilbert_transform(func_data)
	pca_output = pca(func_data, n_comps)
	write_results(level, mask, pca_output, pca_type, 
				  pca_output['loadings'], subj_n,
				  zero_mask, n_vert)


if __name__ == '__main__':
	"""Run main analysis"""
	parser = argparse.ArgumentParser(description='Run PCA or CPCA analysis')
	parser.add_argument('-n', '--n_comps',
						help='<Required> Number of components from PCA',
						required=True,
						type=int)
	parser.add_argument('-l', '--level',
						help='subject or group level analysis',
						default='group',
						choices=['subject', 'group'],
						type=str)
	parser.add_argument('-s', '--subject_n',
						help='subject number for subject level analysis',
						default=None,
						type=int)
	parser.add_argument('-t', '--pca_type',
						help='Calculate complex or real PCA',
						default='real',
						choices=['real', 'complex'],
						type=str)
	parser.add_argument('-c', '--center',
						help='Whether to center along the columns (c) or rows (r)',
						default='c',
						choices=['c','r'],
						type=str)
	parser.add_argument('-p', '--physio_params', 
						help='file path to preprocessing params for physio signals',
						default='physio_proc_params.json',
						type=str)
	parser.add_argument('-m', '--mask', 
						help='file path to brain mask',
						default='masks/MNI152_T1_3mm_brain_mask.nii.gz',
						type=str)
	args_dict = vars(parser.parse_args())
	run_main(args_dict['n_comps'], args_dict['level'], args_dict['subject_n'], 
			 args_dict['pca_type'], args_dict['center'], 
			 args_dict['physio_params'], args_dict['mask'])
