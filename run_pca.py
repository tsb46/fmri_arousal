import argparse
import numpy as np
import fbpca
import pickle

from numpy.linalg import pinv
from scipy.signal import hilbert
from scipy.stats import zscore
from utils.load_write import load_data, write_nifti
from utils.rotation import varimax, promax


def hilbert_transform(input_data):
    input_data = hilbert(input_data, axis=0)
    return input_data.conj()


def pca(input_data, n_comps, n_iter=10):
    n_samples = input_data.shape[0]
    (U, s, Va) = fbpca.pca(input_data, k=n_comps, n_iter=n_iter)
    explained_variance_ = ((s ** 2) / (n_samples - 1)) / input_data.shape[1]
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


def rotation(pca_output, data, rotation):
    if rotation == 'varimax':
        rotated_weights, r_mat = varimax(pca_output['loadings'].T)
        pca_output['r_mat'] = r_mat
    elif rotation == 'promax':
        rotated_weights, r_mat, phi_mat = promax(pca_output['loadings'].T)
        pca_output['r_mat'] = r_mat
        pca_output['phi_mat'] = phi_mat
    # https://stats.stackexchange.com/questions/59213/how-to-compute-varimax-rotated-principal-components-in-r
    projected_scores = data @ pinv(rotated_weights).T
    pca_output['loadings'] = rotated_weights.T
    pca_output['pc_scores'] = projected_scores
    return pca_output


def write_results(dataset, pca_output, pca_type, comp_weights, 
                  zero_mask, n_vert, rotate, params):
    analysis_str = f'{dataset}_pca_group'

    if rotate:
        analysis_str += f'_{rotate}'

    # Write nifti
    if pca_type == 'complex': 
        analysis_str = analysis_str + '_c'
        comp_weights_real = np.real(comp_weights)
        comp_weights_imag = np.imag(comp_weights)
        comp_weights_ang = np.angle(comp_weights)
        comp_weights_abs = np.abs(comp_weights)
        write_nifti(comp_weights_abs, f'{analysis_str}_magnitude', zero_mask, n_vert, params['mask'])
        write_nifti(comp_weights_real, f'{analysis_str}_real', zero_mask, n_vert, params['mask'])
        write_nifti(comp_weights_imag, f'{analysis_str}_imag', zero_mask, n_vert, params['mask'])
        write_nifti(comp_weights_ang, f'{analysis_str}_ang', zero_mask, n_vert, params['mask'])        
    elif pca_type == 'real':
        write_nifti(comp_weights, f'{analysis_str}', zero_mask, n_vert, params['mask'])
    pickle.dump(pca_output, open(f'{analysis_str}_results.pkl', 'wb'))


def run_main(dataset, n_comps, pca_type, center, rotate, regress_global_sig):
    func_data, _, _, zero_mask, n_vert, params = load_data(dataset, physio=None, load_physio=False, 
                                                           regress_global=regress_global_sig) 
    # If specified, center along rows
    if center == 'r':
        func_data -= func_data.mean(axis=1, keepdims=True)
    if pca_type == 'complex':
        print('hilbert')
        # Convert to single float - hack to get HCP data through analysis without memory error
        func_data = np.single(func_data)
        func_data = hilbert_transform(func_data)
    print('pca')
    pca_output = pca(func_data, n_comps)

    if rotate is not None:
        pca_output = rotation(pca_output, func_data, rotate)

    write_results(dataset, pca_output, pca_type, 
                  pca_output['loadings'], zero_mask, n_vert, 
                  rotate, params)


if __name__ == '__main__':
    """Run main analysis"""
    parser = argparse.ArgumentParser(description='Run PCA or CPCA analysis')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to run analysis on',
                        choices=['chang', 'chang_bh', 'chang_cue', 'nki', 'hcp', 
                                 'spreng', 'yale', 'natview'], 
                        required=True,
                        type=str)
    parser.add_argument('-n', '--n_comps',
                        help='<Required> Number of components from PCA',
                        required=True,
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
    parser.add_argument('-r', '--rotate',
                        help='Whether to rotate pca weights',
                        default=None,
                        required=False,
                        choices=['varimax', 'promax'],
                        type=str)
    parser.add_argument('-g', '--regress_global_sig',
                        help='Whether to regress out global signal from functional data',
                        default=0,
                        required=False,
                        type=int)
    args_dict = vars(parser.parse_args())
    run_main(args_dict['dataset'], args_dict['n_comps'], 
             args_dict['pca_type'], args_dict['center'], 
             args_dict['rotate'], args_dict['regress_global_sig'])
