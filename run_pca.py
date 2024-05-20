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
    # hilbert transform
    input_data = hilbert(input_data, axis=0)
    return input_data.conj()


def pca(input_data, n_comps, n_iter=10):
    # compute pca
    # get number of observations
    n_samples = input_data.shape[0]
    # fbpca pca
    (U, s, Va) = fbpca.pca(input_data, k=n_comps, n_iter=n_iter)
    # calc explained variance
    explained_variance_ = ((s ** 2) / (n_samples - 1)) / input_data.shape[1]
    total_var = explained_variance_.sum()
    # compute PC scores
    pc_scores = input_data @ Va.T
    # get loadings from eigenvectors
    loadings =  Va.T @ np.diag(s) 
    loadings /= np.sqrt(input_data.shape[0]-1)
    # package outputs
    output_dict = {'U': U,
                   's': s,
                   'Va': Va,
                   'loadings': loadings.T,
                   'exp_var': explained_variance_,
                   'pc_scores': pc_scores
                   }   
    return output_dict


def rotation(pca_output, data, rotation):
    # rotate PCA weights, if specified, and recompute pc scores
    if rotation == 'varimax':
        rotated_weights, r_mat = varimax(pca_output['loadings'].T)
        pca_output['r_mat'] = r_mat
    elif rotation == 'promax':
        rotated_weights, r_mat, phi_mat = promax(pca_output['loadings'].T)
        pca_output['r_mat'] = r_mat
        pca_output['phi_mat'] = phi_mat
    # https://stats.stackexchange.com/questions/59213/how-to-compute-varimax-rotated-principal-components-in-r
    # recompute pc scores
    projected_scores = data @ pinv(rotated_weights).T
    pca_output['loadings'] = rotated_weights.T
    pca_output['pc_scores'] = projected_scores
    return pca_output


def write_results(dataset, pca_output, pca_type, params, 
                  rotate, out_dir):
    # write out results of pca analysis
    # create output name
    if out_dir is not None:
        analysis_str = f'{out_dir}/{dataset}_pca'
    else:
        analysis_str = f'{dataset}_pca'

    if rotate:
        analysis_str += f'_{rotate}'

    # get loadings
    loadings = pca_output['loadings']
    # Write nifti
    if pca_type == 'complex': 
        analysis_str = analysis_str + '_c'
        loadings_real = np.real(loadings)
        loadings_imag = np.imag(loadings)
        loadings_ang = np.angle(loadings)
        loadings_abs = np.abs(loadings)
        write_nifti(loadings_abs, f'{analysis_str}_magnitude', params)
        write_nifti(loadings_real, f'{analysis_str}_real', params)
        write_nifti(loadings_imag, f'{analysis_str}_imag', params)
        write_nifti(loadings_ang, f'{analysis_str}_ang', params)        
    elif pca_type == 'real':
        write_nifti(loadings, f'{analysis_str}', params)
    pickle.dump(pca_output, open(f'{analysis_str}_results.pkl', 'wb'))


def run_pca(dataset, n_comps, pca_type, rotate, out_dir=None):
    # load dataset
    func_data, _, params = load_data(dataset, physio=None) 
    # if pca_type is complex, compute hilbert transform
    if pca_type == 'complex':
        func_data = hilbert_transform(func_data)

    # compute pca
    pca_output = pca(func_data, n_comps)

    # rotate pca weights, if specified
    if rotate is not None:
        pca_output = rotation(pca_output, func_data, rotate)

    # write out results
    write_results(dataset, pca_output, pca_type, params, 
                  rotate, out_dir)


if __name__ == '__main__':
    """PCA/CPCA Modeling"""
    parser = argparse.ArgumentParser(description='Run PCA or CPCA analysis')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to run analysis on',
                        choices=['chang', 'chang_bh', 'chang_cue', 
                                 'natview', 'nki', 'nki_rest', 'hcp', 
                                 'spreng', 'toronto', 
                                 'yale'], 
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
    parser.add_argument('-r', '--rotate',
                        help='Whether to rotate pca weights',
                        default=None,
                        required=False,
                        choices=['varimax', 'promax'],
                        type=str)
    args_dict = vars(parser.parse_args())
    run_pca(args_dict['dataset'], args_dict['n_comps'], 
            args_dict['pca_type'], args_dict['rotate'])
