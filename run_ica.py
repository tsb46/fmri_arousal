import argparse
import numpy as np
import fbpca
import pickle

from scipy.stats import zscore
from sklearn.decomposition import FastICA
from utils.load_write import load_data, write_nifti


def ica(input_data, n_comps):
    ica = FastICA(whiten=True, n_components=n_comps, max_iter=500)
    ica.fit(input_data)
    sources = ica.transform(input_data)
    return ica.components_, sources


def write_results(dataset, ica_type, spatial_map, ica_ts, zero_mask, n_vert):
    analysis_str = f'{dataset}_ica_{ica_type}_group'
    pickle.dump([spatial_map, ica_ts], open(f'{analysis_str}_results.pkl', 'wb'))
    write_nifti(spatial_map, f'{analysis_str}', zero_mask, n_vert)


def run_main(dataset, n_comps, ica_type, regress_global_sig):
    func_data, _, _, zero_mask, n_vert, _ = load_data(dataset, 'group', physio=None, load_physio=False, 
                                                      regress_global=regress_global_sig) 
    # Normalize data
    if ica_type == 'spatial':
        func_data = zscore(func_data.T)
    elif ica_type == 'temporal':
        func_data = zscore(func_data)

    # Run ICA
    unmixing_matrix, ica_comps = ica(func_data, n_comps)

    if ica_type == 'spatial':
        spatial_map = ica_comps.T
        ts = unmixing_matrix
    elif ica_type == 'temporal':
        spatial_map = unmixing_matrix
        ts = ica_comps

    write_results(dataset, ica_type, spatial_map, ts, zero_mask, n_vert)


if __name__ == '__main__':
    """Run main analysis"""
    parser = argparse.ArgumentParser(description='Run ICA analysis')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to run analysis on',
                        choices=['chang', 'nki', 'yale', 'hcp', 'hcp_fix'], 
                        required=True,
                        type=str)
    parser.add_argument('-n', '--n_comps',
                        help='<Required> Number of components from ICA',
                        required=True,
                        type=int)
    parser.add_argument('-t', '--ica_type',
                        help='Calculate temporal or spatial ICA',
                        default='temporal',
                        choices=['temporal', 'spatial'],
                        type=str)
    parser.add_argument('-g', '--regress_global_sig',
                        help='Whether to regress out global signal from functional data',
                        default=0,
                        required=False,
                        type=int)
    args_dict = vars(parser.parse_args())
    run_main(args_dict['dataset'], args_dict['n_comps'], 
             args_dict['ica_type'], args_dict['regress_global_sig'])
