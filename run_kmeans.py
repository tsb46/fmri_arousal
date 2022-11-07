import argparse
import numpy as np
import pickle

from scipy.stats import zscore
from sklearn.cluster import KMeans
from utils.load_write import load_data, write_nifti


def cluster_voxels(func_data, n_clusters):
    kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(func_data.T)
    return kmeans.cluster_centers_, kmeans.labels_


def write_results(dataset, cluster_centers, cluster_labels, zero_mask, n_vert, params):
    analysis_str = f'{dataset}_kmeans_group'   
    cluster_labels = cluster_labels[np.newaxis, :]+1
    write_nifti(cluster_labels, f'{analysis_str}', zero_mask, n_vert, params['mask'])
    pickle.dump([cluster_centers, cluster_labels], open(f'{analysis_str}_results.pkl', 'wb'))


def run_main(dataset, n_clusters, center, regress_global_sig):
    func_data, _, _, zero_mask, n_vert, params = load_data(dataset, 'group', physio=None, load_physio=False, 
                                                           regress_global=regress_global_sig) 
    func_data = zscore(func_data)
    # If specified, center along rows
    if center == 'r':
        func_data -= func_data.mean(axis=1, keepdims=True)
    print('k-means clustering')
    cluster_centers, cluster_labels = cluster_voxels(func_data, n_clusters)

    write_results(dataset, cluster_centers, cluster_labels, zero_mask, n_vert, params)


if __name__ == '__main__':
    """Run main analysis"""
    parser = argparse.ArgumentParser(description='Run KMeans clustering of voxels')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to run analysis on',
                        choices=['chang', 'chang_bh', 'nki', 'yale', 'hcp_fix', 'spreng'], 
                        required=True,
                        type=str)
    parser.add_argument('-n', '--n_clusters',
                        help='<Required> Number of clusters',
                        required=True,
                        type=int)
    parser.add_argument('-c', '--center',
                        help='Whether to center along the columns (c) or rows (r)',
                        default='c',
                        choices=['c','r'],
                        type=str)
    parser.add_argument('-g', '--regress_global_sig',
                        help='Whether to regress out global signal from functional data',
                        default=0,
                        required=False,
                        type=int)
    args_dict = vars(parser.parse_args())
    run_main(args_dict['dataset'], args_dict['n_clusters'],  
             args_dict['center'], args_dict['regress_global_sig'])
