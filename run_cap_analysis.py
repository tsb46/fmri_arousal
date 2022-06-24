import argparse
import numpy as np
import pickle

from itertools import zip_longest
from scipy.stats import zscore
from sklearn.cluster import KMeans
from utils.load_utils import load_subject_list
from utils.load_write import load_data, write_nifti


def cluster_maps(selected_maps, norm_maps, n_clusters):
    if norm_maps:
        selected_maps = zscore(selected_maps.T).T
    kmeans = KMeans(n_clusters=n_clusters, random_state=0).fit(selected_maps)
    return kmeans.cluster_centers_, kmeans.labels_


def get_suprathreshold_maps(seed_signal, func_data, z_thres, lag):
    seed_ts_len = len(seed_signal)
    seed_signal_z = zscore(seed_signal)
    selected_vals = np.where(seed_signal_z >= z_thres)[0]
    selected_vals += lag
    bound_mask = np.isin(selected_vals, np.arange(seed_ts_len))
    selected_vals = selected_vals[bound_mask]
    return selected_vals, np.squeeze(func_data[selected_vals, :]) 


def run_main(dataset, n_clusters, seed, norm_maps, lag, z_thres, sign, save_suprathres_maps):
    subject_df = load_subject_list(dataset)
    if dataset == 'chang':
        subj_list = subject_df.subject
        scan_list = subject_df.scan
    else:
        subj_list = subject_df.subject
        scan_list = [None]

    supra_thres_maps = []
    supra_thres_tps = []
    # Load through subjects and identify suprathreshold time points
    for subj, scan in zip_longest(subj_list, scan_list):
        print(subj)
        func_data, physio_sig, physio_labels, zero_mask, n_vert, params = load_data(dataset, 'subject', physio=[seed],
                                                                                    load_physio=True, subj_n=subj, 
                                                                                    scan_n=scan, verbose=False, 
                                                                                    group_method='list',
                                                                                    filter_nan_voxels=False) 
        func_data = zscore(func_data[0])
        #Replace NaN voxels with zero
        inds = np.where(np.isnan(func_data))
        func_data[inds] = 0
        # Select suprathreshold maps
        if sign == 'p':
            seed_ts = physio_sig[0]
        elif sign == 'n':
            seed_ts = physio_sig[0]*-1
        selected_tps, selected_maps = get_suprathreshold_maps(seed_ts, func_data, z_thres, lag)
        supra_thres_maps.append(selected_maps)
        supra_thres_tps.append(selected_tps)
    # Concatenate subject suprathreshold maps into matrix
    supra_thres_maps = np.vstack(supra_thres_maps)
    # Cluster maps
    cluster_centroid, cluster_indx = cluster_maps(supra_thres_maps, norm_maps, n_clusters)
    # Write results
    results_list = [cluster_centroid, cluster_indx, supra_thres_tps, supra_thres_maps]
    results_labels = ['c_centroids', 'c_indx', 'supra_thres_tps', 'supra_thres_maps']
    write_results(results_list, results_labels, seed, norm_maps, sign, dataset, 
                  zero_mask, n_vert, save_suprathres_maps)

def write_results(results_list, results_labels, seed, norm, sign, dataset,
                   zero_mask, n_vert, save_suprathres_maps):
    results_dict = {l: r for r, l in zip(results_list, results_labels)}
    analysis_str = f'{dataset}_cap_{seed}'
    if norm:
        analysis_str += f'_norm'
    if sign == 'n':
        analysis_str += f'_neg'

    if ~save_suprathres_maps:
        suprathres_maps = results_dict.pop('supra_thres_maps')
    else:
        write_nifti(results_dict['supra_thres_maps'], f'{analysis_str}_suprathres_maps', 
                        zero_mask, n_vert)

    pickle.dump(results_dict, open(f'{analysis_str}_results.pkl', 'wb'))
    write_nifti(results_dict['c_centroids'], analysis_str, zero_mask, n_vert)
            




if __name__ == '__main__':
    """Run main analysis"""
    parser = argparse.ArgumentParser(description='Run CAP analysis on seed time course')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to run analysis on',
                        choices=['chang', 'nki', 'yale', 'hcp', 'hcp_fix'], 
                        required=True,
                        type=str)
    parser.add_argument('-e', '--seed', 
                        help='<Required> choice of seed ("precuneus" or "superior_parietal")',
                        choices=['precuneus', 'superior_parietal', 'global_sig'], 
                        required=True,
                        type=str)
    parser.add_argument('-n', '--n_clusters',
                        default=2, 
                        help='Number of clusters to estimate',
                        required=False,
                        type=int)
    parser.add_argument('-m', '--norm_maps',
                        default=0, 
                        help='Whether to z-score normalize maps before clustering',
                        required=False,
                        type=int)
    parser.add_argument('-l', '--lag',
                        default=0, 
                        help='Amount of lag (shift) applied to suprathreshold time points',
                        required=False,
                        type=int)
    parser.add_argument('-t', '--z_threshold',
                        help='Z-score threshold to select time points from seed time course',
                        required=False,
                        default=2,
                        type=float)
    parser.add_argument('-sign', '--sign',
                        help='Z-score threshold to select time points from seed time course',
                        required=False,
                        choices=['p', 'n'],
                        default='p',
                        type=str)
    parser.add_argument('-save', '--save_suprathres_maps',
                        help='Save maps corresponding to suprathreshold time points of seed ts',
                        required=False,
                        default=0,
                        type=int)
    args_dict = vars(parser.parse_args())
    run_main(args_dict['dataset'], args_dict['n_clusters'], args_dict['seed'],
             args_dict['norm_maps'], args_dict['lag'], args_dict['z_threshold'], 
             args_dict['sign'], args_dict['save_suprathres_maps'])


