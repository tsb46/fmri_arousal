import argparse
import numpy as np
import pickle

from scipy.stats import zscore
from utils.load_write import load_data, write_nifti        


def create_bins(phase_ts, n_bins): 
    freq, bins = np.histogram(phase_ts, n_bins)
    bin_indx = np.digitize(phase_ts, bins)
    bin_centers = np.mean(np.vstack([bins[0:-1],bins[1:]]), axis=0)
    return bin_indx, bin_centers


def create_dynamic_phase_maps(recon_ts, bin_indx, n_bins):
    bin_timepoints = []
    for n in range(1, n_bins+1):
        ts_indx = np.where(bin_indx==n)[0]
        bin_timepoints.append(np.mean(recon_ts[ts_indx,:], axis=0))
    dynamic_phase_map = np.array(bin_timepoints)
    return dynamic_phase_map


def reconstruct_ts(pca_res, n, real=True, rotation=False):
    if rotation:
        recon_ts = pca_res['pc_scores'][:, n] @ pca_res['loadings'][n, :].conj()
    else:
        U = pca_res['U'][:,n]
        s = np.atleast_2d(pca_res['s'][n])
        Va = pca_res['Va'][n,:].conj()
        recon_ts = U @ s @ Va
    if real:
        recon_ts = np.real(recon_ts)
    else:
        recon_ts = np.imag(recon_ts)
    return recon_ts


def write_results(dataset,recon_comp, n_comp, zero_mask, n_vert, params):
    analysis_str = f'{dataset}_cpca_recon_n{n_comp}'
    write_nifti(recon_comp, analysis_str, zero_mask, n_vert, params['mask'])


def run_main(dataset, cpca_res, n_recon, n_bins, rotation, real=True):
    # Not ideal to load full group data; only needed to get brain_mask to write nifti
    _, _, _, zero_mask, n_vert, params = load_data(dataset, 'group', physio=None, load_physio=False) 
    cpca_res = pickle.load(open(cpca_res, 'rb'))
    bin_indx_all = []
    bin_centers_all = []
    for n in range(n_recon):
        recon_ts = reconstruct_ts(cpca_res, [n], real, rotation)
        phase_ts = np.angle(cpca_res['pc_scores'][:,n]) 
        bin_indx, bin_centers = create_bins(phase_ts, n_bins)
        dynamic_phase_map = create_dynamic_phase_maps(recon_ts, bin_indx, n_bins)
        bin_indx_all.append(bin_indx); bin_centers_all.append(bin_centers)
        write_results(dataset, dynamic_phase_map, n, zero_mask, n_vert, params)
    pickle.dump([bin_indx_all, bin_centers_all], open(f'{dataset}_cpca_recon_results.pkl', 'wb'))


if __name__ == '__main__':
    """Reconstruct cPCA analysis"""
    parser = argparse.ArgumentParser(description='Reconstruct cPCA Components from group cPCA analysis')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to run analysis on',
                        choices=['chang', 'chang_bh', 'nki', 'yale', 'hcp', 'hcp_fix'], 
                        required=True,
                        type=str)
    parser.add_argument('-i', '--input_cpca',
                        help='<Required> File path to cpca results pickle (must match dataset)',
                        required=True,
                        type=str)
    parser.add_argument('-n', '--n_reconstruct',
                        help='<Required> Number of components to reconstruct from cPCA',
                        required=True,
                        type=int)
    parser.add_argument('-b', '--n_bins',
                        help='Number of phase bins',
                        default=30,
                        required=False,
                        type=int)
    parser.add_argument('-r', '--rotation',
                        help='Whether this PCA solution was rotated (reconstruction is different)',
                        default=0,
                        required=False,
                        type=int)
    args_dict = vars(parser.parse_args())
    run_main(args_dict['dataset'], args_dict['input_cpca'], 
             args_dict['n_reconstruct'], args_dict['n_bins'],
             args_dict['rotation'])
