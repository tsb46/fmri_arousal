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


def create_dynamic_phase_maps(data, bin_indx, n_bins, weights):
    bin_timepoints = []
    for n in range(1, n_bins+1):
        ts_indx = np.where(bin_indx==n)[0]
        bin_timepoints.append(
            np.average(data[ts_indx,:], weights=weights[n-1], axis=0)
        )
    dynamic_phase_map = np.array(bin_timepoints)
    return dynamic_phase_map


def reconstruct_ts(pca_res, n, rotation, real=True):
    # reconstruct ts 
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


def write_results(dataset,recon_comp, n_comp, params, out_dir):
    # write results of cpca reconstruction
    if out_dir is not None:
        analysis_str = f'{out_dir}/{dataset}_cpca_recon_n{n_comp}'
    else:
        analysis_str = f'{dataset}_cpca_recon_n{n_comp}'
    write_nifti(recon_comp, analysis_str, params)


def cpca_recon(dataset, cpca_res, n_recon, rotation, recon_method,
               params, data=None, out_dir=None, n_bins=30):
    # reconstruct cpca component 'movies' from cpca results
    bin_indx_all = []
    bin_centers_all = []
    for n in range(n_recon):
        if recon_method == 'proj':
            data = reconstruct_ts(cpca_res, [n], rotation)
            phase_ts = np.angle(cpca_res['pc_scores'][:,n])
        elif recon_method == 'weighted':
            # the complex conjugate of the PC time courses was taken
            # in the pca script for visualization purposes; to align
            # with raw time courses, take the complex conjugate again
            phase_ts = np.angle(cpca_res['pc_scores'][:,n].conj())
        # shift phase delay angles from -pi to pi -> 0 to 2*pi
        phase_ts = np.mod(phase_ts, 2*np.pi)
        # bin time courses into phase bins
        bin_indx, bin_centers = create_bins(phase_ts, n_bins)
        # create weight vector based on CPC amplitude, if 'weighted' method
        if recon_method == 'weighted':
            power_ts = np.abs(cpca_res['pc_scores'][:,n])**2
        weight_vec = []
        for i in range(1, n_bins+1):
            ts_indx = np.where(bin_indx==i)[0]
            if recon_method == 'weighted':
                w = power_ts[ts_indx]
            elif recon_method == 'proj':
                w = np.ones(len(ts_indx))
            weight_vec.append(w)

        # average time courses within bins
        dynamic_phase_map = create_dynamic_phase_maps(data, bin_indx, 
                                                      n_bins, weight_vec)
        bin_indx_all.append(bin_indx); bin_centers_all.append(bin_centers)
        write_results(dataset, dynamic_phase_map, n, params, out_dir)
    if out_dir is not None:
        pkl_out = f'{out_dir}/{dataset}_cpca_recon_results.pkl'
    else:
        pkl_out = f'{dataset}_cpca_recon_results.pkl'
    pickle.dump([bin_indx_all, bin_centers_all], open(pkl_out, 'wb'))
