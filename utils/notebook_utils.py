import json
import nibabel as nb
import numpy as np
import os
import pandas as pd
import pickle
import statsmodels.api as sm
import warnings

from patsy import dmatrix
from scipy.io import loadmat
from scipy.stats import zscore
from utils.load_write import load_subject_list
from utils.glm_utils import onsets_to_block
from utils.signal_utils import butterworth_filter


def lag_basis(var, lag_vec, lag_df):
    # Create Cubic B-spline basis for lag
    lag_splines = dmatrix("cr(x, df=lag_df) - 1", 
                          {"x": lag_vec}, return_type='dataframe')
    # Intialize basis matrix
    basis_lag = np.zeros((len(var), lag_splines.shape[1]))
    # Loop through lag bases and multiply column pairs
    lag_mat = pd.concat([var.shift(l) for l in lag_vec], axis=1)
    for l in np.arange(lag_splines.shape[1]):
        basis_lag[:, l] = np.dot(lag_mat.values, lag_splines.iloc[:,l].values)
    return basis_lag, lag_splines
    

def load_subj(subj, scan, dataset, params, pca_ts, pca_index=None):
    # load subject 
    physio_template = params['physio']
    physio_labels = params['physio_labels']
    physio_dict = {}
    for p in physio_labels:
        if scan is None:
            phys_fp = physio_template.format(subj, p)
        else:
            phys_fp = physio_template.format(subj, scan, p)

        try: 
            physio_dict[p] = np.loadtxt(phys_fp)
        except OSError:
            print(f"warning: {dataset}: {subj}_{scan} does not have {p} signal")

    physio_df = pd.DataFrame(physio_dict)
    if pca_index is None:
        pca_cols = [f'PC{i}' for i in range(pca_ts.shape[1])]
        pca_df = pd.DataFrame(pca_ts, columns=pca_cols)
        physio_df = pd.concat([physio_df, pca_df], axis=1)
    else:
        physio_df[f'PC{pca_index+1}'] = pca_ts[:,pca_index]
    return physio_df


def load_pca_results(dataset, pca_fp, params):
    # load and reorganize PCA time courses into subject time courses
    pca_output = pickle.load(open(pca_fp, 'rb'))

    pc_ts = pca_output['pc_scores']
    
    n_t = 0
    pc_subj = {}
    pc_c_subj = {}

    # load subject list
    subject, scan = load_subject_list(dataset, params['subject_list'])
    # loop through subjects and assign PC time courses, keeping up with time index
    for i, (subj, sc) in enumerate(zip(subject, scan)):
        # get lenght of subject functional scan to keep up with time index
        func_len = nb.load(params['func'].format(subj, sc)).shape[-1]
        # specify subject string
        subj_str = f'{subj}_{sc}'
        # Extract subject PC time series
        pc_subj[subj_str] = pc_ts[n_t:(n_t+func_len),:]
        # push time index forward by length of subject functional scan
        n_t += func_len
    return pc_subj


def xcorr_select_max(c, lags, constrain='abs'):
    if constrain == 'abs':
        max_r = c[np.argsort(np.abs(c))[-1]]
        max_lag = lags[np.argsort(np.abs(c))[-1]]
    elif constrain == 'pos':
        max_r = c[np.argsort(c)[-1]]
        max_lag = lags[np.argsort(c)[-1]]
    elif constrain == 'neg':
        max_r = c[np.argsort(-c)[-1]]
        max_lag = lags[np.argsort(-c)[-1]]
    
    return max_r, max_lag