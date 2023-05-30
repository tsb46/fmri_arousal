import json
import neurokit2 as nk
import numpy as np
import os
import pandas as pd
import pickle
import statsmodels.api as sm


from patsy import dmatrix
from scipy.io import loadmat
from scipy.stats import zscore
from sklearn.model_selection import TimeSeriesSplit
from utils.glm_utils import onsets_to_block
from utils.signal.butterworth_filters import butterworth_filter


# Global variables
fs_chang = 1/2.1
fs_hcp = 1/0.72
fs_nki = 1/1.4
fs_spreng = 1/3.0
fs_yale = 1.0


def commonality_analysis(pc_norm, full_pred, pred_label, pred_label_col_indx):
    sm_fit = sm.OLS(pc_norm, full_pred, hasconst=True).fit()
    full_r2 = sm_fit.rsquared
    common_r2 = full_r2.copy()
    unique_r2 = {}
    for i, label in enumerate(pred_label):
        partial_pred_i = [col_n for col_n, label_i in enumerate(pred_label_col_indx) if label_i != i]
        partial_pred = full_pred[:, partial_pred_i]
        sm_fit_partial = sm.OLS(pc_norm, partial_pred, hasconst=True).fit()
        unique_r2[label] = full_r2 - sm_fit_partial.rsquared
        common_r2 -= unique_r2[label]
    return full_r2, common_r2, unique_r2


def evaluate_model(lag_vec, model, spline_basis, var_eval):
    # Create basis from model evaluation using previously defined design matrix (for model fit)
    basis_lag_pred = dmatrix(spline_basis.design_info, {'x': lag_vec}, return_type='dataframe')
     # Intialize basis matrix
    pred_list = [var_eval * basis_lag_pred.iloc[:, l].values for l in range(basis_lag_pred.shape[1])]
    pred_mat = np.vstack(pred_list).T
    # Get predictions from model
    pred_lags = model.predict(pred_mat)
    return pred_lags


def lag_basis(var, lag_vec, lag_df):
    # Create Cubic B-spline basis for predictor and lag
    lag_splines = dmatrix("cr(x, df=lag_df) - 1", 
                          {"x": lag_vec}, return_type='dataframe')
    # Intialize basis matrix
    basis_lag = np.zeros((len(var), lag_splines.shape[1]))
    # Loop through lag bases and multiply column pairs
    lag_mat = pd.concat([var.shift(l) for l in lag_vec], axis=1)
    for l in np.arange(lag_splines.shape[1]):
        basis_lag[:, l] = np.dot(lag_mat.values, lag_splines.iloc[:,l].values)
    return basis_lag, lag_splines


def timeseries_cv(model, X, y, n_splits, gap, max_train_size, test_size):
    ts_cv = TimeSeriesSplit(
        n_splits=n_splits,
        gap=gap,
        max_train_size=max_train_size,
        test_size=test_size,
    )

    all_splits = list(ts_cv.split(X, y))
    corr_split = []
    for i in range(len(all_splits)):
        train = all_splits[i][0]
        test = all_splits[i][1]
        X_train = X[train, :]
        y_train = y[train]
        X_test = X[test, :]
        y_test = y[test]
        model_fit = model.fit(X_train, y_train)
        y_pred = model.predict(X_test)
        corr_split.append(np.corrcoef(y_test, y_pred)[0,1])
    return corr_split


def load_pca_results(subject_df, pca_fp, n_len, pca_c_fp=None, scan_ind=False, scan_len=None):
    pca_output = pickle.load(open(pca_fp, 'rb'))

    pc_ts = pca_output['pc_scores']

    if pca_c_fp is not None:
        pca_output_c = pickle.load(open(pca_c_fp, 'rb'))
        pc_ts_c = pca_output_c['pc_scores']
        pc_ts_c = pc_ts_c

    n_t = 0
    pc_subj = {}
    pc_c_subj = {}

    if scan_len is None:
        n_len = np.repeat(n_len, subject_df.shape[0])
    else:
        n_len = scan_len
    for i, (subj, subj_len) in enumerate(zip(subject_df.subject, n_len)):
        if scan_ind:
            try:
                scan = int(subject_df.iloc[i,:]['scan'])
            except ValueError:
                scan = subject_df.iloc[i,:]['scan']

            subj_str = f'{subj}_{scan}'
        else:
            subj_str = f'{subj}'
        # Extract subject PC time series
        pc_subj[subj_str] = pc_ts[n_t:(n_t+subj_len),:]
        if pca_c_fp is not None:
            pc_c_subj[subj_str] = pc_ts_c[n_t:(n_t+subj_len),:]

        n_t += subj_len
    return pc_subj, pc_c_subj


def load_raw_physio(subj, dataset, fs=None, fs_resamp=None, scan=None):
    if (dataset == 'chang') | (dataset == 'chang_bh'):
        if scan < 10:
            scan_str = f'000{scan}'
        else:
            scan_str = f'00{scan}'
        if dataset == 'chang':
            p_fp = f'data/dataset_chang/physio/raw/sub_00{subj}-mr_{scan_str}-ecr_echo1_physOUT.mat'
        elif dataset == 'chang_bh':
            p_fp = f'data/dataset_chang_bh/physio/raw/sub_00{subj}-mr_{scan_str}-adb_echo1_physOUT.mat'
        physio_raw = loadmat(p_fp, squeeze_me=True)
        fs = 1/physio_raw['OUT_p']['dt_phys'].item()
        physio = {
          'ppg': physio_raw['OUT_p']['card_dat'].item(), 
          'resp': physio_raw['OUT_p']['resp'].item()['wave'].item()
        }
        trim_n = int(fs*14.7)
        physio = {l: physio[l][trim_n:] for l in physio.keys()}

    elif dataset == 'hcp':
        if fs is None:
            raise Exception('Sampling frequency - fs - must be supplied for HCP data')
        p_fp = f'data/dataset_hcp/physio/raw/{subj}_{scan}1_physio.txt'
        physio_signals = np.loadtxt(p_fp)
        physio = { 
          'resp': physio_signals[:,1],
          'ppg': physio_signals[:,2]
        }

    elif dataset == 'spreng':
        p_fp = f'data/dataset_spreng/physio/raw/{subj}_task-rest_run-02_physio.tsv.gz'
        subj_label = os.path.basename(p_fp).split('_')[0]
        subj_dir = os.path.dirname(p_fp)
        physio_json = json.load(open(f'{subj_dir}/{subj_label}_task-rest_physio.json'))
        physio_df = pd.read_csv(p_fp, compression='gzip', sep='\t', header=None)
        fs = physio_json['SamplingFrequency']
        physio_df.columns = physio_json['Columns']
        physio_df.rename(columns = {'cardiac': 'ppg', 'respiratory': 'resp'}, inplace=True)
        physio = physio_df.to_dict('list')
        # Trim the beginning of the physio (12sec) to account for the removal of the first four functional volumes
        trim_n = fs*12
        physio = {l: physio[l][trim_n:] for l in physio.keys()}

    elif dataset == 'nki':
        fp = f'data/dataset_nki/physio/raw/{subj}_task_breathhold_physio.tsv.gz'
        physio_df = pd.read_csv(fp, compression='gzip', sep='\t', header=None)
        physio_df = physio_df.dropna(axis=1, how='all')
        subj_base = fp.rsplit('.tsv')[0]
        physio_json = json.load(open(f'{subj_base}.json'))
        physio_df.columns = physio_json['Columns']
        physio_df.rename(columns={'cardiac': 'ppg'}, inplace=True)
        physio = physio_df.to_dict('list')
        fs = physio_json['SamplingFrequency']

    if fs_resamp is not None:
        physio = {l: nk.signal_resample(physio[l], sampling_rate=fs, desired_sampling_rate=fs_resamp, method='FFT') 
                  for l in physio.keys()}
    return physio



def load_subj_chang(subj, scan, pc_ts, fs, norm=True, bp_filter=True):
    if scan < 10:
        scan_str = f'000{scan}'
    else:
        scan_str = f'00{scan}'
        
    p_str = f'data/dataset_chang/physio/proc1_physio/sub_00{subj}_mr_{scan_str}_physio.csv'
    physio_df = pd.read_csv(p_str)
    e_str = f'data/dataset_chang/eeg/proc1_fbands/sub_00{subj}_mr_{scan_str}_fbands.csv'
    eeg_bands = pd.read_csv(e_str)
    csf = np.loadtxt(f'data/dataset_chang/physio/proc1_physio/sub_00{subj}_mr_{scan_str}_csf.txt')
    gs = np.loadtxt((f'data/dataset_chang/physio/proc1_physio/sub_00{subj}_mr_{scan_str}_global_sig.txt'))
    df = pd.concat([physio_df, eeg_bands], axis=1)
    df['csf'] = csf
    df['global_sig'] = gs
    infraslow = df.pop('Infraslow')
    vigilance_ad = df['vigilance_ad'].copy()
    vigilance_at = df['vigilance_at'].copy()
    resp_rvt = df['RESP_RVT_NK'].copy()
    if bp_filter:
        df = df.apply(lambda x: butterworth_filter(x, 0.01, 0.1, fs=fs, filter_type='bandpass'), axis=0)
    df['Infraslow'] = infraslow
    df['pc1'] = pc_ts[:,0]
    df['pc2'] = pc_ts[:,1]
    df['pc3'] = pc_ts[:,2]*-1 # sign flip
    df['vigilance_ad'] = vigilance_ad
    df['vigilance_at'] = vigilance_at
    df['vigilance_ad_low'] =  butterworth_filter(vigilance_ad, None, 0.01, fs=fs, filter_type='lowpass')
    df['vigilance_at_low'] = butterworth_filter(vigilance_at, None, 0.01, fs=fs, filter_type='lowpass')
    if norm:
        df = df.apply(zscore, axis=0)
    df.reset_index(inplace=True)
    df = df.rename(columns = {'index':'time'})
    return df


def load_subj_chang_bh(subj, scan, pc_ts, fs, norm=True, bp_filter=True):
    if scan < 10:
        scan_str = f'000{scan}'
    else:
        scan_str = f'00{scan}'
        
    p_str = f'data/dataset_chang_bh/physio/proc1_physio/sub_00{subj}_mr_{scan_str}_physio.csv'
    physio_df = pd.read_csv(p_str)
    e_str = f'data/dataset_chang_bh/eeg/proc1_fbands/sub_00{subj}_mr_{scan_str}_fbands.csv'
    eeg_bands = pd.read_csv(e_str)
    csf = np.loadtxt(f'data/dataset_chang_bh/physio/proc1_physio/sub_00{subj}_mr_{scan_str}_csf.txt')
    gs = np.loadtxt((f'data/dataset_chang_bh/physio/proc1_physio/sub_00{subj}_mr_{scan_str}_global_sig.txt'))
    df = pd.concat([physio_df, eeg_bands], axis=1)
    df['csf'] = csf
    df['gs'] = gs
    infraslow = df.pop('Infraslow')
    vigilance_ad = df.pop('vigilance_ad')
    vigilance_at = df.pop('vigilance_at')
    if bp_filter:
        df = df.apply(lambda x: butterworth_filter(x, 0.01, 0.1, fs=fs, filter_type='bandpass'), axis=0)
    df['Infraslow'] = infraslow
    df['pc1'] = pc_ts[:,0]
    df['pc2'] = pc_ts[:,1]
    df['pc3'] = pc_ts[:,2]
    df['vigilance_ad'] = vigilance_ad
    df['vigilance_at'] = vigilance_at
    df['vigilance_ad_low'] =  butterworth_filter(vigilance_ad, None, 0.01, fs=fs, filter_type='lowpass')
    df['vigilance_at_low'] = butterworth_filter(vigilance_at, None, 0.01, fs=fs, filter_type='lowpass')
    if norm:
        df = df.apply(zscore, axis=0)
    df.reset_index(inplace=True)
    df = df.rename(columns = {'index':'time'})
    return df


def load_subj_chang_cue(subj, scan, pc_ts, fs, norm=True, bp_filter=True):
    if scan < 10:
        scan_str = f'000{scan}'
    else:
        scan_str = f'00{scan}'
        
    p_str = f'data/dataset_chang_cue/physio/proc1_physio/sub_00{subj}_mr_{scan_str}_physio.csv'
    physio_df = pd.read_csv(p_str)
    e_str = f'data/dataset_chang_cue/eeg/proc1_fbands/sub_00{subj}_mr_{scan_str}_fbands.csv'
    eeg_bands = pd.read_csv(e_str)
    csf = np.loadtxt(f'data/dataset_chang_cue/physio/proc1_physio/sub_00{subj}_mr_{scan_str}_csf.txt')
    df = pd.concat([physio_df, eeg_bands], axis=1)
    df['csf'] = csf
    infraslow = df.pop('Infraslow')
    vigilance_ad = df.pop('vigilance_ad')
    vigilance_at = df.pop('vigilance_at')
    if bp_filter:
        df = df.apply(lambda x: butterworth_filter(x, 0.01, 0.1, fs=fs, filter_type='bandpass'), axis=0)
    df['Infraslow'] = infraslow
    df['pc1'] = pc_ts[:,0]
    df['pc2'] = pc_ts[:,1]
    df['pc3'] = pc_ts[:,2]
    df['vigilance_ad'] = vigilance_ad
    df['vigilance_at'] = vigilance_at
    df['vigilance_ad_low'] =  butterworth_filter(vigilance_ad, None, 0.01, fs=fs, filter_type='lowpass')
    df['vigilance_at_low'] = butterworth_filter(vigilance_at, None, 0.01, fs=fs, filter_type='lowpass')
    if norm:
        df = df.apply(zscore, axis=0)
    df.reset_index(inplace=True)
    df = df.rename(columns = {'index':'time'})
    return df


def load_subj_hcp(subj, pc_ts, fs, norm=True, bp_filter=True):
    p_str = f'data/dataset_hcp/physio/proc1_physio/{subj}_physio.csv'
    csf_str = f'data/dataset_hcp/physio/proc1_physio/{subj}_csf.txt'
    global_str = f'data/dataset_hcp/physio/proc1_physio/{subj}_global_sig.txt'
    df = pd.read_csv(p_str)
    csf = np.loadtxt(csf_str)
    global_sig = np.loadtxt(global_str)
    df['csf'] = csf
    df['global_sig'] = global_sig
    ## IMPORTANT! CSF signal has large amplitude spikes at start of scan - set first 10 time points to median
    df['csf'].iloc[:10] = df['csf'].median()
    df['global_sig'].iloc[:10] = df['global_sig'].median()
    if bp_filter:
        df = df.apply(lambda x: butterworth_filter(x, 0.01, 0.1, fs=fs, filter_type='bandpass'), axis=0)
    df['pc1'] = pc_ts[:,0]
    df['pc2'] = pc_ts[:,1]*-1 # sign flip to keep consistent with Chang PCA
    df['pc3'] = pc_ts[:,2]*-1 # sign flip to keep consistent with Chang PCA
    if norm:
        df = df.apply(zscore, axis=0)
    df.reset_index(inplace=True)
    df = df.rename(columns = {'index':'time'})
    return df


def load_subj_nki(subj, pc_ts, fs, norm=True):
    p_str = f'data/dataset_nki/physio/proc1_physio/{subj}_task_breathhold_physio.csv'
    csf_str = f'data/dataset_nki/physio/proc1_physio/{subj}_task_breathhold_physio_csf.txt'
    df = pd.read_csv(p_str)
    csf = np.loadtxt(csf_str)
    df['csf'] = csf
    df = df.apply(lambda x: butterworth_filter(x, 0.01, 0.1, fs=fs, filter_type='bandpass'), axis=0)
    df['pc1'] = pc_ts[:,0]*-1 # sign flip to keep consistent with Chang PCA
    df['pc2'] = pc_ts[:,1]
    if norm:
        df = df.apply(zscore, axis=0)
    df.reset_index(inplace=True)
    df = df.rename(columns = {'index':'time'})
    event_df = pd.read_csv(f'data/dataset_nki/events/{subj}_task_breathhold_events.tsv', sep='\t')
    event_df.loc[event_df.trial_type.str.startswith('H'), 'trial_type'] = 'H'
    event_df.loc[event_df.trial_type.isin(['G', 'Deep','In','Out']), 'trial_type'] = 'B'
    for event in event_df.trial_type.unique():
        trial_df = event_df.loc[event_df.trial_type == event]
        trial_ind = onsets_to_block(trial_df, 186, 1.4)
        block_ts = np.zeros(186)
        for block_i in trial_ind:
            block_ts[block_i] = 1
        df[f'trial_{event}'] = block_ts
    return df


def load_subj_spreng(subj, scan, pc_ts, fs, norm=True):
    p_str = f'data/dataset_spreng/physio/proc1_physio/{subj}_task-rest_{scan}_physio.csv'
    csf_str = f'data/dataset_spreng/physio/proc1_physio/{subj}_task-rest_{scan}_physio_csf.txt'
    df = pd.read_csv(p_str)
    csf = np.loadtxt(csf_str)
    df['csf'] = csf
    df = df.apply(lambda x: butterworth_filter(x, 0.01, 0.1, fs=fs, filter_type='bandpass'), axis=0)
    df['pc1'] = pc_ts[:,0]
    df['pc2'] = pc_ts[:,1]
    df['pc3'] = pc_ts[:,2]*-1 # sign flip to keep consistent with Chang PCA
    if norm:
        df = df.apply(zscore, axis=0)
    df.reset_index(inplace=True)
    df = df.rename(columns = {'index':'time'})
    return df


def load_subj_yale(subj, scan, pc_ts, fs, norm=True):
    p_str = f'data/dataset_yale/physio/raw/{subj}_task-rest_run-0{scan}_et.txt'
    csf_str = f'data/dataset_yale/physio/proc1_physio/{subj}_task-rest_run-0{scan}_bold_csf.txt'
    df = pd.read_csv(p_str, sep=' ', header=None)
    df.columns = ['pupil']
    df['csf'] = np.loadtxt(csf_str)
    df = df.apply(lambda x: butterworth_filter(x, 0.01, 0.1, fs=fs, filter_type='bandpass'), axis=0)
    df['pc1'] = pc_ts[:,0]*-1 # sign flip to keep consistent with Chang PCA
    df['pc2'] = pc_ts[:,1]*-1
    df['pc3'] = pc_ts[:,2]
    if norm:
        df = df.apply(zscore, axis=0)
    df.reset_index(inplace=True)
    df = df.rename(columns = {'index':'time'})
    return df


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