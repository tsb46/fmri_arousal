import statsmodels.api as sm
import numpy as np
import pandas as pd
import pickle

from patsy import dmatrix
from scipy.stats import zscore
from sklearn.model_selection import TimeSeriesSplit
from utils.signal.butterworth_filters import butterworth_filter

# Global variables
fs_chang = 1/2.1
fs_hcp = 1/0.72
fs_nki = 1/1.4
fs_spreng = 1/3.0


def commonality_analysis(pc_norm, full_pred, pred_label, pred_label_col_indx):
    sm_fit = sm.OLS(pc_norm, full_pred, hasconst=False).fit()
    full_r2 = sm_fit.rsquared
    common_r2 = full_r2.copy()
    unique_r2 = {}
    for i, label in enumerate(pred_label):
        partial_pred_i = [col_n for col_n, label_i in enumerate(pred_label_col_indx) if label_i != i]
        partial_pred = full_pred[:, partial_pred_i]
        sm_fit_partial = sm.OLS(pc_norm, partial_pred).fit()
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


def lag_basis(var, lag_vec, lag_nknots):
    # Create Cubic B-spline basis for predictor and lag
    lag_splines = dmatrix("cr(x, df=lag_nknots) - 1", 
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


def load_pca_results(subject_df, pca_fp, pca_p_fp, n_len, scan_ind=False):
    pca_output = pickle.load(open(pca_fp, 'rb'))
    pca_output_p = pickle.load(open(pca_p_fp, 'rb'))

    pc_ts = pca_output['pc_scores']
    pc_ts_p = pca_output_p['pc_scores']

    n_t = 0
    pc_subj = {}
    pc_p_subj = {}
    for i, subj in enumerate(subject_df.subject):
        if scan_ind:
            scan = int(subject_df.iloc[i,:]['scan'])
            subj_str = f'{subj}_{scan}'
        else:
            subj_str = f'{subj}'
        # Extract subject PC time series
        pc_subj[subj_str] = pc_ts[n_t:(n_t+n_len),:]
        # Extract subject PC-Promax time series
        pc_p_subj[subj_str] = pc_ts_p[n_t:(n_t+n_len),:]
        n_t += n_len
    return pc_subj, pc_p_subj
        



def load_subj_chang(subj, scan, pc_ts, pc_ts_p, fs, norm=True, bp_filter=True):
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
    df['gs'] = gs
    infraslow = df.pop('Infraslow')
    vigilance_ad = df.pop('vigilance_ad')
    vigilance_at = df.pop('vigilance_at')
    if bp_filter:
        df = df.apply(lambda x: butterworth_filter(x, 0.01, 0.1, fs=fs, filter_type='bandpass'), axis=0)
    df['Infraslow'] = infraslow
    df['pc1'] = pc_ts[:,0]
    df['pc2'] = pc_ts[:,1]
    df['pc3'] = pc_ts[:,2]*-1 # sign flip
    df['pc1_p'] = pc_ts_p[:,0]
    df['pc2_p'] = pc_ts_p[:,1]
    df['pc3_p'] = pc_ts_p[:,2]*-1 # sign flip
    df['vigilance_ad'] = vigilance_ad
    df['vigilance_at'] = vigilance_at
    df['vigilance_ad_low'] =  butterworth_filter(vigilance_ad, None, 0.01, fs=fs, filter_type='lowpass')
    df['vigilance_at_low'] = butterworth_filter(vigilance_at, None, 0.01, fs=fs, filter_type='lowpass')
    if norm:
        df = df.apply(zscore, axis=0)
    df.reset_index(inplace=True)
    df = df.rename(columns = {'index':'time'})
    return df


def load_subj_chang_bh(subj, scan, pc_ts, pc_ts_p, fs, norm=True, bp_filter=True):
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
    df['pc1_p'] = pc_ts_p[:,0]
    df['pc2_p'] = pc_ts_p[:,1]
    df['pc3_p'] = pc_ts_p[:,2]
    df['vigilance_ad'] = vigilance_ad
    df['vigilance_at'] = vigilance_at
    df['vigilance_ad_low'] =  butterworth_filter(vigilance_ad, None, 0.01, fs=fs, filter_type='lowpass')
    df['vigilance_at_low'] = butterworth_filter(vigilance_at, None, 0.01, fs=fs, filter_type='lowpass')
    if norm:
        df = df.apply(zscore, axis=0)
    df.reset_index(inplace=True)
    df = df.rename(columns = {'index':'time'})
    return df


def load_subj_hcp(subj, pc_ts, pc_ts_p, fs, norm=True):
    p_str = f'data/dataset_hcp/physio/proc1_physio/{subj}_physio.csv'
    csf_str = f'data/dataset_hcp/physio/proc1_physio/{subj}_csf.txt'
    df = pd.read_csv(p_str)
    csf = np.loadtxt(csf_str)
    df['csf'] = csf
    ## IMPORTANT! CSF signal has large amplitude spikes at start of scan - set first 10 time points to median
    breakpoint()
    df['csf'].iloc[:10] = df['csf'].median()
    df = df.apply(lambda x: butterworth_filter(x, 0.01, 0.1, fs=fs, filter_type='bandpass'), axis=0)
    df['pc1'] = pc_ts[:,0]*-1 # sign flip to keep consistent with Chang PCA
    df['pc2'] = pc_ts[:,1]*-1 # sign flip to keep consistent with Chang PCA
    df['pc3'] = pc_ts[:,2]
    df['pc1_p'] = pc_ts_p[:,0]*-1 # sign flip to keep consistent with Chang PCA
    df['pc2_p'] = pc_ts_p[:,1]*-1 # sign flip to keep consistent with Chang PCA
    df['pc3_p'] = pc_ts_p[:,2]
    if norm:
        df = df.apply(zscore, axis=0)
    df.reset_index(inplace=True)
    df = df.rename(columns = {'index':'time'})
    return df


def load_subj_nki(subj, pc_ts, pc_ts_p, fs, norm=True):
    p_str = f'data/dataset_nki/physio/proc1_physio/{subj}_task_breathhold_physio.csv'
    csf_str = f'data/dataset_nki/physio/proc1_physio/{subj}_task_breathhold_physio_csf.txt'
    df = pd.read_csv(p_str)
    csf = np.loadtxt(csf_str)
    df['csf'] = csf
    df = df.apply(lambda x: butterworth_filter(x, 0.01, 0.1, fs=fs, filter_type='bandpass'), axis=0)
    df['pc1'] = pc_ts[:,0]*-1 # sign flip to keep consistent with Chang PCA
    df['pc2'] = pc_ts[:,1]
    df['pc1_p'] = pc_ts_p[:,0]
    df['pc6_p'] = pc_ts_p[:,6]
    df['pc8_p'] = pc_ts_p[:,8]
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


def load_subj_spreng(subj, pc_ts, pc_ts_p, fs, norm=True):
    p_str = f'data/dataset_spreng/physio/proc1_physio/{subj}_ses-1_task-rest_physio.csv'
    csf_str = f'data/dataset_spreng/physio/proc1_physio/{subj}_ses-1_task-rest_physio_csf.txt'
    df = pd.read_csv(p_str)
    csf = np.loadtxt(csf_str)
    df['csf'] = csf
    df = df.apply(lambda x: butterworth_filter(x, 0.01, 0.1, fs=fs, filter_type='bandpass'), axis=0)
    df['pc1'] = pc_ts[:,0]*-1 # sign flip to keep consistent with Chang PCA
    df['pc2'] = pc_ts[:,1]
    df['pc3'] = pc_ts[:,2]
    df['pc1_p'] = pc_ts_p[:,0]
    df['pc2_p'] = pc_ts_p[:,1]*-1 # sign flip to keep consistent with Chang PCA
    df['pc3_p'] = pc_ts_p[:,2]*-1 # sign flip to keep consistent with Chang PCA
    if norm:
        df = df.apply(zscore, axis=0)
    df.reset_index(inplace=True)
    df = df.rename(columns = {'index':'time'})
    return df


def xcorr(x, y, maxlags=30, constrain='abs'):
    Nx = len(x)
    if Nx != len(y):
        raise ValueError('x and y must be equal length')
    c = np.correlate(x, y, mode=2)
    c /= np.sqrt(np.dot(x, x) * np.dot(y, y))
    if maxlags is None:
        maxlags = Nx - 1
    if maxlags >= Nx or maxlags < 1:
        raise ValueError('maglags must be None or strictly '
                         'positive < %d' % Nx)
    lags = np.arange(-maxlags, maxlags + 1)
    c = c[Nx - 1 - maxlags:Nx + maxlags]
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