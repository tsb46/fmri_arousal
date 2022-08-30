import numpy as np
import pandas as pd

from patsy import dmatrix
from scipy.stats import gamma, zscore
from scipy.signal import correlate
from sklearn.linear_model import LinearRegression, Ridge


def construct_design_matrix(model_formula, df, omit_intercept=True):
    # Modify formula to omit intercept (sklearn includes by default)
    if omit_intercept:
        model_formula = '0 + ' + model_formula
    design_mat = dmatrix(model_formula, df, 
                         return_type='dataframe')
    return design_mat


def construct_lag_splines(physio_var, p_nlags, n_n_lags, nknots):
    # Define model formula
    # Create lag sequence array (include lag of 0!)
    seq_lag = np.arange(-n_n_lags, p_nlags+1)
    # Create lag splines
    lag_splines = dmatrix("cr(x, df=nknots) - 1", 
                          {"x": seq_lag}, return_type='dataframe')
    # Create design matrix 
    basis_lag = np.zeros((len(physio_var), lag_splines.shape[1]))
    # Loop through lag bases and multiply column pairs
    lag_mat = pd.concat([physio_var.shift(l) for l in seq_lag], axis=1)
    for l in np.arange(lag_splines.shape[1]):
        basis_lag[:, l] = np.dot(lag_mat.values, lag_splines.iloc[:,l].values)
    return basis_lag, lag_splines


def construct_tensor_spline(v1, v2, nknots):
    # Define model formula
    spline_basis = dmatrix("te(cr(x, df=nknots), cr(y, df=nknots)) - 1", 
                          {"x": v1, 'y': v2}, 
                          return_type='dataframe')
    return spline_basis


def convolve_hrf(hrf, ts):
    n_drop = len(hrf) - 1
    convolved_events = np.convolve(ts, hrf)
    return convolved_events[:-n_drop]


def evaluate_tensor_spline(p1_eval, p2_eval, model, spline_basis):
    pred_maps = []
    for p1_e in p1_eval:
        p1_pred_maps = []
        for p2_e in p2_eval:
            # Create basis from model evaluation using previously defined design matrix (for model fit)
            pred_mat = dmatrix(spline_basis.design_info, {'x': p1_e, 'y': p2_e}, 
                                     return_type='dataframe')
            # Get predictions from model
            pred_map = model.predict(pred_mat)
            p1_pred_maps.append(pred_map)
        pred_maps.append(np.vstack(p1_pred_maps))
    return pred_maps


def get_hrf(t, tr, type):
    t_steps = np.arange(0, t, tr)
    if type == 'canonical':
        hrf = hrf_double_gamma(t_steps)
    elif type == 'rvt':
        hrf = hrf_rvt(t_steps)
    elif type == 'hr':
        hrf = hrf_hr(t_steps)
    return hrf


def hrf_double_gamma(t, dip=0.35):
    # http://www.jarrodmillman.com/rcsds/lectures/convolution_background.html
    gamma_peak = gamma.pdf(t, 6)
    gamma_under = gamma.pdf(t, 12)
    gamma_double = gamma_peak - dip * gamma_under
    return gamma_double / np.max(gamma_double) * 0.6


def hrf_hr(t):
    return (0.6*t**2.7)*np.exp(-t/1.6) - (16/np.sqrt(2*np.pi*9)*(np.exp((-1/2)*((t-12)**2)/9)))


def hrf_rvt(t):
    return (0.6*t**2.1)*np.exp(-t/1.6) - ((0.0023*t**3.54)*np.exp(-t/4.25))


def lag(arr, num, fill_value=0):
    # https://stackoverflow.com/questions/30399534/shift-elements-in-a-numpy-array
    result = np.empty_like(arr)
    if num > 0:
        result[:num] = fill_value
        result[num:] = arr[:-num]
    elif num < 0:
        result[num:] = fill_value
        result[:num] = arr[-num:]
    else:
        result[:] = arr
    return result


def lag_and_convolve_physio(physio_signals, physio_labels, n_subj, time_lag, convolve, tr):
    physio_sig_proc = []
    for subj_n in range(n_subj):
        subj_phys = []
        for p, p_label in zip(physio_signals, physio_labels):
            if convolve and ((p_label == 'rv') or (p_label == 'egg')):
                hrf = get_hrf(30, tr, 'rvt')
                subj_sig = convolve_hrf(hrf, p[subj_n])
            elif convolve and (p_label == 'hr'):
                hrf = get_hrf(30, tr, 'hr')
                subj_sig = convolve_hrf(hrf, p[subj_n])
            elif convolve:
                hrf = get_hrf(30, tr, 'canonical')
                subj_sig = convolve_hrf(hrf ,p[subj_n])
            else:
                subj_sig = p[subj_n]
            physio_sig_lag = lag(subj_sig, time_lag)
            subj_phys.append(physio_sig_lag)
        physio_sig_proc.append(np.stack(subj_phys, axis=1))
    return physio_sig_proc


def linear_regression(design_mat, func_data, labels=None, return_model=False, 
                      intercept=True, norm=True, reg=None):
    if norm:
        func_data = zscore(func_data)
        design_mat = zscore(design_mat)
    # Simple OLS - no mixed/multilevel model
    if reg == 'ridge':
        lin_reg = Ridge(fit_intercept=intercept, solver='lsqr')
    else:
        lin_reg = LinearRegression(fit_intercept=intercept)
    lin_reg.fit(design_mat, func_data)
    if return_model:
        return lin_reg
    elif labels is not None:
        betas = []
        for l in labels:
            l_indx = np.where(lin_reg.feature_names_in_ == l)[0][0]
            betas.append(lin_reg.coef_[:, l_indx])
        return betas
    else:
        raise Exception('labels must be provided, if fitted model is not returned')


def mask_voxels(func_data):
    # Identify and mask out voxels that have zero variance (i.e. are constant values, usually 0's)
    mask = np.where(np.std(func_data, axis=0) > 0)[0]
    func_data_masked = func_data[:, mask]
    beta_map = np.empty(func_data.shape[1])
    beta_map.fill(np.nan)
    return func_data_masked, mask, beta_map


def onsets_to_block(df, scan_len, tr):
    block_ts = []
    for onset, dur in zip(df.onset, df.duration):
        tr_event = int(np.floor(onset/tr))
        tr_dur = int(np.ceil(dur/tr))
        block_ts.append(np.arange(tr_event, tr_event+tr_dur))

    return block_ts


def xcorr(x, y, maxlags=30):
    # adjusted cross-correlation between two (equal-length) signals
    # https://www.statsmodels.org/dev/generated/statsmodels.tsa.stattools.ccf.html
    n = len(x)
    xo = x - x.mean()
    yo = y - y.mean()
    lags = np.arange(-maxlags, maxlags + 1)
    c_cov = correlate(xo, yo, "full", method='fft')/ n
    c_corr = c_cov/(np.std(x) * np.std(y))
    return lags, c_corr[n - 1 - maxlags:n + maxlags]




