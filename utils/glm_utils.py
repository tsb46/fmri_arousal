import numpy as np
import pandas as pd

from patsy import dmatrix
from scipy.stats import gamma, zscore
from sklearn.linear_model import LinearRegression

def construct_design_matrix(model_formula, df, omit_intercept=True):
    # Modify formula to omit intercept (sklearn includes by default)
    if omit_intercept:
        model_formula = '0 + ' + model_formula
    design_mat = dmatrix(model_formula, df, 
                         return_type='dataframe')
    return design_mat


def convolve_hrf(hrf, ts):
    n_drop = len(hrf) - 1
    convolved_events = np.convolve(ts, hrf)
    return convolved_events[:-n_drop]


def create_interaction_maps(beta_s1, beta_i, 
                            plot_rng=[-4,-3,-2,-1,0,1,2,3,4]):
    # We ASSUME that the variables have been z-scored, so that a value of 1 represents 1 std. above the mean
    # Calculate the effect of V1 (beta_s1) at different levels of V2 (beta_s2)
    beta_simple = [beta_s1 + s2_val*beta_i for s2_val in plot_rng]
    return np.array(beta_simple)


def double_gamma_hrf(t, tr, dip=0.35):
    # http://www.jarrodmillman.com/rcsds/lectures/convolution_background.html
    n_steps = np.arange(0, t, tr)
    gamma_peak = gamma.pdf(n_steps, 6)
    gamma_under = gamma.pdf(n_steps, 12)
    gamma_double = gamma_peak - dip * gamma_under
    return gamma_double / np.max(gamma_double) * 0.6


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


def linear_regression(design_mat, func_data, labels, norm=True):
    if norm:
        func_data = zscore(func_data)
        design_mat = zscore(design_mat)
    # Simple OLS - no mixed/multilevel model
    lin_reg = LinearRegression()
    lin_reg.fit(design_mat, func_data)
    betas = []
    for l in labels:
        l_indx = np.where(lin_reg.feature_names_in_ == l)[0][0]
        betas.append(lin_reg.coef_[:, l_indx])
    return betas


def onsets_to_block(df, scan_len, tr):
    block_ts = np.zeros(scan_len)
    for onset, dur in zip(df.onset, df.duration):
        tr_event = int(np.round(onset/tr))
        tr_dur = int(np.round(dur))
        block_ts[tr_event:(tr_event+tr_dur)] = 1
    return block_ts


def parse_interaction_string(i_str):
    # parse interaction string from patsy model formula
    v1, v2 = i_str.split(':')
    return v1, v2
