import numpy as np
import pandas as pd

from patsy import dmatrix
from scipy.stats import gamma, zscore
from sklearn.linear_model import LinearRegression, Ridge


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


def create_quadratic_maps(beta_s1, beta_quad, 
                          plot_rng=[-4,-3,-2,-1,0,1,2,3,4]):
    # We ASSUME that the variables have been z-scored, so that a value of 1 represents 1 std. above the mean
    # Calculate the effect of V1 (beta_s1) at different levels of V2 (beta_s2)
    beta_simple = [beta_s1 + s1_val*beta_quad for s1_val in plot_rng]
    return np.array(beta_simple)


def get_hrf(t, tr, type):
    t_steps = np.arange(0, t, tr)
    if type == 'canonical':
        hrf = hrf_double_gamma(t_steps)
    elif type == 'rvt':
        hrf = hrf_rvt(t_steps)
    elif type == 'hr':
        hrf = hrf_hr(t_steps)
    return hrf


def get_interaction_map(interaction_map_str, design_cols, subj_beta_maps):
    if interaction_map_str not in design_cols:
            raise Exception('Interaction string specified for option -i does not match any string in model formula')
    else:
        i_index = design_cols.tolist().index(interaction_map_str)
    v1v2 = interaction_map_str
    v1_i, v2_i = parse_interaction_string(v1v2)
    # Remember, the beta maps are ordered according to  the order of the columns in the design_mat dataframe
    avg_beta_inter = np.mean([bmap[i_index] for bmap in subj_beta_maps], axis=0) 
    return avg_beta_inter, v1_i, v2_i


def get_quadratic_map(quadratic_map_str, design_cols, subj_beta_maps):
    design_cols_strip = [col.replace(' ', '') for col in design_cols]
    if quadratic_map_str not in design_cols_strip:
            raise Exception('Quadratic string specified for option -q does not match any string in model formula')
    else:
        i_index = design_cols_strip.index(quadratic_map_str)
    quad_term = quadratic_map_str
    simple_var = parse_quadratic_string(quad_term)
    # Remember, the beta maps are ordered according to the order of the columns in the design_mat dataframe
    avg_beta_quad = np.mean([bmap[i_index] for bmap in subj_beta_maps], axis=0) 
    return avg_beta_quad, simple_var


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
    block_ts = np.zeros(scan_len)
    for onset, dur in zip(df.onset, df.duration):
        tr_event = int(np.floor(onset/tr))
        tr_dur = int(np.ceil(dur/tr))
        block_ts[tr_event:(tr_event+tr_dur)] = 1
    return block_ts


def parse_quadratic_string(q_str):
    # parse interaction string from patsy model formula
    var = q_str.split('(')[1].split('**')[0]
    return var


def parse_interaction_string(i_str):
    # parse interaction string from patsy model formula
    v1, v2 = i_str.split(':')
    return v1, v2



