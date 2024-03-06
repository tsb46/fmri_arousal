import json
import numpy as np
import pandas as pd

from patsy import dmatrix
from scipy.stats import gamma, zscore
from scipy.signal import correlate
from sklearn.linear_model import LinearRegression, Ridge

# Global parameters
chang_block_len = 15 # length of chang task blocks (30s post auditory tone)
chang_trim = 14.7 # remove first 14.7s to account for first 7 volumes removed

# load params file (used to get TR and subject lists)
params = json.load(open('analysis_params.json', 'rb'))

# get nki and chang_bh TR
nki_tr = params['nki']['tr']
chang_tr = params['chang_bh']['tr'] # chang_bh and chang_cue have same TR


def construct_task_blocks(dataset, events, expand):
    # get TR indices of task blocks
    if dataset in ['chang_bh', 'chang_cue']:
        if dataset == 'chang_cue':
            events = events/1000 # events ms to s
        events = events - chang_trim # 14.7 secs
        # Convert event secs to TRs
        events_tr = events/chang_tr
        events_tr = np.round(events_tr).astype(int)
        # Construct blocks
        event_blocks = []
        for event in events_tr:
            event_blocks.append(np.arange(event, event+chang_block_len+expand))
        if dataset == 'chang_bh':
            event_blocks = [('breath', event_blocks)]
        elif dataset == 'chang_cue':
            event_blocks = [('cue', event_blocks)]
    elif dataset == 'nki':
        event_df = group_nki_blocks(events)
        group_func = {'onset': 'first', 'duration': 'sum'}
        # Construct 'deep breath' blocks
        breath_df = event_df.loc[event_df.breath]
        breath_df_grp = breath_df.groupby(['breath_indx']).agg(group_func)
        breath_blocks = onsets_to_block(breath_df_grp, nki_tr, expand)
        # Construct 'breath hold' blocks
        hold_df = event_df.loc[event_df.hold]
        hold_df_grp = hold_df.groupby(['hold_indx']).agg(group_func)
        hold_blocks = onsets_to_block(hold_df_grp, nki_tr, expand)
        event_blocks = [('rest-breath', breath_blocks), ('hold', hold_blocks)]
    return event_blocks


def group_nki_blocks(event_df, merge_rest_breath=True):
    # group consecutive trials of nki breathhold task
    event_df['hold'] = event_df.trial_type.str.startswith('H')
    event_df['hold_indx'] = (event_df.hold != event_df.hold.shift(1)).cumsum()
    if merge_rest_breath:
        rest_breath_block = ['R', 'G', 'Deep', 'In', 'Out']
        event_df['breath'] = event_df.trial_type.isin(rest_breath_block)
        event_df['breath_indx'] = (event_df.breath != event_df.breath.shift(1)).cumsum()
    else:
        breath_block = ['Deep', 'In', 'Out']
        event_df['breath'] = event_df.trial_type.isin(breath_block)
        event_df['breath_indx'] = (event_df.breath != event_df.breath.shift(1)).cumsum()
        rest_block = ['R', 'G']
        event_df['rest'] = event_df.trial_type.isin(rest_block)
        event_df['rest_indx'] = (event_df.rest != event_df.rest.shift(1)).cumsum()
    return event_df


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


def load_chang_compliance():
     # Load trial level compliance file
    compliance = pd.read_csv('data/dataset_chang_bh/compliance.csv')
    trial_cols = [f'trial{n+1}' for n in range(9)]
    # Convert to dict
    compliance_dict = {}
    for i, (subj, scan) in enumerate(zip(compliance.subject, compliance.scan)): 
        if scan < 10:
            scan_str = f'000{scan}'
        else:
            scan_str = f'00{scan}'
        compliance_dict[f'{subj}_{scan_str}'] = compliance.iloc[i, :][trial_cols].values
    return compliance_dict


def onsets_to_block(df, tr, expand=None):
    # convert blocks onsets and duration to TR indices of blocks
    block_ts = []
    for onset, dur in zip(df.onset, df.duration):
        tr_event = int(np.floor(onset/tr))
        tr_dur = int(np.ceil(dur/tr))
        if expand is not None:
            tr_dur += expand
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




