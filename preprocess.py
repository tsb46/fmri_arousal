import argparse
import json
import mne
import nibabel as nb
import neurokit2 as nk
import numpy as np
import os
import pandas as pd
import shutil
import warnings


from itertools import repeat
from multiprocessing import Pool
from scipy.io import loadmat
from scipy.stats import zscore
from utils.fsl_utils import (
    apply_mask, apply_transform_mask, bet, concat_transform, 
    coregister, fast, first_vol, flirt, fnirt, 
    invert_transform, mcflirt, resample_func, robustfov, 
    reorient, slicetime, spatial_smooth, trim_vol, wm_thres, 
    warp_func
 )
from utils.load_write import (
    convert_2d, convert_4d,
    get_fp_base, load_subject_list
)
from utils.physio_utils import (
    extract_ecg_signals, extract_eeg_signals, 
    extract_gsr_signals, extract_ppg_signals, 
    extract_resp_signals
)
from utils.signal_utils import butterworth_filter, clip_spikes

# path to the analysis parameter .json file - contains dataset tr
params_fp='analysis_params.json'

# datasets
datasets = [
    'chang', 
    'chang_bh', 
    'chang_cue', 
    'nki', 
    'nki_rest', 
    'hcp', 
    'spreng', 
    'toronto', 
    'yale', 
    'natview'
]

# eeg channel selections for chang and natview datasets
chang_eeg_chan = ['P3', 'P4', 'Pz', 'O1', 'O2', 'Oz']
natview_eeg_chan = ['P3', 'P4', 'P7', 'P8', 'Pz', 'POz',
                    'P1', 'P2', 'PO3', 'PO4', 'P5', 'P6',
                    'PO7', 'PO8', 'O1', 'O2', 'Oz']


def anat_proc(fp, subj, output_dict, crop, bet_frac):
    # anatomical preprocessing pipeline
    fp = fp.format(subj)
    fp_base = get_fp_base(fp)
    # reorient structural to standard (just in case)
    fp_in = f"{output_dict['anat']['raw']}/{fp}"
    fp_reorient = f"{output_dict['anat']['reorient']}/{fp}"
    reorient(fp_in, fp_reorient)
    if crop:
        robustfov(fp_reorient, fp_reorient)
    # brain extraction
    fp_bet = f"{output_dict['anat']['bet']}/{fp}"
    bet(fp_reorient, fp_bet, bet_frac)
    # fast tissue segmentation - return white matter segment
    fp_fast = f"{output_dict['anat']['fast']}/{fp}"
    fp_wm_base = fast(fp_bet, fp_fast)
    # threshold white matter segmentation for coregistration
    fp_wm_thres = f"{fp_wm_base}_thres.nii.gz"
    fp_wm = f'{fp_wm_base}.nii.gz'
    wm_thres(fp_wm, fp_wm_thres)
    # FLIRT affine registration to MNI template
    fp_flirt = f"{output_dict['anat']['flirt']}/{fp}"
    fp_flirt_mat = f"{output_dict['anat']['flirt']}/{fp_base}_flirt.mat"
    flirt(fp_bet, fp_flirt, fp_flirt_mat)
    # FNIRT non-linear registration to MNI
    fp_fnirt = f"{output_dict['anat']['fnirt']}/{fp}"
    fp_fnirt_coef = \
    f"{output_dict['anat']['fnirt']}/{fp_base}_fn_coef.nii.gz"
    fnirt(fp_reorient, fp_flirt_mat, fp_fnirt, fp_fnirt_coef)
    # return head, brain, white matter seg,
    # flirt affine mat and fnirt warp for functional preproc

    anat_out = {
        'reorient': fp_reorient,
        'bet': fp_bet,
        'wm': fp_wm_thres,
        'flirt_mat': fp_flirt_mat,
        'fnirt_coef': fp_fnirt_coef
    }
    return (subj, anat_out)


def bandpass_filter(fp, fp_out, tr, mask, cut_low=0.01, cut_high=0.1):
    fs = 1/tr
    # custom bandpass filter function
    img = nb.load(fp)
    data = img.get_fdata()
    data_2d = convert_2d(mask, data)
    data_filt = butterworth_filter(data_2d, cut_low, cut_high, fs, 'bandpass')
    # zscore along time dimension
    data_filt = zscore(data_filt)
    data_filt_4d = convert_4d(mask, data_filt)
    img_out = nb.Nifti1Image(data_filt_4d, img.affine,
                             img.header)
    img_out.to_filename(fp_out)


def create_directories(dataset, p_type, eeg, slice_time=None, trim=None):
    # create directories for preprocessing
    if (p_type == 'full') | (p_type == 'multiecho'):
        output_dict = {
            'anat': {
                'raw': f'data/dataset_{dataset}/anat/raw',
                'reorient': f'data/dataset_{dataset}/anat/proc1_reorient',
                'bet': f'data/dataset_{dataset}/anat/proc2_bet',
                'fast': f'data/dataset_{dataset}/anat/proc3_fast',
                'flirt': f'data/dataset_{dataset}/anat/proc4_flirt',
                'fnirt': f'data/dataset_{dataset}/anat/proc5_fnirt'
            }
        }
        if p_type == 'full':
            output_dict['func'] = {
                'raw': f'data/dataset_{dataset}/func/raw',
                'trim': f'data/dataset_{dataset}/func/procA_trim',
                'slicetime': f'data/dataset_{dataset}/func/procB_slicetime',
                'mcflirt': f'data/dataset_{dataset}/func/proc1_mcflirt',
                'func2struct': f'data/dataset_{dataset}/func/proc2_func2struct',
                'standard': f'data/dataset_{dataset}/func/proc3_standard',
                'smooth': f'data/dataset_{dataset}/func/proc4_mask_smooth',
                'bandpass': f'data/dataset_{dataset}/func/proc5_bandpass'
             }
        elif p_type == 'multiecho':
            output_dict['func'] = {
                'raw': f'data/dataset_{dataset}/func/raw',
                'trim': f'data/dataset_{dataset}/func/procA_trim',
                'slicetime': f'data/dataset_{dataset}/func/procB_slicetime',
                'mcflirt': f'data/dataset_{dataset}/func/proc1_mcflirt',
                'tedana': f'data/dataset_{dataset}/func/proc2_tedana',
                'func2struct': f'data/dataset_{dataset}/func/proc3_func2struct',
                'standard': f'data/dataset_{dataset}/func/proc4_standard',
                'smooth': f'data/dataset_{dataset}/func/proc5_mask_smooth',
                'bandpass': f'data/dataset_{dataset}/func/proc6_bandpass'
             }
    elif p_type == 'minimal':
        output_dict = {
            'func': {
                'raw': f'data/dataset_{dataset}/func/raw',
                'resample': f'data/dataset_{dataset}/func/proc1_resample',
                'smooth': f'data/dataset_{dataset}/func/proc2_smooth_mask',
                'bandpass': f'data/dataset_{dataset}/func/proc3_bandpass'
             }
        }

    # set preprocessed physio directory
    output_dict['physio'] = {
        'raw': f'data/dataset_{dataset}/physio/raw',
        'proc': f'data/dataset_{dataset}/physio/proc1_physio'
    }
    if eeg:
        output_dict['eeg'] = {
            'raw': f'data/dataset_{dataset}/eeg/raw',
            'proc': f'data/dataset_{dataset}/eeg/proc1_eeg'
        }
    for f_type in output_dict:
        for out in output_dict[f_type]:
            os.makedirs(output_dict[f_type][out], exist_ok=True)
    return output_dict


def eeglab_natview_preprocess(fp, output_dir):
    """
    load and preprocess raw natview EEG data using EEGLAB
    need to use the system command to run matlab script in 'utils' folder
    ensure that matlab is callable from the command line
    the eeg/raw directory for natview is added to path so just need to 
    supply filename to function. Handling file paths is a mess w/ matlab
    from the command line - couple hacks to make things work.
    """
    # import MATLAB Engine API for Python
    import matlab.engine
    import io
    # start matlab to python engine and add path to 'eeglab_preprocess.m'
    eng = matlab.engine.start_matlab()
    eng.addpath('utils', nargout=0)
    eng.feval('eeglab_preprocess', fp, output_dir, nargout=0, stdout=io.StringIO())
    eng.quit()


def epi_mask(fp_in, fp_mask):
    # create mask from EPI for tedana pipeline
    from nilearn.masking import compute_epi_mask
    nii = nb.load(fp_in)
    nii_mask = compute_epi_mask(nii, lower_cutoff=0.01, upper_cutoff=0.6, 
                                opening=True)
    nb.save(nii_mask, fp_mask)


def extract_physio(ts, phys_label, sf):
    # extract features from physiological signals using Neurokit
    # extraction functions return pandas dataframe
    if phys_label == 'ecg':
        phys_ts = extract_ecg_signals(ts, sf)
    elif phys_label == 'eeg':
        phys_ts = extract_eeg_signals(ts, sf)
    elif phys_label == 'gsr':
        phys_ts = extract_gsr_signals(ts, sf)
    elif phys_label == 'resp':
        phys_ts = extract_resp_signals(ts, sf)
    elif phys_label == 'ppg':
        phys_ts = extract_ppg_signals(ts, sf)
    elif phys_label == 'pupil':
        # not implemented, pupil signals are already extracted
        phys_ts = pd.DataFrame({'PUPIL': ts})

    return phys_ts


def func_full_proc(fp, subj, scan, anat_out_dict, output_dict, tr,
                   mask, slice_timing=None, trim=None):
    # full preprocessing pipeline starting with raw functional
    fp = fp.format(subj, scan)
    fp_base = get_fp_base(fp)
    # get anatomical filepaths
    anat_out = anat_out_dict[subj]
    # get starting raw functional scan
    fp_in = f"{output_dict['func']['raw']}/{fp}"
    # if specified, trim first N volumes
    if trim is not None:
        fp_trim = f"{output_dict['func']['trim']}/{fp}"
        trim_vol(fp_in, fp_trim, trim)
        fp_in = fp_trim
    # if specified, slicetime correct
    if slice_timing is not None:
        fp_slicetime = f"{output_dict['func']['slicetime']}/{fp}"
        slicetime(fp_in, fp_slicetime, slice_timing, tr)
        fp_in = fp_slicetime
    # apply mcflirt motion correction
    fp_mcflirt = f"{output_dict['func']['mcflirt']}/{fp}"
    fp_mean = mcflirt(fp_in, fp_mcflirt)
    # epi coregistration
    fp_coreg = f"{output_dict['func']['func2struct']}/{fp}"
    fp_func2struct = coregister(
        fp_mean, anat_out['reorient'], anat_out['bet'], anat_out['wm'], fp_coreg
    )
    # apply warp to get functional to MNI
    fp_warp = f"{output_dict['func']['standard']}/{fp}"
    warp_func(fp_mcflirt, fp_func2struct, anat_out['fnirt_coef'], fp_warp, mask)
    # apply func minimal preprocessing pipeline
    func_minimal_proc(fp, subj, scan, output_dict, tr, mask,
                      resample=False, smooth=True)


def func_minimal_proc(fp, subj, scan, output_dict, tr, mask,
                      resample=True, smooth=True):
    # minimal preprocessing pipeline starting with MNI-registred and 
    # preprocessed functional data
    fp = fp.format(subj, scan)
    fp_base = get_fp_base(fp)
    if resample:
        # resample functional to 3mm MNI space
        fp_in = f"{output_dict['func']['raw']}/{fp}"
        fp_resample = f"{output_dict['func']['resample']}/{fp}"
        resample_func(fp_in, fp_resample, mask)
        fp_in = fp_resample
    else:
        fp_in = f"{output_dict['func']['standard']}/{fp}"
    # mask and, if specified, smooth (5mm FWHM)
    fp_smooth = f"{output_dict['func']['smooth']}/{fp}"
    if smooth:
        spatial_smooth(fp_in, fp_smooth)
        fp_in = fp_smooth
    apply_mask(fp_in, fp_smooth, mask)
    # bandpass filter
    # Load mask
    mask_bin = nb.load(mask).get_fdata() > 0
    fp_bandpass = f"{output_dict['func']['bandpass']}/{fp}"
    bandpass_filter(fp_smooth, fp_bandpass, tr, mask_bin)


def func_me_proc(fp_me, echo_times, subj, scan, anat_out_dict, output_dict,
                  tr, mask, slice_timing=None, trim=None):
    # full multiecho preprocessing pipeline starting with raw functional
    # get combined functional name
    fp_func = fp_me['func'].format(subj, scan)
    fp_func_base = get_fp_base(fp_func)
    # get functional names of echos (assume index starts at 1)
    fp_echo = []
    for i in range(len(echo_times)):
        if scan is None:
            fp = fp_me['echo'].format(subj, i+1)
        else:
            fp = fp_me['echo'].format(subj, scan, i+1)
        fp_echo.append(fp)
    # get anatomical filepaths
    anat_out = anat_out_dict[subj]
    # get full path to raw functional files
    fp_in_dir = output_dict['func']['raw']
    # if specified, trim first N volumes
    if trim is not None:
        for fp in fp_echo:
            fp_in = f"{fp_in_dir}/{fp}"
            fp_trim = f"{output_dict['func']['trim']}/{fp}"
            trim_vol(fp_in, fp_trim, trim)
        fp_in_dir = output_dict['func']['trim']
    # if specified, slicetime correct
    if slice_timing is not None:
        for fp in fp_echo:
            fp_in = f"{fp_in_dir}/{fp}"
            fp_slicetime = f"{output_dict['func']['slicetime']}/{fp}"
            slicetime(fp_in, fp_slicetime, slice_timing, tr)
        fp_in_dir = output_dict['func']['slicetime']
    # apply mcflirt motion correction on first echo
    fp_in = f'{fp_in_dir}/{fp_echo[0]}'
    fp_mcflirt = f"{output_dict['func']['mcflirt']}/{fp_echo[0]}"
    fp_echo1_mean = mcflirt(fp_in, fp_mcflirt, save_mats=True)
    # apply mcflirt transform to rest of echos
    for fp in fp_echo[1:]:
        fp_in = f"{fp_in_dir}/{fp}"
        fp_mcflirt = f"{output_dict['func']['mcflirt']}/{fp}"
        fp_mcflirt_mat = f"{output_dict['func']['mcflirt']}/{fp_echo[0]}.mat"
        # applyxfm4d is not available in nipype, run from terminal
        os.system(
            f'applyxfm4D {fp_in} {fp_echo1_mean} {fp_mcflirt} {fp_mcflirt_mat} -fourdigit'
        )
    # epi coregistration on mean func
    fp_coreg = f"{output_dict['func']['func2struct']}/{fp_func}"
    fp_func2struct = coregister(
        fp_echo1_mean, anat_out['reorient'], anat_out['bet'], anat_out['wm'], fp_coreg
    )
    fp_func2struct_mat = f"{get_fp_base(fp_func2struct)}_init.mat"
    # tedana denoising
    # create directory for tedana output
    tedana_out_dir = f"{output_dict['func']['tedana']}/{fp_func_base}"
    os.makedirs(tedana_out_dir, exist_ok=True)
    # put mask from BET into to functional space
    fp_coreg_mat = f"{get_fp_base(fp_func2struct)}_init_invert.mat"
    invert_transform(fp_func2struct_mat, fp_coreg_mat)
    fp_bet_base = f"{os.path.basename(get_fp_base(anat_out['bet']))}_mask"
    fp_bet_mask = f"{output_dict['anat']['bet']}/{fp_bet_base}.nii.gz"
    fp_bet_func = f"{tedana_out_dir}/{fp_bet_base}_func.nii.gz"
    apply_transform_mask(fp_bet_mask, fp_bet_func, fp_echo1_mean, fp_coreg_mat)    
    # execute tedana pipeline
    fps_denoise = []
    for fp in fp_echo:
        fps_denoise.append(f"{output_dict['func']['mcflirt']}/{fp}")
    tedana_denoise(fps_denoise, echo_times, fp_bet_func, tedana_out_dir, fp_func_base)
    # define path to denoised output
    fp_tedana = f"{tedana_out_dir}/{fp_func_base}_desc-optcomDenoised_bold.nii.gz"
    # apply warp to get functional to MNI
    fp_warp = f"{output_dict['func']['standard']}/{fp_func}"
    warp_func(fp_tedana, fp_func2struct, anat_out['fnirt_coef'], fp_warp, mask)
    # apply func minimal preprocessing pipeline
    func_minimal_proc(fp_func, subj, scan, output_dict, tr, mask,
                      resample=False, smooth=True)


def get_anat_fp(fp, subj_list, output_dict):
    # get file paths to anatomical outputs necessary
    # for preprocessing functional scans
    anat_out = {}
    for subj in subj_list:
        fp_s = str(fp)
        fp_s = fp_s.format(subj)
        fp_base = get_fp_base(fp_s)
        fp_fast_base = get_fp_base(f"{output_dict['anat']['fast']}/{fp_s}")
        anat_out_subj = {
            'reorient': f"{output_dict['anat']['reorient']}/{fp_s}",
            'bet': f"{output_dict['anat']['bet']}/{fp_s}",
            'wm': f'{fp_fast_base}_pve_2_thres.nii.gz',
            'flirt_mat': f"{output_dict['anat']['flirt']}/{fp_base}_flirt.mat",
            'fnirt_coef': f"{output_dict['anat']['fnirt']}/{fp_base}_fn_coef.nii.gz"
        }
        anat_out[subj] = anat_out_subj
    return anat_out


def get_fp(dataset):
    # set filepath templates for raw data
    if dataset == 'chang':
        func = {
            'echo': 'sub_00{0}-mr_{1}-ecr_echo{2}.nii.gz',
            'func': 'sub_00{0}-mr_{1}-ecr.nii.gz',
        }
        anat = 'sub_00{0}-mprage.nii.gz'
        physio = {
            'physio': 'sub_00{0}-mr_{1}-ecr_echo1_physOUT.mat',
            'eeg': 'sub_00{0}-mr_{1}_eeg_pp.mat',
            'out': 'sub_00{0}-mr_{1}-ecr_physio'
        }
    elif dataset == 'chang_bh':
        func = {
            'echo': 'sub_00{0}-mr_{1}-adb_echo{2}.nii.gz',
            'func': 'sub_00{0}-mr_{1}-adb.nii.gz',
        }
        anat = 'sub_00{0}-mprage.nii.gz'
        physio = {
            'physio': 'sub_00{0}-mr_{1}-adb_echo1_physOUT.mat',
            'eeg': 'sub_00{0}-mr_{1}-adb_echo1_EEG_pp.mat',
            'out': 'sub_00{0}-mr_{1}-adb_physio'
        }
    elif dataset == 'chang_cue':
        func = {
            'echo': 'sub_{0}-mr_{1}-ect_echo{2}.nii.gz',
            'func': 'sub_{0}-mr_{1}-ect.nii.gz',
        }
        anat = 'sub_{0}-mprage.nii.gz'
        physio = {
            'physio': 'sub_{0}-mr_{1}-ect_echo1_physOUT.mat',
            'eeg': 'sub_{0}-mr_{1}-ect_echo1_EEG_pp.mat',
            'out': 'sub_{0}-mr_{1}-ect_physio'
        }
    elif dataset == 'hcp':
        func = '{0}_{1}_hp2000_clean.nii.gz'
        anat = None
        physio = {
          'physio': '{0}_{1}_physio.txt',
          'out': '{0}_{1}_physio'
        }
    elif dataset == 'natview':
        func='sub-{0}_ses-0{1}_task-rest_bold.nii.gz'
        anat = 'sub-{0}_T1w1_denoise.nii.gz'
        physio = {
            'eeg': 'sub-{0}_ses-0{1}_task-rest_eeg.set',
            'eye': 'sub-{0}_ses-0{1}_task-rest_recording-eyetracking_physio.tsv.gz',
            'resp': 'sub-{0}_ses-0{1}_task-rest_recording-respiratory_physio.tsv.gz',
            'ecg': 'sub-{0}_ses-0{1}_task-rest_ecg.set',
            'out': 'sub-{0}_ses-0{1}_task-rest_physio'
        }
    elif dataset == 'nki': 
        func = 'sub-{0}_ses-{1}_task-BREATHHOLD_acq-1400_bold.nii.gz'
        anat = 'sub-{0}_T1w.nii.gz'
        physio = {
            'physio': 'sub-{0}_ses-{1}_task-BREATHHOLD_acq-1400_physio.tsv.gz',
            'out': 'sub-{0}_ses-{1}_task-BREATHHOLD_acq-1400_physio'
        }
    elif dataset == 'nki_rest': 
        func = 'sub-{0}_ses-{1}_task-rest_acq-1400_bold.nii.gz'
        anat = 'sub-{0}_T1w.nii.gz'
        physio = {
            'physio': 'sub-{0}_ses-{1}_task-rest_acq-1400_physio.tsv.gz',
            'out': 'sub-{0}_ses-{1}_task-rest_acq-1400_physio'
        }
    elif dataset == 'spreng':
        func = {
            'echo': '{0}_{1}_task-rest_echo-{2}_bold.nii.gz',
            'func': '{0}_{1}_task-rest_bold.nii.gz',
        }
        anat = '{0}_T1w.nii.gz'
        physio = {
            'physio': '{0}_{1}_task-rest_physio.tsv.gz',
            'out': '{0}_{1}_task-rest_physio'
        }
    elif dataset == 'toronto': 
        func = '{0}_rest_bold_{1}.nii.gz'
        anat = '{0}_T1w.nii.gz'
        physio = {
            'physio': '{0}_rest_bold_{1}.mat',
            'out': '{0}_rest_bold_{1}_physio'
        }
    elif dataset == 'yale':
        func = '{0}_task-rest_run-0{1}_bold.nii.gz'
        anat = '{0}_T1w.nii.gz'
        physio = {
            'physio': '{0}_task-rest_run-0{1}_et.tsv',
            'out': '{0}_task-rest_run-0{1}_physio'
        }
    return func, anat, physio


def load_proc_eeg(fp, dataset, resample, trim, no_eeglab, eeglab_dir=None):
    # load eeg data and return MNE object
    if dataset in ['chang', 'chang_bh', 'chang_cue']:
        eeg_mat = loadmat(fp, squeeze_me=True)
        sf = eeg_mat['EEG']['srate'].item()
        n_chan = eeg_mat['EEG']['nbchan'].item()
        chan_labels = [chan[0] for chan in eeg_mat['EEG']['chanlocs'].item()]
        info = mne.create_info(chan_labels, ch_types='eeg', sfreq=sf)
        data = np.vstack(eeg_mat['EEG']['data'].item())
        eeg_mne = mne.io.RawArray(data, info, verbose=False)
        # Resample from 250 to 100Hz
        eeg_mne = eeg_mne.resample(resample)
        # Trim off first 14.7s to align w/ functional (first 7 TRs were trimmed from functional)
        eeg_mne.crop(tmin=trim)
        # select channels
        eeg_mne.pick(chang_eeg_chan)
        return eeg_mne, resample

    elif dataset == 'natview':
        if not no_eeglab:
            # Run EEGLAB Preprocessing script
            eeglab_natview_preprocess(fp, eeglab_dir)
        # load in preprocessed EEG and trim to align with functional
        eeg_base = os.path.basename(fp)
        eeg_in = f'{eeglab_dir}/{eeg_base}'
        eeg_mat = loadmat(eeg_in, squeeze_me=True)
        # 'R128' is the trigger label
        trigger_indx = [e['latency'] for e in eeg_mat['urevent'] 
                        if e['type'] == 'R128']
        start_indx = int(round(trigger_indx[0])) # first trigger is the first TR
        # we remove the last two triggers to align with functions
        end_indx = int(round(trigger_indx[-2]))
        eeg_data = eeg_mat['data'][:, start_indx:end_indx]
        # Create MNE object
        chan_labels = eeg_mat['chanlocs']['labels'].tolist()
        sf = eeg_mat['srate']
        info = mne.create_info(chan_labels, ch_types='eeg', sfreq=sf)
        eeg_mne = mne.io.RawArray(eeg_data, info, verbose=False)
        # select channels
        eeg_mne.pick(natview_eeg_chan)
        # package start and end indx for trimming ECG data that was 
        # contained in EEG data and loaded later in pipeline
        trim_indx = (start_indx, end_indx)
        return eeg_mne, sf, trim_indx


def load_physio(fp, subj, scan, dataset, output_dict, resample, trim, no_eeglab):
    # load physiological signals for each dataset and return dictionary with 
    # physio labels as keys.
    # each dataset has unique file formats that must be accounted for  
    # load physio for chang datasets  
    if dataset in ['chang', 'chang_bh', 'chang_cue']:
        fp_p = fp['physio'].format(subj, scan)
        fp_eeg = fp['eeg'].format(subj, scan)
        fp_p_in = f"{output_dict['physio']['raw']}/{fp_p}"
        fp_eeg_in = f"{output_dict['eeg']['raw']}/{fp_eeg}"
        # Load eeg .mat file
        physio_raw = loadmat(fp_p_in, squeeze_me=True)
        sf = 1/physio_raw['OUT_p']['dt_phys'].item()
        # Pull physio data into dict
        ppg = physio_raw['OUT_p']['card_dat'].item()
        resp = physio_raw['OUT_p']['resp'].item()['wave'].item()
        # Load EEG for chang datasets
        eeg, sf_eeg = load_proc_eeg(fp_eeg_in, dataset, resample, trim, no_eeglab)
        # write out eeg object
        fp_eeg_out = get_fp_base(fp_eeg)
        eeg.save(f"{output_dict['eeg']['proc']}/{fp_eeg_out}.raw.fif", overwrite=True)
        physio = {'ppg': ppg, 'resp': resp, 'eeg': eeg}
        sf_dict = {'ppg': sf, 'resp': sf, 'eeg': sf_eeg}
    # load hcp physio
    elif dataset == 'hcp':
        sf = 400 # hcp physio sampling frequency
        fp_p = fp['physio'].format(subj, scan)
        fp_p_in = f"{output_dict['physio']['raw']}/{fp_p}"
        physio_raw = np.loadtxt(fp_p_in)
        physio = {
            'resp': physio_raw[:,1], 
            'ppg': physio_raw[:,2]
        }
        sf_dict = {'resp': sf, 'ppg': sf}
    # load natview physio
    elif dataset == 'natview':
        # get file paths
        fp_eye = fp['eye'].format(subj, scan)
        fp_resp = fp['resp'].format(subj, scan)
        fp_eeg = fp['eeg'].format(subj, scan)
        fp_ecg = fp['ecg'].format(subj, scan)
        fp_eye_in = f"{output_dict['physio']['raw']}/{fp_eye}"
        fp_resp_in = f"{output_dict['physio']['raw']}/{fp_resp}"
        fp_eeg_in = f"{output_dict['eeg']['raw']}/{fp_eeg}"
        # ECG file doesn't exist yet, created in EEGLAB preprocessing
        fp_ecg_in = f"{output_dict['eeg']['proc']}/{fp_ecg}"
        # Load eye tracking file and json metadata
        eye = pd.read_csv(fp_eye_in, compression='gzip', delimiter='\t', header=None)
        eye_json = json.load(open(f'{get_fp_base(fp_eye_in)}.json', 'rb'))
        eye.columns = eye_json['Columns']
        sf_eye = eye_json['SamplingFrequency']
        # eye signals trimming to align w/ functional
        # create a indicator of start/end of task - and trim to start of task
        eye['task_span'] = eye.Task_Start_End_Trigger.cumsum() 
        eye_trim = eye.loc[eye.task_span >= 1].copy()
        # Create a running sum of volume triggers.
        """
        Important Note: 
        the eye recordings have 285 triggers (three short of the functional - 
        the functional has 288 volumes). One trigger at the beginning is missing 
        because it occurs 50 milliseconds before the task start trigger, so we can assume 
        that the start of the task is aligned with the first (missing) trigger (i.e. first 
        volume of the functional). Thus, we can trim all signals before the start of the task 
        trigger in the eye recordings. We trim all signals in the eye recordings after the 
        last trigger (285) recorded. Because there are two (missing) fMRI triggers after the 
        final trigger recorded in the eye recordings, we need to remove the last two volumes
        from the fMRI and EEG recordings. 
        """
        eye_trim['trigger_sum'] = eye_trim.fMRI_Volume_Trigger.cumsum()
        last_indx = eye_trim.loc[eye_trim.trigger_sum == 285].index[0]+1
        pupil_trim = eye_trim.loc[:last_indx]['Pupil_Area'].copy()
        # Load respiratory data
        """
        Important Note:
        We have little metadata on the respiratory belt signals. Still waiting to hear
        back from NATVIEW team on alignment with the rest of the recordings. For now,
        assume that is already aligned to the start and end of the functional scan. Thus,
        we should shave off 4.2 sec (2 TRs) from the end to match the trimming off 
        the end of the functional scan by 2TRs. 
        """
        resp = pd.read_csv(fp_resp_in, header=None, delimiter='\t', compression='gzip')
        resp_json = json.load(open(f'{get_fp_base(fp_resp_in)}.json', 'rb'))
        resp.columns = resp_json['Columns']
        sf_resp = resp_json['SamplingFrequency']
        # Trim off 4.2 sec (2 TRs) to align with functional (TR = 2.1)
        indx = np.floor((2*2.1)/(1/sf_resp)).astype(int)
        resp_trim = resp.iloc[:-indx]['respiratory']
        sf_resp = resp_json['SamplingFrequency']
        # package up physio signals and their sampling frequency
        physio = {'pupil': pupil_trim, 'resp': resp_trim}
        sf_dict = {'pupil': sf_eye, 'resp': sf_resp}
        # If the user specifies no_eeglab, do not run eeglab matlab script.
        if no_eeglab:
            # check whether the eeglab preprocessing script has been run
            # separately without matlabengine, if so move forward with
            # further eeg and ecg preprocessing
            fp_eeg_in_proc = f"{output_dict['eeg']['proc']}/{fp_eeg}"
            if os.path.isfile(fp_eeg_in_proc):
                eeg_p_check = True
            else:
                eeg_p_check = False
                warnings.warn("""
                    EEGLAB preprocessing was chosen not to run, and a preprocessed EEG 
                    file doesn't look to exist. No EEG and ECG signals will be extracted.
                """)
        else:
            eeg_p_check = True
        # if preprocessed eeg object exists, move forward with further preprocessing
        if eeg_p_check:
            # Load eeg data and, if specified, preprocess with EEGLAB (need MATLAB)
            eeg, sf_eeg, trim_indx = load_proc_eeg(
                fp_eeg_in, dataset, resample, trim, no_eeglab, output_dict['eeg']['proc']
            )
            # write out eeg object
            fp_eeg_out = get_fp_base(fp_eeg)
            eeg.save(f"{output_dict['eeg']['proc']}/{fp_eeg_out}.raw.fif", overwrite=True)
            # load in ECG data
            ecg = loadmat(fp_ecg_in, squeeze_me=True)
            ecg_trim = ecg['data'][trim_indx[0]:trim_indx[1]][:, np.newaxis]
            # insert into dictionary
            physio['ecg'] = ecg_trim; physio['eeg'] = eeg;
            sf_dict['ecg'] = sf_eeg; sf_dict['eeg'] = sf_eeg 
 
        return physio, sf_dict
    # load nki or spreng physio
    elif (dataset == 'nki') | (dataset == 'nki_rest') | (dataset == 'spreng'):
        # get file path and load
        fp_p = fp['physio'].format(subj, scan)
        fp_p_in = f"{output_dict['physio']['raw']}/{fp_p}"
        physio_df = pd.read_csv(fp_p_in, compression='gzip', sep='\t', header=None)
        physio_df = physio_df.dropna(axis=1, how='all')
        # metadata jsons for each dataset need separate handling
        fp_json = get_fp_base(fp_p_in)
        # load json and set columns
        physio_json = json.load(open(f'{fp_json}.json'))
        physio_df.columns = physio_json['Columns']
        # If nki dataset, trim off extra signals after last trigger
        if (dataset == 'nki') | (dataset == 'nki_rest'):
            last_trigger = physio_df.loc[physio_df['trigger'] == 5].index[-1]
            physio_df = physio_df.iloc[:last_trigger, :].copy()
        physio = {
            'ppg': physio_df['cardiac'].values, 
            'resp': physio_df['respiratory'].values
        }
        sf = physio_json["SamplingFrequency"]
        sf_dict = {'ppg': sf, 'resp': sf}
        if (dataset == 'nki') | (dataset == 'nki_rest'):
            physio['gsr'] = physio_df['gsr'].values
            sf_dict['gsr'] = sf
    # load toronto physio
    elif dataset == 'toronto':
        sf = 200 # toronto physio sampling frequency
        resp_label = 'RSP: LPF'
        fp_p = fp['physio'].format(subj, scan)
        fp_p_in = f"{output_dict['physio']['raw']}/{fp_p}"
        # Load eeg .mat file
        physio_raw = loadmat(fp_p_in, squeeze_me=True)
        physio_df = pd.DataFrame(physio_raw['data'], 
                                 columns=physio_raw['labels'])
        # trim to start and end of trigger
        trig_start = physio_df.loc[physio_df['TRIG'] == 5].index[0]
        trig_end = physio_df.loc[physio_df['TRIG'] == 5].index[-1]
        physio = {
            'resp': physio_df.loc[trig_start:trig_end][resp_label].values
        }
        sf_dict = {'resp': sf}
    # load yale physio
    elif dataset == 'yale':
        sf = 1 # already resampled to functional scan TR
        fp_p = fp['physio'].format(subj, scan)
        fp_p_in = f"{output_dict['physio']['raw']}/{fp_p}"
        physio_df = pd.read_csv(fp_p_in, sep='\t', header=None)
        physio = {
            'pupil': physio_df.iloc[:,0].values
        }
        sf_dict = {'pupil': sf}

    # loop through physio signals and trim/resample (if specified)
    for p in physio.keys():
        if p != 'eeg':
            # some dataset physio signals need to be trimmed to align w/ functional
            if trim is not None:
                sf_trim = sf_dict[p]
                trim_n = int(sf_trim*trim)
                physio[p] = physio[p][trim_n:]
            # to ease computational burden, some high-frequency 
            # physio signals are downsampled before pre-processing
            if resample is not None:
                sf_resamp = sf_dict[p]
                physio[p] = nk.signal_resample(
                           physio[p], 
                           sampling_rate=sf_resamp, 
                           desired_sampling_rate=resample, 
                           method='FFT'
                        )
                # set resample freq as new freq (sf)
                sf_dict[p] = resample                
        
    return physio, sf_dict


def physio_proc(fp, subj, scan, dataset, fp_func, 
                output_dict, resample, trim, 
                resample_to_func, no_eeglab):
    # preprocess physio signals - save out unfiltered and 
    # bandpass filtered preprocessed signals
    # load physio signals
    physio, sf_dict = load_physio(fp, subj, scan, dataset, output_dict, 
                                  resample, trim, no_eeglab)
    # get n of time points of functional scan for aligning w/ physio
    # need to load in subj functional nifti header
    fp_func_in = f"{output_dict['func']['bandpass']}/{fp_func.format(subj, scan)}"
    func_n = nb.load(fp_func_in).shape[-1]
    # loop through physio signals and extract, clip, filter and resample
    physio_out = []
    physio_out_filt = []
    for p in physio:
        p_df = extract_physio(physio[p], p, sf_dict[p])
        # clip potential spikes - abs(z) >= 5
        p_df = p_df.apply(clip_spikes, axis=0)
        # Bandpass filter signals (0.01-0.1Hz)
        p_df_filt = p_df.apply(
           butterworth_filter, lowcut=0.01, highcut=0.1, fs=sf_dict[p], 
           filter_type='bandpass', axis=0
        )
        # whether to resample to functional - some physio signals are already 
        # preprocessed and resampled to functional
        if resample_to_func:
            # Resample unfiltered physio signals to length of functional scan
            p_df_r = p_df.apply(
                nk.signal_resample, desired_length=func_n, method='FFT', axis=0
            )
            # Resample filtered physio signals to length of functional scan
            # bandpass filter
            p_df_filt_r = p_df_filt.apply(
                nk.signal_resample, desired_length=func_n, method='FFT', axis=0
            )
            physio_out.append(p_df_r)
            physio_out_filt.append(p_df_filt_r)
        else:
            physio_out.append(p_df)
            physio_out_filt.append(p_df_filt)

    # concatenate all preproc signals into one dataframe
    physio_out = pd.concat(physio_out, axis=1)
    physio_out_filt = pd.concat(physio_out_filt, axis=1)
    fp_out = fp['out'].format(subj, scan)
    # write out columns as separate text files
    fp_out = f"{output_dict['physio']['proc']}/{fp_out}"
    for col in physio_out.columns:
        np.savetxt(f'{fp_out}_{col}.txt', physio_out[col].values)
        np.savetxt(f'{fp_out}_{col}_filt.txt', physio_out_filt[col].values)


def preprocess(dataset, n_cores, anat_skip, func_skip, physio_skip,
                no_eeglab):
    # master function for preprocessing datasets
    print(f'preprocessing {dataset}')
    # load analysis_params.json to get dataset tr
    params_json = json.load(open(params_fp, 'rb'))
    # Set dataset preprocessing parameters
    if dataset in ['chang', 'chang_bh', 'chang_cue']:
        params_dataset = params_json[dataset]
        params = {
            'p_type': 'multiecho', # minimal, multiecho or full preprocessing pipeline
            'mask': params_dataset['mask'], # path to binary brain mask
            'robustfov': True, # whether to crop anatomical image
            'bet_frac': 0.5, # bet fractional intensity threshold (0 - 1): higher -> more aggresive
            'slicetime': 'data/dataset_chang/slicetiming_chang.txt', # path to slice timing file or boolean (read from header)
            'smooth': True, # whether to smooth (5mm fwhm) functional scan
            'trim': 7, # number of volumes to trim from begin of functional scan (if negative, trim from end)
            'n_cores': n_cores, 
            'eeg': True, # whether eeg is collected in this dataset
            'tr': params_dataset['tr'], # functional TR
            'resample_physio': 100, # resample frequency for physio,
            'trim_physio': 14.7, # time (in secs) to trim off front of physio signals,
            'resample_to_func': True, # whether to resample physio to functional scan length,
            'echo_times': params_dataset['echotimes'] # echo times for multiecho scan
        }

    elif dataset == 'hcp':
        params_dataset = params_json[dataset]
        params = {
            'p_type': 'minimal',
            'mask': params_dataset['mask'],
            'robustfov': False,
            'bet_frac': None,
            'slicetime': None,
            'smooth': True,
            'trim': None,
            'n_cores': n_cores,
            'eeg': False,
            'tr': params_dataset['tr'],
            'resample_physio': None,
            'trim_physio': None,
            'resample_to_func': True 
        }
    elif dataset == 'natview':
        params_dataset = params_json[dataset]
        params = {
            'p_type': 'full',
            'mask': params_dataset['mask'],
            'robustfov': True,
            'bet_frac': 0.25,
            'slicetime': 'data/dataset_natview/slicetiming_natview.txt',
            'smooth': True,
            'trim': -2,
            'n_cores': n_cores,
            'eeg': True,
            'tr': params_dataset['tr'],
            'resample_physio': None,
            'trim_physio': None,
            'resample_to_func': True 
        }
    elif dataset == 'nki':
        params_dataset = params_json[dataset]
        params = {
            'p_type': 'full',
            'mask': params_dataset['mask'],
            'robustfov': False,
            'bet_frac': 0.25,
            'slicetime': None,
            'smooth': True,
            'trim': 7,
            'n_cores': n_cores,
            'eeg': False,
            'tr': params_dataset['tr'],
            'resample_physio': None,
            'trim_physio': 9.8,
            'resample_to_func': True 
        }
    elif dataset == 'nki_rest':
        params_dataset = params_json[dataset]
        params = {
            'p_type': 'full',
            'mask': params_dataset['mask'],
            'robustfov': False,
            'bet_frac': 0.25,
            'slicetime': None,
            'smooth': True,
            'trim': 7,
            'n_cores': n_cores,
            'eeg': False,
            'tr': params_dataset['tr'],
            'resample_physio': None,
            'trim_physio': 9.8,
            'resample_to_func': True 
        }
    elif dataset == 'spreng':
        params_dataset = params_json[dataset]
        params = {
            'p_type': 'multiecho',
            'mask': params_dataset['mask'],
            'robustfov': False,
            'bet_frac': 0.25,
            'slicetime': f'data/dataset_spreng/slicetiming_spreng.txt',
            'smooth': True,
            'trim': 4,
            'n_cores': n_cores,
            'eeg': False,
            'tr': params_dataset['tr'],
            'resample_physio': None,
            'trim_physio': 12,
            'resample_to_func': True,
            'echo_times': params_dataset['echotimes'] 
        }
    elif dataset == 'toronto':
        params_dataset = params_json[dataset]
        params = {
            'p_type': 'full',
            'mask': params_dataset['mask'],
            'robustfov': True,
            'bet_frac': 0.5,
            'slicetime': None,
            'smooth': True,
            'trim': 26,
            'n_cores': n_cores,
            'eeg': False,
            'tr': params_dataset['tr'],
            'resample_physio': None,
            'trim_physio': 9.88,
            'resample_to_func': True 
        }
    elif (dataset == 'yale'):
        params_dataset = params_json[dataset]
        params = {
            'p_type': 'full',
            'mask': params_dataset['mask'],
            'robustfov': False,
            'bet_frac': 0.25,
            'slicetime': None, 
            'smooth': True,
            'trim': 10,
            'n_cores': n_cores,
            'eeg': False,
            'tr': params_dataset['tr'],
            'resample_physio': None,
            'trim_physio': None,
            'resample_to_func': False
        }

    # Pull subject and scan labels for each subject
    subj, scan = load_subject_list(dataset, params_dataset['subject_list'])
    # get output directories
    output_dict = create_directories(dataset, params['p_type'], params['eeg'])
    # get filepaths to functional and anatomical scans
    params['func'], params['anat'], params['physio'] = get_fp(dataset)
    # apply preprocessing pipeline (possibly in parallel)
    preprocess_map(
        subj, scan, params, output_dict, dataset, anat_skip, 
        func_skip, physio_skip, no_eeglab
    )


def preprocess_map(subj, scan, params, output_dict, dataset, 
                   anat_skip, func_skip, physio_skip, no_eeglab):
    # apply preprocessing pipeline to each subject in parallel
    pool = Pool(processes=params['n_cores'])
    # Full preprocessing pipeline - starting from raw
    if (params['p_type'] == 'full') | (params['p_type'] == 'multiecho'):
        # anatomical pipeline
        # get unique subj ids while preserving order
        subj_unq = list(dict.fromkeys(subj))
        # skip, if specified
        if anat_skip:
            anat_out_dict = get_anat_fp(params['anat'], subj_unq, output_dict)
        else:
            # Apply anatomical pipeline to structural scans (possibly in parallel)
            anat_iter = zip(repeat(params['anat']), subj_unq, repeat(output_dict),
                            repeat(params['robustfov']), repeat(params['bet_frac']))
            anat_out = pool.starmap(anat_proc, anat_iter)
            # convert anat output to dict with subj id as keys
            anat_out_dict = {a[0]: a[1] for a in anat_out}
        # functional pipeline
        if not func_skip:
            if params['p_type'] == 'full':
                func_iter = zip(
                 repeat(params['func']), subj, scan, repeat(anat_out_dict),
                 repeat(output_dict), repeat(params['tr']), 
                 repeat(params['mask']), repeat(params['slicetime']), 
                 repeat(params['trim']),
                )
                pool.starmap(func_full_proc, func_iter)
            elif params['p_type'] == 'multiecho':
                func_iter = zip(
                 repeat(params['func']), repeat(params['echo_times']),
                 subj, scan, repeat(anat_out_dict),
                 repeat(output_dict), repeat(params['tr']), 
                 repeat(params['mask']), repeat(params['slicetime']),
                  repeat(params['trim']),
                )
                pool.starmap(func_me_proc, func_iter)
                # func_me_proc(
                #     params['func'], params['echo_times'],
                #     subj[0], scan[0],
                #     anat_out_dict, output_dict,
                #     params['tr'], params['mask'], params['slicetime'],
                #     params['trim']
                # )
     # Minimal preprocessing pipeline - starting from preprocessed
    elif params['p_type'] == 'minimal':
        if not func_skip:
            func_iter = zip(repeat(params['func']), subj, scan, 
                            repeat(output_dict), repeat(params['tr']), 
                            repeat(params['mask']), repeat(True), 
                            repeat(params['smooth'])
                            )
            pool.starmap(func_minimal_proc, func_iter)

    # Physio preprocessing
    if not physio_skip:
        if params['p_type'] == 'multiecho':
            func_template = params['func']['func']
        else:
            func_template = params['func'] 

        physio_iter = zip(
          repeat(params['physio']), subj, scan, repeat(dataset), 
          repeat(func_template), repeat(output_dict), 
          repeat(params['resample_physio']), repeat(params['trim_physio']), 
          repeat(params['resample_to_func']), repeat(no_eeglab)
        )
        pool.starmap(physio_proc, physio_iter)
        # physio_proc(params['physio'], subj[0], scan[0], dataset,
        #             func_template, output_dict, params['resample_physio'], 
        #             params['trim_physio'], params['resample_to_func'], no_eeglab)

def tedana_denoise(fps_in, echo_times, mask, out_dir, out_prefix, 
                   fittype='curvefit'):
    # function-scoped import so tedana need not be installed unless needed
    from tedana.workflows import tedana_workflow
    # run tedana workflow
    tedana_workflow(data=fps_in, tes=echo_times, mask=mask, fittype=fittype, 
                    prefix=out_prefix, out_dir=out_dir, overwrite=True,
                    tedpca='kundu-stabilize')


if __name__ == '__main__':
    """Run preprocessing"""
    parser = argparse.ArgumentParser(description='Preprocess datasets')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to preprocess - '
                        'to run all datasets use the arg "all"',
                        choices=['all', 'chang', 'chang_bh', 'chang_cue', 
                                 'natview', 'nki', 'nki_rest', 
                                 'hcp', 'spreng', 'toronto', 'yale'
                                 ], 
                        required=True,
                        type=str)
    parser.add_argument('-n', '--n_cores',
                        help='number of cores to use for parallel processing',
                        default = 1,
                        required=False,
                        type=int)
    parser.add_argument('-anat_skip', '--anat_skip',
                        help='whether to run preprocessing '
                        'without having to re-run anatomical pipeline',
                        action='store_true')
    parser.add_argument('-func_skip', '--func_skip',
                        help='whether to run preprocessing without'
                        ' having to re-run functional pipeline',
                        action='store_true')
    parser.add_argument('-physio_skip', '--physio_skip',
                        help='whether to run preprocessing without'
                        ' having to re-run physio pipeline',
                        action='store_true')
    parser.add_argument('-no_eeglab', '--no_eeglab',
                        help='Whether to skip eeglab preprocessing for natview dataset',
                        action='store_true')

    args_dict = vars(parser.parse_args())
    # if dataset == 'all', run through all datasets and preprocess
    if args_dict['dataset'] == 'all':
        for d in datasets:
            preprocess(d, args_dict['n_cores'], args_dict['anat_skip'], 
                       args_dict['func_skip'], args_dict['physio_skip'], 
                       args_dict['no_eeglab'])
    else:
        preprocess(args_dict['dataset'], args_dict['n_cores'], 
                   args_dict['anat_skip'], args_dict['func_skip'], 
                   args_dict['physio_skip'], args_dict['no_eeglab'])



