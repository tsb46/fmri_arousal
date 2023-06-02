import argparse
import json
import mne
import nibabel as nb
import neurokit2 as nk
import numpy as np
import os
import pandas as pd
import shutil

from nipype.interfaces import fsl
from nipype.interfaces.utility import Function
from itertools import repeat
from multiprocessing import Pool
from scipy.io import loadmat
from scipy.stats import zscore
from utils.load_write import convert_2d, convert_4d
from utils.physio_utils import (
    extract_ecg_signals, extract_eeg_signals,
    extract_gsr_signals, extract_ppg_signals, 
    extract_resp_signals
)
from utils.signal_utils import butterworth_filter, clip_spikes

# Ensure output is .nii.gz
fsl.FSLCommand.set_default_output_type('NIFTI_GZ')

# Dilated mask that includes sinuses slightly outside gray matter tissue
mask="masks/MNI152_T1_3mm_brain_mask_dilated.nii.gz"

# path to the analysis parameter .json file - contains dataset tr
params_fp='analysis_params.json'


# datasets
datasets = ['chang', 'chang_bh', 'chang_cue', 
            'nki', 'hcp', 'spreng', 'yale', 
            'natview']

# eeg channel selections for chang and natview datasets
chang_eeg_chan = ['P3', 'P4', 'Pz', 'O1', 'O2', 'Oz']
natview_eeg_chan = []


def anat_proc(fp, subj, output_dict):
    fp = fp.format(subj)
    fp_base = get_fp_base(fp)
    # reorient structural to standard (just in case)
    fp_in = f"{output_dict['anat']['raw']}/{fp}"
    fp_reorient = f"{output_dict['anat']['reorient']}/{fp}"
    reorient(fp_in, fp_reorient)
    # brain extraction
    fp_bet = f"{output_dict['anat']['bet']}/{fp}"
    bet(fp_reorient, fp_bet)
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


def bet(fp, fp_out):
    # BET - Skullstrip anatomical Image
    bet_anat = fsl.BET(frac=0.25, robust=True, mask=True)
    bet_anat.inputs.in_file = fp
    bet_anat.inputs.out_file = fp_out
    bet_anat_res = bet_anat.run()


def concat_transform(fp_func2struct, fp_flirt, fp_out):
    # Concatenate affine transform matrices (func2struct & struct2MNI)
    convertxfm = fsl.ConvertXFM(concat_xfm=True)
    convertxfm.inputs.in_file = fp_func2struct
    convertxfm.inputs.in_file2=fp_flirt
    convertxfm.inputs.out_file=fp_out
    convertxfm_res = convertxfm.run()


def coregister(fp_mean, fp_reorient, fp_bet, fp_wm, fp_out):
    fp_out_base = get_fp_base(fp_out)
    # Coregister functional with T1w
    epireg = fsl.EpiReg()
    epireg.inputs.epi = fp_mean
    epireg.inputs.t1_head=fp_reorient
    epireg.inputs.t1_brain=fp_bet
    # epireg expects the wmseg output as a suffix to the epi image (weird)
    # rename for now
    wmseg = f'{fp_out_base}_fast_wmseg.nii.gz'
    shutil.copyfile(fp_wm, wmseg)
    epireg.inputs.wmseg = wmseg
    epireg.inputs.out_base = fp_out_base
    epireg_res = epireg.run()
    fp_func2struct = epireg_res.outputs.epi2str_mat
    return fp_func2struct


def create_directories(dataset, p_type, eeg, slicetime=None, trim=None):
    # create directories for preprocessing
    if p_type == 'full':
        output_dict = {
            'anat': {
                'raw': f'data/dataset_{dataset}/anat/raw',
                'reorient': f'data/dataset_{dataset}/anat/proc1_reorient',
                'bet': f'data/dataset_{dataset}/anat/proc2_bet',
                'fast': f'data/dataset_{dataset}/anat/proc3_fast',
                'flirt': f'data/dataset_{dataset}/anat/proc4_flirt',
                'fnirt': f'data/dataset_{dataset}/anat/proc5_fnirt'
            },
            'func': {
                'raw': f'data/dataset_{dataset}/func/raw',
                'trim': f'data/dataset_{dataset}/func/procA_trim',
                'slicetime': f'data/dataset_{dataset}/func/procB_slicetime',
                'mcflirt': f'data/dataset_{dataset}/func/proc1_mcflirt',
                'func2struct': f'data/dataset_{dataset}/func/proc2_func2struct',
                'standard': f'data/dataset_{dataset}/func/proc3_standard',
                'smooth': f'data/dataset_{dataset}/func/proc4_mask_smooth',
                'bandpass': f'data/dataset_{dataset}/func/proc5_bandpass'
             }
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


def fast(fp, fp_out):
    # FAST - Image Segmentation
    fast = fsl.FAST()
    fast.inputs.in_files = fp
    fast.inputs.out_basename = fp_out
    # Nipype FAST issue with writing out tissue_class_map - Ignore
    # https://github.com/nipy/nipype/issues/3311
    try: 
        fast_res = fast.run()
    except FileNotFoundError:
        fast_out = fp_out
    fp_out_base = get_fp_base(fp_out)
    fp_out_wm = f'{fp_out_base}_pve_2'
    return fp_out_wm


def flirt(fp, fp_out, fp_out_mat):
    # FLIRT affine registration to MNI template
    flirt = fsl.FLIRT()
    flirt.inputs.in_file = fp
    flirt.inputs.reference = f'{os.environ["FSLDIR"]}/data/standard/MNI152_T1_2mm_brain.nii.gz'
    flirt.inputs.out_file = fp_out
    flirt.out_matrix_file = fp_out_mat
    flirt_res = flirt.run()
    # Flirt saves output matrix in base directory (seems to be an issue related to the FAST issue above), 
    # move to results directory
    os.rename(flirt_res.outputs.out_matrix_file, fp_out_mat)


def fnirt(fp, fp_affine, fp_out, fp_out_coef):
    # FNIRT non-linear registration
    fnirt = fsl.FNIRT()
    fnirt.inputs.in_file = fp
    fnirt.inputs.ref_file = f'{os.environ["FSLDIR"]}/data/standard/MNI152_T1_2mm.nii.gz'
    fnirt.inputs.affine_file = fp_affine
    fnirt.inputs.config_file='T1_2_MNI152_2mm'
    fnirt.inputs.warped_file=fp_out
    fnirt.inputs.fieldcoeff_file = fp_out_coef
    fp_out_base = get_fp_base(fp_out)
    fnirt.inputs.log_file = f'{fp_out_base}_log.txt'
    fnirt_res = fnirt.run()


def func_full_proc(fp, subj, scan, anat_out_dict, output_dict, tr,
                   slicetime=None, trim=None):
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
    if slicetime:
        fp_slicetime = f"{output_dict['func']['slicetime']}/{fp}"
        slicetime(fp_in, fp_slicetime, slicetime)
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
    warp_func(fp_mcflirt, fp_func2struct, anat_out['fnirt_coef'], fp_warp)
    # apply func minimal preprocessing pipeline
    func_n = func_minimal_proc(fp, subj, scan, output_dict, tr, resample=False)
    return func_n


def func_minimal_proc(fp, subj, scan, output_dict, tr, resample=True):
    # minimal preprocessing pipeline starting with MNI-registred and 
    # preprocessed functional data
    fp = fp.format(subj, scan)
    fp_base = get_fp_base(fp)
    if resample:
        # resample functional to 3mm MNI space
        fp_in = f"{output_dict['func']['raw']}/{fp}"
        fp_resample = f"{output_dict['func']['resample']}/{fp}"
        resample_func(fp_in, fp_resample)
        fp_in = fp_resample
    else:
        fp_in = f"{output_dict['func']['standard']}/{fp}"
    # spatial smoothing (and mask)
    # Get mask
    fp_smooth = f"{output_dict['func']['smooth']}/{fp}"
    smooth(fp_in, fp_smooth)
    # bandpass filter
    # Load mask
    mask_bin = nb.load(mask).get_fdata() > 0
    fp_bandpass = f"{output_dict['func']['bandpass']}/{fp}"
    bandpass_filter(fp_smooth, fp_bandpass, tr, mask_bin)
    # get number of time points from final functional
    func_n = nb.load(fp_bandpass).shape[-1]


def get_fp(dataset):
    # set filepath templates for raw data
    if dataset == 'chang':
        func = 'sub_00{0}-mr_{1}-ecr_echo1_w_dspk_blur3mm.nii.gz'
        anat = None
        physio = {
            'physio': 'sub_00{0}-mr_{1}-ecr_echo1_physOUT.mat',
            'eeg': 'sub_00{0}-mr_{1}_eeg_pp.mat',
            'out': 'sub_00{0}-mr_{1}-ecr_physio'
        }
    elif dataset == 'chang_bh':
        func = 'sub_00{0}-mr_{1}-adb_echo1_w_dspk_blur3mm.nii.gz'
        anat = None
        physio = {
            'physio': 'sub_00{0}-mr_{1}-adb_echo1_physOUT.mat',
            'eeg': 'sub_00{0}-mr_{1}-adb_echo1_EEG_pp.mat',
            'out': 'sub_00{0}-mr_{1}-adb_physio'
        }
    elif dataset == 'chang_cue':
        func = 'sub_00{0}-mr_{1}-ectp_echo1_w_dspk_dtr_blur3mm.nii.gz'
        anat = None
        physio = {
            'physio': 'sub_00{0}-mr_{1}-ectp_echo1_physOUT.mat',
            'eeg': 'sub_00{0}-mr_{1}-ectp_echo1_EEG_pp.mat',
            'out': 'sub_00{0}-mr_{1}-ectp_physio'
        }
    elif dataset == 'hcp':
        func = '{0}_{1}_hp2000_clean.nii.gz'
        anat = None
        physio = {
          'physio': '{0}_{1}_physio.txt',
          'out': '{0}_{1}_physio'
        }
    elif dataset == 'natview':
        func='sub-{0}_ses-0{1}_func_mc.nii.gz'
        anat = None
        physio = {
            'eeg': 'sub-{0}_ses-0{1}_task-rest_eeg.set',
            'eye': 'sub-{0}_ses-0{1}_task-rest_recording-eyetracking_physio.tsv.gz',
            'resp': 'sub-{0}_ses-0{1}_task-rest_recording-respiratory_physio.tsv.gz',
            'out': 'sub-{0}_ses-0{1}_task-rest_physio'
        }
    elif dataset == 'nki': 
        func = '{0}_task_breathhold.nii.gz'
        anat = '{0}_T1w.nii.gz'
        physio = {
            'physio': '{0}_task_breathhold_physio.tsv.gz',
            'out': '{0}_task_breathhold_physio'
        }
    elif dataset == 'spreng':
        func = '{0}_task-rest_{1}_echo-123_bold_medn_afw.nii.gz'
        anat = None
        physio = {
            'physio': '{0}_task-rest_{1}_physio.tsv.gz',
            'out': '{0}_task-rest_{1}_physio'
        }
    elif dataset == 'yale':
        func = '{0}_task-rest_run-0{1}_bold.nii.gz'
        anat = '{0}_T1w.nii.gz'
        physio = {
            'physio': '{0}_task-rest_run-0{1}_et.tsv',
            'out': '{0}_task-rest_run-0{1}_physio'
        }
    return func, anat, physio


def get_fp_base(fp):
    # get nifti file path without extenstion
    fp_split = os.path.splitext(fp)
    if fp_split[1] == '.gz':
        fp_base = os.path.splitext(fp_split[0])[0]
    else:
        fp_base = fp_split[0]
    return fp_base


def load_subject_list(dataset, subject_list_fp):
    subj_df = pd.read_csv(subject_list_fp)
    # load subject list for a dataset
    if dataset == 'chang':
        subj = subj_df.subject.tolist()
        scan = [f'000{s}' if s <10 else f'00{s}' for s in subj_df.scan] 
    elif dataset == 'chang_bh':
        subj = subj_df.subject.tolist()
        scan = [f'000{s}' if s <10 else f'00{s}' for s in subj_df.scan]
    elif dataset == 'chang_cue':
        subj = subj_df.subject.tolist()
        scan = [f'000{s}' if s <10 else f'00{s}' for s in subj_df.scan]
    elif dataset == 'hcp':
        subj = subj_df.subject.tolist()
        scan = subj_df.lr.tolist()
    elif dataset == 'natview':
        subj = [f'0{s}' if s <10 else f'{s}' for s in subj_df.subject]
        scan = subj_df.scan.tolist()
    elif dataset == 'nki': 
        subj = subj_df.subject.tolist()
        scan = [None] * len(subj)
    elif dataset == 'spreng':
        subj = subj_df.subject.tolist()
        scan = subj_df.scan.tolist()
    elif dataset == 'yale':
        subj = subj_df.subject.tolist()
        scan = subj_df.scan.tolist()
    return subj, scan


def load_eeg(fp, dataset, resample, trim):
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

    return eeg_mne


def load_physio(fp, subj, scan, dataset, output_dict, resample, trim):
    # load physiological signals for each dataset and return dictionary with 
    # physio labels as keys.
    # each dataset has unique file formats that must be accounted for    
    if dataset in ['chang', 'chang_bh', 'chang_cue']:
        fp_p = fp['physio'].format(subj, scan)
        fp_eeg = fp['eeg'].format(subj, scan)
        fp_p_in = f"{output_dict['physio']['raw']}/{fp_p}"
        fp_eeg_in = f"{output_dict['eeg']['raw']}/{fp_eeg}"
        # Trim off first 14.7s to align w/ functional (first 7 TRs were trimmed from functional)
        physio_raw = loadmat(fp_p_in, squeeze_me=True)
        sf = 1/physio_raw['OUT_p']['dt_phys'].item()
        # Pull physio data into dict
        ppg = physio_raw['OUT_p']['card_dat'].item()
        resp = physio_raw['OUT_p']['resp'].item()['wave'].item()
        eeg = load_eeg(fp_eeg_in, dataset, resample, trim)
        # write out eeg object
        fp_eeg_out = get_fp_base(fp_eeg)
        eeg.save(f"{output_dict['eeg']['proc']}/{fp_eeg_out}.raw.fif", overwrite=True)
        physio = {'ppg': ppg, 'resp': resp, 'eeg': eeg}

    elif (dataset == 'nki') | (dataset == 'spreng'):
        fp_p = fp['physio'].format(subj, scan)
        fp_p_in = f"{output_dict['physio']['raw']}/{fp_p}"
        physio_df = pd.read_csv(fp_p_in, compression='gzip', sep='\t', header=None)
        physio_df = physio_df.dropna(axis=1, how='all')
        # metadata jsons for each dataset need separate handling
        if dataset == 'nki':
            fp_json = get_fp_base(fp_p_in)
        elif dataset == 'spreng':
            fp_j = get_fp_base(fp['physio'].format(subj, '').replace('__','_'))
            fp_json = f"{output_dict['physio']['raw']}/{fp_j}"
        # load json and set columns
        physio_json = json.load(open(f'{fp_json}.json'))
        physio_df.columns = physio_json['Columns']
        physio = {
            'ppg': physio_df['cardiac'].values, 
            'resp': physio_df['respiratory'].values
        }
        if dataset == 'nki':
            physio['gsr'] = physio_df['gsr'].values
        sf = physio_json["SamplingFrequency"]
    elif dataset == 'hcp':
        sf = 400 # hcp physio sampling frequency
        fp_p = fp['physio'].format(subj, scan)
        fp_p_in = f"{output_dict['physio']['raw']}/{fp_p}"
        physio_raw = np.loadtxt(fp_p_in)
        physio = {
            'resp': physio_raw[:,1], 
            'ppg': physio_raw[:,2]
        }

    elif dataset == 'yale':
        sf = 1 # already resampled to functional scan TR
        fp_p = fp['physio'].format(subj, scan)
        fp_p_in = f"{output_dict['physio']['raw']}/{fp_p}"
        physio_df = pd.read_csv(fp_p_in, sep='\t', header=None)
        physio = {
            'pupil': physio_df.iloc[:,0].values
        }

    # loop through physio signals and trim/resample (if specified)
    for p in physio.keys():
        if p != 'eeg':
            # some dataset physio signals need to be trimmed to align w/ functional
            if trim is not None:
                trim_n = int(sf*trim)
                physio[p] = physio[p][trim_n:]
            # to ease computational burden, some high-frequency 
            # physio signals are downsampled before pre-processing
            if resample is not None:
                physio[p] = nk.signal_resample(
                               physio[p], 
                               sampling_rate=sf, 
                               desired_sampling_rate=resample, 
                               method='FFT'
                            )
    if resample is not None:
        # set resample freq as new freq (sf)
        sf = resample
    return physio, sf


def mcflirt(fp, fp_out):
    # McFLIRT Motion Correction
    fp_out_base = get_fp_base(fp_out)
    func_file_meanvol = f'{fp_out_base}_mean.nii.gz'
    mcflirt = fsl.MCFLIRT(mean_vol=True, save_plots=True)
    mcflirt.inputs.in_file = fp
    mcflirt.inputs.out_file = fp_out
    mcflirt_res = mcflirt.run()
    # weird renaming of mean vol, rename
    os.rename(mcflirt_res.outputs.mean_img, func_file_meanvol)
    return func_file_meanvol


def physio_proc(fp, subj, scan, dataset, physio_labels, 
                fp_func, output_dict, resample, trim, 
                resample_to_func):
    # preprocess physio signals - save out unfiltered and 
    # bandpass filtered preprocessed signals
    # load physio signals
    physio, sf = load_physio(fp, subj, scan, dataset,
                             output_dict, resample, trim)
    # get n of time points of functional scan for aligning w/ physio
    # need to load in subj functional nifti header
    fp_func_in = f"{output_dict['func']['bandpass']}/{fp_func.format(subj, scan)}"
    func_n = nb.load(fp_func_in).shape[-1]
    # loop through physio signals and extract, clip, filter and resample
    physio_out = []
    physio_out_filt = []
    for p in physio_labels:
        p_df = extract_physio(physio[p], p, sf)
        # clip potential spikes - abs(z) >= 5
        p_df = p_df.apply(clip_spikes, axis=0)
        # Bandpass filter signals (0.01-0.1Hz)
        p_df_filt = p_df.apply(
           butterworth_filter, lowcut=0.01, highcut=0.1, fs=sf, 
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


def preprocess(dataset, n_cores):
    # master function for preprocessing datasets
    print(f'preprocessing {dataset}')
    # load analysis_params.json to get dataset tr
    params_json = json.load(open(params_fp, 'rb'))
    # Set dataset preprocessing parameters
    if dataset in ['chang', 'chang_bh', 'chang_cue']:
        params_dataset = params_json[dataset]
        params = {
            'p_type': 'minimal', # minimal or full preprocessing pipeline
            'slicetime': None, # filepath to slice order file
            'trim': None, # number of volumes to trim from begin of functional scan
            'n_cores': n_cores, 
            'tr': params_dataset['tr'], # functional TR
            'signals': ['resp', 'ppg', 'eeg'], # physio signals in dataset
            'resample_physio': 100, # resample frequency for physio,
            'trim_physio': 14.7, # time (in secs) to trim off front of physio signals,
            'resample_to_func': True # whether to resample physio to functional scan length
        }

    elif dataset == 'hcp':
        params_dataset = params_json[dataset]
        params = {
            'p_type': 'minimal',
            'slicetime': None, 
            'trim': None,
            'n_cores': n_cores,
            'tr': params_dataset['tr'],
            'signals': ['resp', 'ppg'],
            'resample_physio': None,
            'trim_physio': None,
            'resample_to_func': True 
        }
    elif dataset == 'natview':
        params_dataset = params_json[dataset]
        params = {
            'p_type': 'minimal',
            'slicetime': None, 
            'trim': None,
            'n_cores': n_cores,
            'tr': params_dataset['tr'],
            'signals': ['resp', 'pupil', 'eeg'],
            'resample_physio': None,
            'trim_physio': None,
            'resample_to_func': True 
        }
    if dataset == 'nki':
        params_dataset = params_json[dataset]
        params = {
            'p_type': 'full',
            'slicetime': None, 
            'trim': None,
            'n_cores': n_cores,
            'tr': params_dataset['tr'],
            'signals': ['resp', 'ppg', 'gsr'],
            'resample_physio': None,
            'trim_physio': None,
            'resample_to_func': True 
        }
    if dataset == 'spreng':
        params_dataset = params_json[dataset]
        params = {
            'p_type': 'minimal',
            'slicetime': None, 
            'trim': None,
            'n_cores': n_cores,
            'tr': params_dataset['tr'],
            'signals': ['resp', 'ppg'],
            'resample_physio': None,
            'trim_physio': 12,
            'resample_to_func': True
        }
    if (dataset == 'yale') | (dataset == 'all'):
        params_dataset = params_json[dataset]
        params = {
            'p_type': 'full',
            'slicetime': None, 
            'trim': 10,
            'n_cores': n_cores,
            'tr': params_dataset['tr'],
            'signals': ['pupil'],
            'resample_physio': None,
            'trim_physio': None,
            'resample_to_func': False
        }

    # Pull subject and scan labels for each subject
    subj, scan = load_subject_list(dataset, params_dataset['subject_list'])
    # get output directories
    if 'eeg' in params['signals']:
        eeg = True
    else:
        eeg = False
    output_dict = create_directories(dataset, params['p_type'], eeg)
    # get filepaths to functional and anatomical scans
    params['func'], params['anat'], params['physio'] = get_fp(dataset)
    # apply preprocessing pipeline (possibly in parallel)
    preprocess_map(
        subj, scan, params, output_dict, dataset
    )


def preprocess_map(subj, scan, params, output_dict, dataset):
    # apply preprocessing pipeline to each subject in parallel
    pool = Pool(processes=params['n_cores'])
    # Full preprocessing pipeline - starting from raw
    if params['p_type'] == 'full':
        # anatomical pipeline
        # get unique subj ids while preserving order
        subj_unq = list(dict.fromkeys(subj))
        # Apply anatomical pipeline to structural scans (possibly in parallel)
        anat_iter = zip(repeat(params['anat']), subj_unq, repeat(output_dict))
        anat_out = pool.starmap(anat_proc, anat_iter)
        # convert anat output to dict with subj id as keys
        anat_out_dict = {a[0]: a[1] for a in anat_out}
        # functional pipeline
        func_iter = zip(
            repeat(params['func']), subj, scan, repeat(anat_out_dict), 
            repeat(output_dict), repeat(params['tr']), 
            repeat(params['slicetime']), repeat(params['trim'])
        )
        func_n = pool.starmap(func_full_proc, func_iter)

    # Minimal preprocessing pipeline - starting from preprocessed
    elif params['p_type'] == 'minimal':
        func_iter = zip(repeat(params['func']), subj, scan, repeat(output_dict), 
                        repeat(params['tr']))
        pool.starmap(func_minimal_proc, func_iter)

    # Physio preprocessing
    physio_iter = zip(repeat(params['physio']), subj, scan, 
                      repeat(dataset), repeat(params['signals']),
                      repeat(params['func']), repeat(output_dict), 
                      repeat(params['resample_physio']),
                      repeat(params['trim_physio']),
                      repeat(params['resample_to_func'])
                      )


def resample_func(fp, fp_out):
    # FLIRT resample functional scan to 3mm MNI
    flirt = fsl.FLIRT()
    flirt.inputs.in_file = fp
    flirt.inputs.reference = mask
    flirt.inputs.out_file = fp_out
    flirt.inputs.apply_xfm = True
    flirt.inputs.uses_qform = True
    flirt_res = flirt.run()
    # Flirt saves output matrix in base directory (seems to be an issue related to the FAST issue above), 
    # move to results directory
    os.remove(flirt_res.outputs.out_matrix_file)


def reorient(fp, fp_out):
    # Reorient 2 standard
    reorient = fsl.utils.Reorient2Std()
    reorient.inputs.in_file = fp
    reorient.inputs.out_file = fp_out
    reorient_res = reorient.run()


def slicetime(fp, fp_out, st_fp):
    # # Slice time correction
    slicetimer = fsl.SliceTimer(custom_timings=st_fp, 
                                time_repetition=tr) 
    slicetimer.inputs.in_file = fp
    slicetimer.inputs.out_file = fp_out
    slicetimer_res = slicetimer.run()


def smooth(fp, fp_out, fwhm=5.0):
    # 3mm FWHM isotropic smoothing
    smooth = fsl.Smooth(fwhm=fwhm)
    smooth.inputs.in_file = fp
    smooth.inputs.smoothed_file=fp_out
    smooth_res = smooth.run()


def trim_vol(fp, fp_out, n_trim):
    # Trim first N volumes
    trim = fsl.ExtractROI(t_min=n_trim, t_size=-1)
    trim.inputs.in_file = fp
    trim.inputs.roi_file = fp_out
    trim_res = trim.run()


def wm_thres(fp, fp_out):
    # Threshold white matter partial volume
    wm_thres = fsl.Threshold(thresh=0.5, args='-bin')
    wm_thres.inputs.in_file = fp
    wm_thres.inputs.out_file = fp_out
    wm_thres_res = wm_thres.run()


def warp_func(fp, fp_affine, fp_coef, fp_out):
    # Warp functional to MNI space
    applywarp = fsl.ApplyWarp()
    applywarp.inputs.in_file = fp
    applywarp.inputs.ref_file = mask
    applywarp.inputs.premat=fp_affine
    applywarp.inputs.field_file=fp_coef
    applywarp.inputs.out_file=fp_out
    applywarp_res = applywarp.run()


if __name__ == '__main__':
    """Run main analysis"""
    parser = argparse.ArgumentParser(description='Preprocess datasets')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to preprocess - '
                        'to run all datasets use the arg "all"',
                        choices=['all', 'chang', 'chang_bh', 'chang_cue', 
                                 'nki', 'hcp', 'spreng', 'yale', 
                                 'natview'], 
                        required=True,
                        type=str)
    parser.add_argument('-n', '--n_cores',
                        help='number of cores to use for parallel processing',
                        default = 1,
                        required=False,
                        type=int)
    args_dict = vars(parser.parse_args())
    if args_dict['dataset'] == 'all':
        for d in datasets:
            preprocess(d, args_dict['n_cores'])
    else:
        preprocess(args_dict['dataset'], args_dict['n_cores'])



