# include parent directory in interpeter path
import sys
sys.path.append('..')

import argparse
import json
import nibabel as nb
import numpy as np
import os

from itertools import repeat
from multiprocessing import Pool
from preprocess import (
    get_anat_fp, get_fp, func_minimal_proc
)
from tedana import io
from tedana.decay import fit_decay_ts
from tedana.utils import make_adaptive_mask, unmask
from utils.fsl_utils import warp_func
from utils.load_write import (
    get_fp_base, load_subject_list,
    convert_2d, convert_4d
)

# path to the analysis parameter .json file
params_fp='../analysis_params.json'


def compute_me_maps(fp_echos, echo_times, mask, out_file):
    # compute t2* and S0 images from multiecho data
    # load data into 3d array (S x E [x T]) array_like
    # where `S` is samples, `E` is echos, and `T` is
    catd, ref_img = io.load_data(fp_echos, n_echos=3)
    # get adaptive mask
    mask_denoise, masksum_denoise = make_adaptive_mask(
        catd,
        mask=mask,
        getsum=True,
        threshold=1,
    )
    # estimate t2* and S0 images
    t2s_limited_ts, s0_limited_ts, t2s_full_ts, s0_full_ts = fit_decay_ts(
        catd, echo_times, mask_denoise, masksum_denoise, 'loglin'
    )
    t2_out_file = f'{out_file}_t2.nii.gz'
    write_map(t2s_full_ts, mask, t2_out_file)
    s0_out_file = f'{out_file}_s0.nii.gz'
    write_map(s0_full_ts, mask, s0_out_file)
    return t2_out_file, s0_out_file


def create_directories(dataset):
    # directories for preprocessing
    output_dict = {
        'anat': {
            'raw': f'../data/dataset_{dataset}/anat/raw',
            'reorient': f'../data/dataset_{dataset}/anat/proc1_reorient',
            'bet': f'../data/dataset_{dataset}/anat/proc2_bet',
            'fast': f'../data/dataset_{dataset}/anat/proc3_fast',
            'flirt': f'../data/dataset_{dataset}/anat/proc4_flirt',
            'fnirt': f'../data/dataset_{dataset}/anat/proc5_fnirt'
        },
        'func': {
            'mcflirt': f'../data/dataset_{dataset}/func/proc1_mcflirt',
            'tedana': f'../data/dataset_{dataset}/func/proc2_tedana',
            'func2struct': f'../data/dataset_{dataset}/func/proc3_func2struct',
            'multiechomaps': f'../data/dataset_{dataset}/func/proc_me/proc4_memaps',
            'standard': f'../data/dataset_{dataset}/func/proc_me/proc5_standard',
            'smooth': f'../data/dataset_{dataset}/func/proc_me/proc6_mask_smooth',
            'bandpass': f'../data/dataset_{dataset}/func/proc_me/proc7_bandpass'
        }
    }
    return output_dict


def func_me_proc(fp_me, echo_times, subj, scan, anat_out_dict, 
                 output_dict, tr, mask):
    # estimate preprocessed t2* and S0 images
    # get functional file name
    fp_func = fp_me['func'].format(subj, scan)
    fp_func_base = get_fp_base(fp_func)
    # get anatomical filepaths
    anat_out = anat_out_dict[subj]
    # get file paths of echos
    mcflirt_dir = f"{output_dict['func']['mcflirt']}"
    fp_echos = [
        f"{mcflirt_dir}/{fp_me['echo'].format(subj, scan, i+1)}"
        for i in range(len(echo_times))
    ]
    # define path to booelan mask generated by BET in anat pipeline
    tedana_out_dir = f"{output_dict['func']['tedana']}/{fp_func_base}"
    fp_bet_base = f"{os.path.basename(get_fp_base(anat_out['bet']))}_mask"
    fp_bet_mask = f"{tedana_out_dir}/{fp_bet_base}_func.nii.gz"
    # define path to coregistration affine mat
    fp_func2struct = f"{output_dict['func']['func2struct']}/{fp_func_base}.mat"
    # define output file paths for t2* and S0 images
    fp_me_maps = f"{output_dict['func']['multiechomaps']}/{fp_func_base}"
    fp_t2_full, fp_s0_full = compute_me_maps(
        fp_echos, echo_times, fp_bet_mask, fp_me_maps
    )
    fp_t2, fp_s0 = [
        os.path.basename(f) 
        for f in (fp_t2_full, fp_s0_full)
    ]
    # apply minimal functional preprocessing pipeline to each image
    for fp, fp_full in zip((fp_t2, fp_s0), (fp_t2_full, fp_s0_full)):
        # apply warp to get functional to MNI
        fp_warp = f"{output_dict['func']['standard']}/{fp}"
        warp_func(fp_full, fp_func2struct, anat_out['fnirt_coef'], fp_warp, mask)
        # apply func minimal preprocessing pipeline
        func_minimal_proc(fp, subj, scan, output_dict, tr, mask, 
                          resample=False, smooth=True)


def preprocess_multiecho(subj, scan, params, output_dict):
    # create directories for multiecho output
    for out in output_dict['func']:
        os.makedirs(output_dict['func'][out], exist_ok=True)
    # apply preprocessing pipeline to each subject in parallel
    pool = Pool(processes=params['n_cores'])
    # get unique subj ids while preserving order
    subj_unq = list(dict.fromkeys(subj))
    # get file paths to anatomical preprocessing output
    anat_out_dict = get_anat_fp(params['anat'], subj_unq, output_dict)
    func_iter = zip(
     repeat(params['func']), repeat(params['echo_times']),
     subj, scan, repeat(anat_out_dict),
     repeat(output_dict), repeat(params['tr']), 
     repeat(params['mask']),
    )
    pool.starmap(func_me_proc, func_iter)
    # func_me_proc(
    #     params['func'], params['echo_times'],
    #     subj[0], scan[0], anat_out_dict,
    #     output_dict, params['tr'], 
    #     params['mask'],
    # )


def write_map(data, ref_img, out_file):
    """
    write 2d time series to 4d nifti
    """
    vox_type = data.dtype
    if vox_type == np.int64:
        data = np.int32(data)
    elif vox_type == np.float64:
        data = np.float32(data)
    # Make new img and save
    img = io.new_nii_like(ref_img, data)
    img.to_filename(out_file)


if __name__ == '__main__':
    """
    Apply preprocessing pipeline to t2* and S0 images computed from
    multi-echo data from the spreng dataset. We start with the AFNI multiecho
    preprocessing output from the main preprocessing pipeline (which means that
    the spreng dataset must have been previously preprocessed with the main pipeline).
    Tedana is used to compute the T2* and S0* time series, which are then preprocessed
    in subsequent steps similar to that of the main preprocessing pipeline (normalization,
    smoothing and filtering)
    """
    parser = argparse.ArgumentParser(description='Multiecho preprocess datasets')
    parser.add_argument('-d', '--dataset',
                        help='<Required> Dataset to preprocess - '
                        'to run all datasets use the arg "all"',
                        choices=['chang', 'chang_bh', 'chang_cue', 'spreng'], 
                        required=True,
                        type=str)
    parser.add_argument('-n', '--n_cores',
                        help='number of cores to use for parallel processing',
                        default = 1,
                        required=False,
                        type=int)
    # load analysis_params.json
    params_json = json.load(open(params_fp, 'rb'))
    # parse arguments
    args_dict = vars(parser.parse_args())
    dataset = args_dict['dataset']
    # set file path and analysis parameters
    params = {}
    params['func'], params['anat'], _ = get_fp(dataset)
    params['echo_times'] = params_json[dataset]['echotimes']
    params['tr'] = params_json[dataset]['tr']
    params['mask'] = f"../{params_json[dataset]['mask']}"
    params['n_cores'] = args_dict['n_cores']
    # load subject list
    subj, scan= load_subject_list(
        dataset, f"../{params_json[dataset]['subject_list']}"
    )
    # map directories for input/output files
    output_dict = create_directories(dataset)
    # preprocess multiecho data
    preprocess_multiecho(subj, scan, params, output_dict)
