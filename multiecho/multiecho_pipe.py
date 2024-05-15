# include parent directory in interpeter path
import sys
sys.path.append('..')

import argparse
import json
import nibabel as nb
import numpy as np
import pandas as pd
import os

from itertools import repeat
from multiprocessing import Pool
from preprocess import (
    get_anat_fp, get_fp, func_minimal_proc
)
from tedana.workflows import ica_reclassify_workflow
from utils.fsl_utils import warp_func
from utils.load_write import (
    get_fp_base, load_subject_list
)

# path to the analysis parameter .json file
params_fp='../analysis_params.json'


def compute_me_maps(subj_func, tedana_dir, metric, out_file):
    # compute high kappa (t2*) or high rho (S0) images from multiecho data
    # define file paths to tedana output
    registry_fp = f'{tedana_dir}/{subj_func}_desc-tedana_registry.json'
    metrics_meta = json.load(
        open(f'{tedana_dir}/{subj_func}_desc-ICACrossComponent_metrics.json', 'rb')
    )
    metrics_df = pd.read_csv(
        f'{tedana_dir}/{subj_func}_desc-tedana_metrics.tsv',
        delimiter='\t'
    )
    # get elbows from kappa and rho metrics
    kappa_elbow = metrics_meta['kappa_elbow_kundu']
    rho_elbow = metrics_meta['rho_elbow_kundu']
    # estimate high kappa time courses
    accept, reject = select_components(
        metrics_df, metric, kappa_elbow, rho_elbow
    )
    ica_reclassify_workflow(
        registry_fp, accept, reject, out_file, 
        no_reports=True, overwrite=True
    )


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


def func_me_proc(fp_me, subj, scan, anat_out_dict, 
                 output_dict, tr, mask):
    # estimate preprocessed t2* and S0 images
    # get functional file name
    fp_func = fp_me['func'].format(subj, scan)
    fp_func_base = get_fp_base(fp_func)
    # get anatomical filepaths
    anat_out = anat_out_dict[subj]
    # define path to coregistration affine mat
    fp_func2struct = f"{output_dict['func']['func2struct']}/{fp_func_base}.mat"
    # define tedana dir
    fp_tedana = f"{output_dict['func']['tedana']}/{fp_func_base}" 
    # preprocess high kappa and high rho time courses
    for m in ['kappa', 'rho']:
        fp_me_dir = f"{output_dict['func']['multiechomaps']}/{fp_func_base}_{m}" 
        compute_me_maps(fp_func_base, fp_tedana, m, fp_me_dir)
        fp_me_out = f'{fp_me_dir}/desc-optcomAccepted_bold.nii.gz'
        # apply warp to get functional to MNI
        fp = f'{fp_func_base}_{m}.nii.gz'
        fp_warp = f"{output_dict['func']['standard']}/{fp}"
        warp_func(fp_me_out, fp_func2struct, anat_out['fnirt_coef'], fp_warp, mask)
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
     repeat(params['func']),subj, scan, 
     repeat(anat_out_dict), repeat(output_dict), 
     repeat(params['tr']), repeat(params['mask']),
    )
    # pool.starmap(func_me_proc, func_iter)
    func_me_proc(
        params['func'],
        subj[0], scan[0], anat_out_dict,
        output_dict, params['tr'], 
        params['mask'],
    )

def select_components(metrics_df, metric, k_elbow, r_elbow):
    # select components based on kappa and rho elbow
    if metric == 'kappa':
        k_mask = metrics_df.kappa > k_elbow
        r_mask = metrics_df.rho < r_elbow
    elif metric == 'rho':
        k_mask = metrics_df.kappa < k_elbow
        r_mask = metrics_df.rho > r_elbow

    mask = k_mask & r_mask
    # mask metrics df based on criteria
    accept_df = metrics_df.loc[mask]
    accept_ic = accept_df.Component.str.replace('ICA_', '').astype(int)
    reject_df = metrics_df.loc[~mask]
    reject_ic = reject_df.Component.str.replace('ICA_', '').astype(int)
    return accept_ic.values.tolist(), reject_ic.values.tolist()


if __name__ == '__main__':
    """
    Apply preprocessing pipeline to high kappa (T2) and high rho (S0) volumes computed from
    multi-echo data from multiecho datasets. We start with the tedana output and reclassify
    the components according to kappa and rho criteria, so the main preprocessing script
    must have already been run. The high kappa and high rho volumes are then preprocessed
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
