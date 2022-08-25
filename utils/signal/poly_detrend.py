import argparse
import neurokit2 as nk
import nibabel as nb
import numpy as np
import os

from utils.load_write import convert_2d, convert_4d



def run_main(nifti, order, mask, output_dir):
    nii = nb.load(nifti)
    nifti_data = nii.get_fdata()
    mask = nb.load(mask).get_fdata() > 0
    nifti_data = convert_2d(mask, nifti_data)
    nifti_data_detrend = np.apply_along_axis(nk.signal_detrend, 0, nifti_data)
    nifti_data_detrend = convert_4d(mask, nifti_data_detrend)
    nb.Nifti1Image(nifti_data_detrend, nii.affine, nii.header).to_filename(output_dir)


if __name__ == '__main__':
    """Trim first N volumes from nifti """
    parser = argparse.ArgumentParser(description='Polynomial detrending of data via least squares')
    parser.add_argument('-f', '--nifti',
                        help='<Required> path to nifti file',
                        required=True,
                        type=str)
    parser.add_argument('-p', '--order',
                        help='order of polynomial fit to the data (and removed); 1 for linear detrend',
                        required=True,
                        type=int)
    parser.add_argument('-m', '--mask',
                        help='path to nifti mask file',
                        required=True,
                        type=str)
    parser.add_argument('-o', '--output_dir',
                        help='output directory',
                        required=False,
                        default=os.getcwd(),
                        type=str)
    args_dict = vars(parser.parse_args())
    run_main(args_dict['nifti'], args_dict['order'], args_dict['mask'], args_dict['output_dir'])
