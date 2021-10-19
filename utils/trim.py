import argparse
import nibabel as nb
import numpy as np
import os


def run_main(nifti, n_trim, output_dir):
    nii = nb.load(nifti)
    data = nii.get_fdata()[:,:,:,n_trim:]
    nb.Nifti1Image(data, nii.affine, nii.header).to_filename(output_dir)


if __name__ == '__main__':
    """Trim first N volumes from nifti """
    parser = argparse.ArgumentParser(description='Normalize and bandpass filter cifti files')
    parser.add_argument('-f', '--nifti',
                        help='<Required> path to nifti file',
                        required=True,
                        type=str)
    parser.add_argument('-n', '--n_trim',
                        help='how many volumes (from the first) to trim off (index starts at 0)',
                        required=True,
                        type=int)
    parser.add_argument('-o', '--output_dir',
                        help='output directory',
                        required=False,
                        default=os.getcwd(),
                        type=str)
    args_dict = vars(parser.parse_args())
    run_main(args_dict['nifti'], args_dict['n_trim'],
             args_dict['output_dir'])
