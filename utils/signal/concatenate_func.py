import argparse
import neurokit2 as nk
import nibabel as nb
import numpy as np
import os

from scipy.stats import zscore



def run_main(niftis, norm, output_dir):
    nii_concat = []
    for nifti in niftis:
        nii = nb.load(nifti)
        nifti_data = nii.get_fdata()
        # zscore along time axis
        if norm:
            nifti_data = zscore(nifti_data, axis=3)
        nii_concat.append(nifti_data)
    # concatenate along time axis
    nii_concat = np.concatenate(nii_concat, axis=3)
    # write out concatenated functional
    nb.Nifti1Image(nii_concat, nii.affine, nii.header).to_filename(output_dir)


if __name__ == '__main__':
    """Trim first N volumes from nifti """
    parser = argparse.ArgumentParser(description='Concatenate')
    parser.add_argument('-f', '--niftis',
                        help='<Required> path to nifti files (appended as a list)',
                        required=True,
                        type=str,
                        action='append')
    parser.add_argument('-n', '--norm',
                        help='<Required> whether to normalize (zscore) functional time series before concatenation',
                        required=False,
                        default=1,
                        type=int)
    parser.add_argument('-o', '--output_dir',
                        help='output directory',
                        required=False,
                        default=os.getcwd(),
                        type=str)
    args_dict = vars(parser.parse_args())
    run_main(args_dict['niftis'], args_dict['norm'], args_dict['output_dir'])
