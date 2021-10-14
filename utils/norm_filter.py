import argparse
import nibabel as nb
import numpy as np
import os

from butterworth_filters import butterworth_filter
from scipy.stats import zscore


def convert_2d(mask, nifti_data):
    nonzero_indx = np.nonzero(mask)
    nifti_2d = nifti_data[nonzero_indx]
    return nifti_2d.T


def convert_4d(mask, nifti_data):
    nifti_4d = np.zeros(mask.shape + (nifti_data.shape[0],), 
                        dtype=nifti_data.dtype)
    nifti_4d[mask, :] = nifti_data.T
    return nifti_4d


def output_nifti(output_file, nifti_file, nifti_4d):
    # for some reason, nibabel is telling me to change the 'size of hdr' field to 540, seems to be fine
    nifti_file.header['sizeof_hdr'] = 540
    nifti_out = nb.Nifti2Image(nifti_4d, nifti_file.affine, nifti_file.header)
    nb.save(nifti_out, output_file)


def run_main(file, nifti, mask, cut_high, cut_low, tr, output_file):
    # Get sampling rate
    fs = 1 / tr
    if nifti:
        # Load mask and functional scan
        nifti = nb.load(file)
        data = nifti.get_fdata()
        nifti.uncache()
        mask = nb.load(mask).get_fdata() > 0
        data = convert_2d(mask, data)
    else:
        data = np.loadtxt(file)
    # Demean and norm data
    data = zscore(data)
    # Lowpass filter with Butterworth
    if cut_low is None:
        data = butterworth_filter(data, None, cut_high, fs, 'lowpass')
    # Bandpass filter with Butterworth
    else:
        data = butterworth_filter(data, cut_low, cut_high, fs, 'bandpass')    

    # Write output
    if nifti:
        nifti_data_4d = convert_4d(mask, data)
        output_nifti(output_file, nifti, nifti_data_4d)
    else:
        np.savetxt(output_file, data)



if __name__ == '__main__':
    """Temporal filter w/ Butterworth Filter """
    parser = argparse.ArgumentParser(description='Normalize and temporal filter nifti files')
    parser.add_argument('-f', '--file_path',
                        help='<Required> path to file',
                        required=True,
                        type=str)
    parser.add_argument('-n', '--nifti',
                        help='whether file is a nifti file, otherwise signal .txt file',
                        default=1,
                        required=False,
                        type=int)
    parser.add_argument('-m', '--mask',
                        help='path to nifti mask file',
                        required=False,
                        type=str)
    parser.add_argument('-ch', '--cut_high',
                        help='<Required> cutoff of the lowpass filter',
                        required=True,
                        type=float)
    parser.add_argument('-cl', '--cut_low',
                        help='cutoff of the highpass filter',
                        required=False,
                        type=float)
    parser.add_argument('-t', '--tr',
                        help='the repetition time of the data',
                        required=False,
                        default=2.1,
                        type=float)
    parser.add_argument('-o', '--output_file',
                        help='output file path',
                        required=False,
                        default=os.getcwd(),
                        type=str)
    args_dict = vars(parser.parse_args())
    run_main(args_dict['file_path'], args_dict['nifti'], args_dict['mask'],
             args_dict['cut_high'], args_dict['cut_low'], 
             args_dict['tr'], args_dict['output_file'])
