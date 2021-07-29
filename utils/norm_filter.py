import argparse
import nibabel as nb
import numpy as np
import os

from scipy.signal import butter, sosfiltfilt, sosfreqz
from scipy.stats import zscore

def butter_bandpass(lowcut, highcut, fs, order=5):
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq
    sos = butter(order, [low, high], analog=False, btype='band', output='sos')
    return sos


def butter_bandpass_filter(data, lowcut, highcut, fs, order=5):
    sos = butter_bandpass(lowcut, highcut, fs, order=order)
    # Use filtfilt to avoid phase delay
    data_filt = sosfiltfilt(sos, data, axis=0)
    return data_filt


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


def run_main(nifti_file, mask, low, high, tr, output_file):
    # Get sampling rate
    fs = 1 / tr
    # Load mask and functional scan
    nifti = nb.load(nifti_file)
    nifti_data = nifti.get_fdata()
    nifti.uncache()
    mask = nb.load(mask).get_fdata() > 0
    nifti_data = convert_2d(mask, nifti_data)
    # Demean and norm data
    nifti_data = zscore(nifti_data)
    # Bandpass filter with Butterworth
    nifti_data = butter_bandpass_filter(nifti_data, low, high, fs)
    nifti_data_4d = convert_4d(mask, nifti_data)
    # Write output
    output_nifti(output_file, nifti, nifti_data_4d)


if __name__ == '__main__':
    """Bandpass filter w/ Butterworth Filter """
    parser = argparse.ArgumentParser(description='Normalize and bandpass filter cifti files')
    parser.add_argument('-n', '--nifti',
                        help='<Required> path to nifti file',
                        required=True,
                        type=str)
    parser.add_argument('-m', '--mask',
                        help='<Required> path to nifti file',
                        required=True,
                        type=str)
    parser.add_argument('-l', '--low_cut',
                        help='<Required> Lower bound of the bandpass filter',
                        required=True,
                        type=float)
    parser.add_argument('-u', '--high_cut',
                        help='<Required> Higher bound of the bandpass filter',
                        required=True,
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
    run_main(args_dict['nifti'], args_dict['mask'],
             args_dict['low_cut'], args_dict['high_cut'], 
             args_dict['tr'], args_dict['output_file'])
