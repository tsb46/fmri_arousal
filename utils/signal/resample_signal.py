import argparse
import numpy as np
import nibabel as nb
import os

from scipy.signal import resample


def run_main(file, func_scan, output_file):
    # Get number of time points from func_scan
    func = nb.load(func_scan)
    n_resample = func.header.get_data_shape()[3]
    # Load signal txt file
    data = np.loadtxt(file)
    data_resamp = resample(data, n_resample)
    np.savetxt(output_file, data_resamp)


if __name__ == '__main__':
    """Temporal filter w/ Butterworth Filter """
    parser = argparse.ArgumentParser(description='Resample (downsample) signal')
    parser.add_argument('-f', '--file_path',
                        help='<Required> path to file',
                        required=True,
                        type=str)
    parser.add_argument('-s', '--func_scan',
                        help='functional scan, used to determine number of time points to resample to',
                        required=True,
                        type=str)
    parser.add_argument('-o', '--output_file',
                        help='output file path',
                        required=False,
                        default=os.getcwd(),
                        type=str)
    args_dict = vars(parser.parse_args())
    run_main(args_dict['file_path'], args_dict['func_scan'], args_dict['output_file'])
