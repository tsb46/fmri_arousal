import argparse
import numpy as np
import os

from scipy.signal import resample


def run_main(file, n_resample, output_file):
    # Load signal txt file
    data = np.loadtxt(file)
    data_resamp = resample(data, n_resample)
    np.savetxt(output_file, data_resamp)


if __name__ == '__main__':
    """Temporal filter w/ Butterworth Filter """
    parser = argparse.ArgumentParser(description='Normalize and temporal filter nifti files')
    parser.add_argument('-f', '--file_path',
                        help='<Required> path to file',
                        required=True,
                        type=str)
    parser.add_argument('-n', '--n_resample',
                        help='number of time points to resample to',
                        required=True,
                        type=int)
    parser.add_argument('-o', '--output_file',
                        help='output file path',
                        required=False,
                        default=os.getcwd(),
                        type=str)
    args_dict = vars(parser.parse_args())
    run_main(args_dict['file_path'], args_dict['n_resample'], args_dict['output_file'])
