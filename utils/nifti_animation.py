import argparse
import matplotlib.pyplot as plt
import nibabel as nb 
import numpy as np
import os
import subprocess as sbp

from pathlib import Path
from nilearn import plotting, image

# Inspired by:
# https://gist.github.com/tknapen/0ebf22353130635d899eafba75de31f5


def split_string(xyz_str):
    return [int(x) for x in xyz_str.split(',')]


def run_main(input_nifti, write_directory, color_map, color_max, xyz_coords, figure_dimensions, tr_label, frame_rate):
    x, y, z = split_string(xyz_coords)
    fig_x, fig_y = split_string(figure_dimensions)
    if tr_label is not None:
        tr_0, tr_max = split_string(tr_label)
        labels = np.arange(tr_0, tr_max+1) 

    input_nifti_base = Path(input_nifti).stem
    nifti_img = nb.load(input_nifti)
    for i, img in enumerate(image.iter_img(nifti_img)):
        if tr_label is not None:
            fig_title = f'TR: {labels[i]}'
        else:
            fig_title = None
        f = plt.figure(figsize=(fig_x,fig_y))
        img_indx = str(i).zfill(3)
        plotting.plot_stat_map(img,
                               bg_img=None,
                               colorbar=False, 
                               threshold=0,
                               vmax=color_max,
                               draw_cross=False,
                               title=fig_title,
                               cut_coords=(x,y,z),
                               output_file=f'{write_directory}/{input_nifti_base}_{img_indx}.png',
                               figure=f)
    # Write image to video
    write_animation_ffmpeg(input_nifti_base, write_directory, i, frame_rate)

def write_animation_ffmpeg(file_base, write_directory, img_n, frame_rate):
    # Write file
    cmd = """
    ffmpeg -y -r {0} -i {2}/{1}_%03d.png -vcodec libx264 -crf 25 -acodec aac -pix_fmt yuv420p {2}/{1}.mp4
    """.format(frame_rate, file_base, write_directory)
    sbp.call(cmd, shell=True)

    # Delete PNG files in python, don't trust my self with shell :) 
    for i in range(img_n+1):
        img_indx = str(i).zfill(3)
        os.remove(f'{write_directory}/{file_base}_{img_indx}.png')




if __name__ == '__main__':
    """create nifti animation"""
    parser = argparse.ArgumentParser(description='create nifti animation')
    parser.add_argument('-n', '--input_nifti',
                        help='<Required> file path to nifti file for animation',
                        required=True,
                        type=str)
    parser.add_argument('-w', '--write_directory',
                        help='path to write static pics and animation',
                        required=True,
                        type=str)
    parser.add_argument('-c', '--color_map',
                        help='matplotlib colormap specified as string',
                        default='cold_hot',
                        type=str)
    parser.add_argument('-cmax', '--color_max',
                        help='max and min (symmetric) colorbar limit',
                        default=None,
                        type=float)
    parser.add_argument('-xyz', '--xyz_coords',
                        help='x, y and z coordinates of 3D slices to display. comma separated list with no spaces',
                        default='0,-18,15', # orthogonal slices for 3D MNI template
                        type=str)
    parser.add_argument('-f', '--figure_dimensions',
                        help='x and y figure dimension size. comma separated list with no spaces',
                        default='16,6', # orthogonal slices for 3D MNI template
                        type=str)
    parser.add_argument('-t', '--tr_label',
                        help='two integers, comma separated list with no spaces',
                        default=None, 
                        type=str)
    parser.add_argument('-r', '--frame_rate',
                        help='frame rate',
                        default=10,
                        type=float)

    args_dict = vars(parser.parse_args())
    run_main(args_dict['input_nifti'], args_dict['write_directory'], 
             args_dict['color_map'], args_dict['color_max'], 
             args_dict['xyz_coords'], args_dict['figure_dimensions'],
             args_dict['tr_label'], args_dict['frame_rate'])
