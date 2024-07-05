import argparse
import matplotlib.pyplot as plt
import nibabel as nb 
import numpy as np
import os
import pickle
import subprocess as sbp

from pathlib import Path
from nilearn import plotting, image

# Inspired by:
# https://gist.github.com/tknapen/0ebf22353130635d899eafba75de31f5


def split_string(xyz_str):
    return [int(x) for x in xyz_str.split(',')]


def run_main(input_nifti, write_directory, color_map, cmax, cmin, 
             xyz_coords, figure_dimensions, tr_label, cpca_res, 
             frame_rate):
    x, y, z = split_string(xyz_coords)
    fig_x, fig_y = split_string(figure_dimensions)
    if tr_label is not None:
        tr_0, tr_max = split_string(tr_label)
        labels = np.arange(tr_0, tr_max+1) 
        fig_title = [f'TR: {l}' for l in labels]
    elif cpca_res is not None:
        cpca_res = pickle.load(open(cpca_res, 'rb'))
        t_phase = cpca_res[1][0]
        fig_title = [f'$\phi$ = {round(p,1)}' for p in t_phase]
    else:
        fig_title = None

    input_nifti_base = Path(input_nifti).stem
    nifti_img = nb.load(input_nifti)
    for i, img in enumerate(image.iter_img(nifti_img)):
        if fig_title is None:
            title = None
        else:
            title = fig_title[i]
        f = plt.figure(figsize=(fig_x,fig_y))
        img_indx = str(i).zfill(3)
        plotting.plot_img(
            img,
            bg_img=None,
            colorbar=True, 
            threshold=0,
            vmax=cmax,
            vmin=cmin,
            draw_cross=False,
            cmap=color_map,
            title=title,
            cut_coords=(x,y,z),
            output_file=f'{write_directory}/{input_nifti_base}_{img_indx}.png',
            figure=f
        )
    # Write image to video
    write_animation_ffmpeg(input_nifti_base, write_directory, i, frame_rate)

def write_animation_ffmpeg(file_base, write_directory, img_n, frame_rate):
    # Write file
    cmd = """
    ffmpeg -y -r {0} -i {2}/{1}_%03d.png -vcodec libx264 -crf 25 -acodec aac -pix_fmt yuv420p {2}/{1}.mp4
    """.format(frame_rate, file_base, write_directory)
    sbp.call(cmd, shell=True)

    # Delete PNG files in python
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
                        help='max colorbar limit',
                        default=None,
                        type=float)
    parser.add_argument('-cmin', '--color_min',
                        help='min colorbar limit',
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
                        help='two integers (first TR, last TR), comma separated list with no spaces',
                        default=None, 
                        type=str)
    parser.add_argument('-cpca', '--cpca_res',
                        help='if cpca reconstruction is provided to input_nifti, the '
                        'cpca recon results .pkl file for labeling the phase samples. If '
                        'tr_label is provided, this argument is ignored.',
                        default=None,
                        type=str)
    parser.add_argument('-r', '--frame_rate',
                        help='frame rate',
                        default=5,
                        type=float)

    args_dict = vars(parser.parse_args())
    run_main(args_dict['input_nifti'], args_dict['write_directory'], 
             args_dict['color_map'], args_dict['color_max'], 
             args_dict['color_min'],args_dict['xyz_coords'],
             args_dict['figure_dimensions'], args_dict['tr_label'], 
             args_dict['cpca_res'], args_dict['frame_rate'])
