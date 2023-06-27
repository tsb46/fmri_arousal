import os
import nibabel as nb
import shutil

from nipype.interfaces import fsl
from nipype.interfaces.utility import Function
from utils.load_write import get_fp_base

# Ensure output is .nii.gz
fsl.FSLCommand.set_default_output_type('NIFTI_GZ')


def apply_mask(fp, fp_out, mask):
    # apply mask to functional image 
    applymask = fsl.ApplyMask()
    applymask.inputs.in_file = fp
    applymask.inputs.mask_file = mask
    applymask.inputs.out_file = fp_out
    applymask_res = applymask.run()


def apply_transform_mask(fp, fp_out, ref, mat):
    # apply affine transform to functional
    applyxfm = fsl.ApplyXFM()
    applyxfm.inputs.in_file = fp
    applyxfm.inputs.in_matrix_file = mat 
    applyxfm.inputs.reference = ref
    applyxfm.inputs.out_file = fp_out 
    applyxfm.inputs.apply_xfm = True
    applyxfm.inputs.interp = 'nearestneighbour' 
    applyxfm_res = applyxfm.run()
    # Flirt saves output matrix in base directory 
    # move to results directory
    os.remove(applyxfm_res.outputs.out_matrix_file)

def bet(fp, fp_out):
    # BET - Skullstrip anatomical Image
    bet_anat = fsl.BET(frac=0.25, robust=True, mask=True)
    bet_anat.inputs.in_file = fp
    bet_anat.inputs.out_file = fp_out
    bet_anat_res = bet_anat.run()


def concat_transform(fp_func2struct, fp_flirt, fp_out):
    # Concatenate affine transform matrices (func2struct & struct2MNI)
    convertxfm = fsl.ConvertXFM(concat_xfm=True)
    convertxfm.inputs.in_file = fp_func2struct
    convertxfm.inputs.in_file2=fp_flirt
    convertxfm.inputs.out_file=fp_out
    convertxfm_res = convertxfm.run()


def coregister(fp_mean, fp_reorient, fp_bet, fp_wm, fp_out):
    fp_out_base = get_fp_base(fp_out)
    # Coregister functional with T1w
    epireg = fsl.EpiReg()
    epireg.inputs.epi = fp_mean
    epireg.inputs.t1_head=fp_reorient
    epireg.inputs.t1_brain=fp_bet
    # epireg expects the wmseg output as a suffix to the epi image (weird)
    # rename for now
    wmseg = f'{fp_out_base}_fast_wmseg.nii.gz'
    shutil.copyfile(fp_wm, wmseg)
    epireg.inputs.wmseg = wmseg
    epireg.inputs.out_base = fp_out_base
    epireg_res = epireg.run()
    fp_func2struct = epireg_res.outputs.epi2str_mat
    return fp_func2struct


def fast(fp, fp_out):
    # FAST - Image Segmentation
    fast = fsl.FAST()
    fast.inputs.in_files = fp
    fast.inputs.out_basename = fp_out
    # Nipype FAST issue with writing out tissue_class_map - Ignore
    # https://github.com/nipy/nipype/issues/3311
    try: 
        fast_res = fast.run()
    except FileNotFoundError:
        fast_out = fp_out
    fp_out_base = get_fp_base(fp_out)
    fp_out_wm = f'{fp_out_base}_pve_2'
    return fp_out_wm


def first_vol(fp, fp_out):
    first_vol = fsl.ExtractROI()
    first_vol.inputs.in_file = fp
    first_vol.inputs.roi_file = fp_out
    first_vol.inputs.t_min=0
    first_vol.inputs.t_size=1
    first_vol.run()


def flirt(fp, fp_out, fp_out_mat):
    # FLIRT affine registration to MNI template
    flirt = fsl.FLIRT()
    flirt.inputs.in_file = fp
    flirt.inputs.reference = f'{os.environ["FSLDIR"]}/data/standard/MNI152_T1_2mm_brain.nii.gz'
    flirt.inputs.out_file = fp_out
    flirt.out_matrix_file = fp_out_mat
    flirt_res = flirt.run()
    # Flirt saves output matrix in base directory (seems to be an issue related to the FAST issue above), 
    # move to results directory
    os.rename(flirt_res.outputs.out_matrix_file, fp_out_mat)


def fnirt(fp, fp_affine, fp_out, fp_out_coef):
    # FNIRT non-linear registration
    fnirt = fsl.FNIRT()
    fnirt.inputs.in_file = fp
    fnirt.inputs.ref_file = f'{os.environ["FSLDIR"]}/data/standard/MNI152_T1_2mm.nii.gz'
    fnirt.inputs.affine_file = fp_affine
    fnirt.inputs.config_file='T1_2_MNI152_2mm'
    fnirt.inputs.warped_file=fp_out
    fnirt.inputs.fieldcoeff_file = fp_out_coef
    fp_out_base = get_fp_base(fp_out)
    fnirt.inputs.log_file = f'{fp_out_base}_log.txt'
    fnirt_res = fnirt.run()


def invert_transform(in_mat, out_mat):
    # invert affine transform
    invertxfm = fsl.ConvertXFM()
    invertxfm.inputs.in_file = in_mat
    invertxfm.inputs.invert_xfm = True
    invertxfm.inputs.out_file = out_mat
    invertxfm.run()

def mcflirt(fp, fp_out):
    # McFLIRT Motion Correction
    fp_out_base = get_fp_base(fp_out)
    func_file_meanvol = f'{fp_out_base}_mean.nii.gz'
    mcflirt = fsl.MCFLIRT(mean_vol=True, save_plots=True)
    mcflirt.inputs.in_file = fp
    mcflirt.inputs.out_file = fp_out
    mcflirt_res = mcflirt.run()
    # weird renaming of mean vol, rename
    os.rename(mcflirt_res.outputs.mean_img, func_file_meanvol)
    return func_file_meanvol


def resample_func(fp, fp_out, mask):
    # FLIRT resample functional scan to 3mm MNI
    flirt = fsl.FLIRT()
    flirt.inputs.in_file = fp
    flirt.inputs.reference = mask
    flirt.inputs.out_file = fp_out
    flirt.inputs.apply_xfm = True
    flirt.inputs.uses_qform = True
    flirt_res = flirt.run()
    # Flirt saves output matrix in base directory (seems to be an issue related to the FAST issue above), 
    # move to results directory
    os.remove(flirt_res.outputs.out_matrix_file)


def robustfov(fp, fp_out):
    crop = fsl.RobustFOV()
    crop.inputs.in_file = fp
    crop.inputs.out_roi = fp_out
    crop.inputs.out_transform = f'{get_fp_base(fp)}_roi.mat'
    crop_res = crop.run()


def reorient(fp, fp_out):
    # Reorient 2 standard
    reorient = fsl.utils.Reorient2Std()
    reorient.inputs.in_file = fp
    reorient.inputs.out_file = fp_out
    reorient_res = reorient.run()


def slicetime(fp, fp_out, st_fp, tr):
    # # Slice time correction
    slicetimer = fsl.SliceTimer(custom_timings=st_fp, 
                                time_repetition=tr) 
    slicetimer.inputs.in_file = fp
    slicetimer.inputs.out_file = fp_out
    slicetimer_res = slicetimer.run()


def spatial_smooth(fp, fp_out, fwhm=5.0):
    # 5mm FWHM isotropic smoothing
    smooth = fsl.Smooth(fwhm=fwhm)
    smooth.inputs.in_file = fp
    smooth.inputs.smoothed_file=fp_out
    smooth_res = smooth.run()


def trim_vol(fp, fp_out, n_trim):
    # Trim first (+) or last (-) N volumes
    """
    if the integer, n, provided is positive, trim off the first 
    n volumes. If negative, trim off the last n volumes.
    """
    if n_trim >= 0:
        trim = fsl.ExtractROI(t_min=n_trim, t_size=-1)
    else:
        n_end = nb.load(fp).shape[-1]
        n_end -= abs(n_trim)
        trim = fsl.ExtractROI(t_min=0, t_size=n_end)
    trim.inputs.in_file = fp
    trim.inputs.roi_file = fp_out
    trim_res = trim.run()


def wm_thres(fp, fp_out):
    # Threshold white matter partial volume
    wm_thres = fsl.Threshold(thresh=0.5, args='-bin')
    wm_thres.inputs.in_file = fp
    wm_thres.inputs.out_file = fp_out
    wm_thres_res = wm_thres.run()


def warp_func(fp, fp_affine, fp_coef, fp_out, mask):
    # Warp functional to MNI space
    applywarp = fsl.ApplyWarp()
    applywarp.inputs.in_file = fp
    applywarp.inputs.ref_file = mask
    applywarp.inputs.premat=fp_affine
    applywarp.inputs.field_file=fp_coef
    applywarp.inputs.out_file=fp_out
    applywarp_res = applywarp.run()

