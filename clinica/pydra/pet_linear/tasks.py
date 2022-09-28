import typing as ty
import nibabel as nib
import numpy as np
from pydra.mark import annotate, task


@task
@annotate({"return": {"transforms_list": ty.List[str]}})
def concatenate_transforms(pet_to_t1w_tranform: str, t1w_to_mni_tranform: str) -> ty.List[str]:
    """Concatenate two input transformation files into a list.

    Parameters
    ----------
    pet_to_t1w_tranform : str
        First transformation to apply.

    t1w_to_mni_tranform : str
        Second transformation to apply.

    Returns
    -------
    List[str] :
        Both transform files path in a list.
    """
    return [t1w_to_mni_tranform, pet_to_t1w_tranform]


@task
@annotate({"return": {"output_img": nib.nifti1.Nifti1Image}})
def suvr_normalization(input_img: str, norm_img: str, ref_mask: str) -> nib.nifti1.Nifti1Image:
    """Normalize the input image according to the reference region.
            
    It uses nilearn `resample_to_img` and scipy `trim_mean` functions.
    This function is different than the one in other PET pipelines
    because there is a downsampling step.

    Parameters
    ----------
    input_img : str
        Path to the image to be processed.

    norm_img : str
        Path to the image used to compute the mean of the reference region.

    ref_mask : str
        Path to the mask of the reference region.

    Returns
    -------
    output_img : Nifti1Image
        Normalized nifty image
    """
    import os
    from nilearn.image import resample_to_img
    from scipy.stats import trim_mean

    pet = nib.load(input_img)
    norm = nib.load(norm_img)
    mask = nib.load(ref_mask)
    
    # Downsample the pet image used for normalization so we can multiply it with the mask
    ds_img = resample_to_img(norm, mask, interpolation="nearest")
    
    # Compute the mean of the region
    region = np.multiply(ds_img.get_fdata(), mask.get_fdata(dtype="float32"))
    array_region = np.where(region != 0, region, np.nan).flatten()
    region_mean = trim_mean(array_region[~np.isnan(array_region)], 0.1)

    from clinica.utils.stream import cprint
    cprint(region_mean)
    
    # Divide the value of the image voxels by the computed mean
    data = pet.get_fdata(dtype="float32") / region_mean

    # Create and save the normalized image
    output_img = os.path.join(
            os.getcwd(),
            os.path.basename(input_img).split(".nii")[0] + "_suvr_normalized.nii.gz",
    )
    normalized_img = nib.Nifti1Image(data, pet.affine, header=pet.header)
    normalized_img.to_filename(output_img)
    return output_img


@task
@annotate({"return": {"output_img": nib.nifti1.Nifti1Image}})
def crop_nifti(input_img: str, ref_crop: str):
    """Crop input image based on the reference.

    It uses nilearn `resample_to_img` function.
                
    Parameters
    ----------
    input_img : str
        Image to be processed.

    ref_img : str
        Template used to crop the image.

    Returns
    -------
    output_img : Nifti1Image
        Cropped image on disk.
    """
    import os
    from nilearn.image import resample_to_img
    
    basedir = os.getcwd()
    # resample the individual MRI into the cropped template image
    crop_img = resample_to_img(input_img, ref_crop, force_resample=True)
    output_img = os.path.join(
            basedir, os.path.basename(input_img).split(".nii")[0] + "_cropped.nii.gz"
    )
    crop_img.to_filename(output_img)
    return output_img
