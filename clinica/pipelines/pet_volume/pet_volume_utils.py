def init_input_node(pet_nii):
    import nibabel as nib

    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.stream import cprint
    from clinica.utils.ux import print_begin_image

    # Extract image ID
    image_id = get_subject_id(pet_nii)

    # Check that the PET file is a 3D volume
    img = nib.load(pet_nii)
    if len(img.shape) == 4:
        error_msg = f"Clinica does not handle 4D volumes for {image_id.replace('_', ' | ')}"
        cprint(error_msg, lvl="error")
        raise NotImplementedError(error_msg)

    # Print begin message
    print_begin_image(image_id)

    return pet_nii


def create_binary_mask(tissues, threshold=0.3):
    from os import getcwd
    from os.path import basename, join

    import nibabel as nib
    import numpy as np

    if len(tissues) == 0:
        raise RuntimeError(
            "The length of the list of tissues must be greater than zero."
        )

    img_0 = nib.load(tissues[0])
    shape = list(img_0.get_fdata(dtype="float32").shape)

    data = np.zeros(shape=shape)

    for image in tissues:
        data = data + nib.load(image).get_fdata(dtype="float32")

    data = (data > threshold) * 1.0
    out_mask = join(getcwd(), basename(tissues[0]) + "_brainmask.nii")

    mask = nib.Nifti1Image(data, img_0.affine, header=img_0.header)
    nib.save(mask, out_mask)
    return out_mask


def apply_binary_mask(image, binary_mask):
    from os import getcwd
    from os.path import basename, join

    import nibabel as nib

    original_image = nib.load(image)
    mask = nib.load(binary_mask)

    data = original_image.get_fdata(dtype="float32") * mask.get_fdata(dtype="float32")

    masked_image_path = join(getcwd(), "masked_" + basename(image))
    masked_image = nib.Nifti1Image(
        data, original_image.affine, header=original_image.header
    )
    nib.save(masked_image, masked_image_path)
    return masked_image_path


def create_pvc_mask(tissues):

    from os import getcwd
    from os.path import join

    import nibabel as nib
    import numpy as np

    if len(tissues) == 0:
        raise RuntimeError(
            "The length of the list of tissues must be greater than zero."
        )

    img_0 = nib.load(tissues[0])
    shape = img_0.get_fdata(dtype="float32").shape
    background = np.zeros(shape=shape)

    shape += tuple([len(tissues) + 1])
    data = np.empty(shape=shape, dtype=np.float64)

    for i in range(len(tissues)):
        image = nib.load(tissues[i])
        data[..., i] = np.array(image.get_fdata(dtype="float32"))
        background = background + image.get_fdata(dtype="float32")

    background = 1.0 - background
    data[..., len(tissues)] = np.array(background)

    out_mask = join(getcwd(), "pvc_mask.nii")
    mask = nib.Nifti1Image(data, img_0.affine, header=img_0.header)
    nib.save(mask, out_mask)
    return out_mask


def pet_pvc_name(pet_image, pvc_method):
    from os.path import basename

    pet_pvc_path = "pvc-" + pvc_method.lower() + "_" + basename(pet_image)
    return pet_pvc_path


def normalize_to_reference(pet_image, region_mask):
    from os import getcwd
    from os.path import basename, join

    import nibabel as nib
    import numpy as np

    pet = nib.load(pet_image)
    ref = nib.load(region_mask)

    region = np.multiply(pet.get_fdata(dtype="float32"), ref.get_fdata(dtype="float32"))
    region_mean = np.nanmean(np.where(region != 0, region, np.nan))

    data = pet.get_fdata(dtype="float32") / region_mean

    suvr_pet_path = join(getcwd(), "suvr_" + basename(pet_image))

    suvr_pet = nib.Nifti1Image(data, pet.affine, header=pet.header)
    nib.save(suvr_pet, suvr_pet_path)

    return suvr_pet_path


def atlas_statistics(in_image, in_atlas_list):
    """Generate regional measure from atlas_list in TSV files.

    For each atlas name provided it calculates for the input image the mean
    for each region in the atlas and saves it to a TSV file.

    Args:
        in_image: A Nifti image
        in_atlas_list: List of names of atlas to be applied

    Returns:
        List of paths to TSV files
    """
    from os import getcwd
    from os.path import abspath, join

    from nipype.utils.filemanip import split_filename

    from clinica.utils.atlas import AtlasAbstract
    from clinica.utils.statistics import statistics_on_atlas

    orig_dir, base, ext = split_filename(in_image)
    atlas_classes = AtlasAbstract.__subclasses__()
    atlas_statistics_list = []
    for atlas in in_atlas_list:
        for atlas_class in atlas_classes:
            if atlas_class.get_name_atlas() == atlas:
                out_atlas_statistics = abspath(
                    join(getcwd(), base + "_space-" + atlas + "_statistics.tsv")
                )
                statistics_on_atlas(in_image, atlas_class(), out_atlas_statistics)
                atlas_statistics_list.append(out_atlas_statistics)
                break

    return atlas_statistics_list


def get_from_list(in_list, index):
    return in_list[index]
