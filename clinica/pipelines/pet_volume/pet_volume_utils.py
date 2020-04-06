# coding: utf8


def read_psf_information(psf_tsv, subject_ids, session_ids):
    import os
    from pandas.io.parsers import read_csv

    if not os.path.isfile(psf_tsv):
        raise FileNotFoundError('Could not find the psf_tsv file %s' % psf_tsv)
    try:
        psf_df = read_csv(psf_tsv, sep='\t')
    except (IOError, UnicodeDecodeError):
        raise RuntimeError('An error while reading %s happened' % psf_tsv)

    if psf_df.shape[0] != len(subject_ids):
        raise ValueError('The number of rows in fwhm_tsv file must match the number of subject-session pairs.')

    if any(elem not in ['participant_id', 'session_id', 'fwhm_x', 'fwhm_y', 'fwhm_z'] for elem in list(psf_df.columns)):
        raise IOError('The file %s must contain the following columns (separated by tabulations):\n'
                      'participant_id, session_id, fwhm_x, fwhm_y, fwhm_z\n'
                      'But we found:\n'
                      '%s\n'
                      'Pay attention to the spaces (there should be none).' %
                      (psf_tsv, str(list(psf_df.columns))))

    subjects_fwhm = list(psf_df.participant_id)
    sessions_fwhm = list(psf_df.session_id)
    idx_reordered = []
    for i, sub in enumerate(subject_ids):
        current_ses = session_ids[i]
        idx_sub = [j for j in range(len(subjects_fwhm)) if sub == subjects_fwhm[j] and current_ses == sessions_fwhm[j]]
        if len(idx_sub) == 0:
            raise RuntimeError('Subject %s with session %s that you want to proceed was not found '
                               'in the TSV file containing PSF specifications (%s).' %
                               (sub, current_ses, psf_tsv))
        if len(idx_sub) > 1:
            raise RuntimeError('Subject %s with session %s that you want to proceed was found multiple times '
                               'in the TSV file containing PSF specifications (%s).' %
                               (sub, current_ses, psf_tsv))
        idx_reordered.append(idx_sub[0])

    fwhm_x = list(psf_df.fwhm_x)
    fwhm_y = list(psf_df.fwhm_y)
    fwhm_z = list(psf_df.fwhm_z)
    iterables_fwhm = [[fwhm_x[i], fwhm_y[i], fwhm_z[i]] for i in idx_reordered]
    return iterables_fwhm


def init_input_node(pet_nii):
    import datetime
    import nibabel as nib
    from colorama import Fore
    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.stream import cprint
    from clinica.utils.ux import print_begin_image

    # Extract image ID
    image_id = get_subject_id(pet_nii)

    # Check that the PET file is a 3D volume
    img = nib.load(pet_nii)
    if len(img.shape) == 4:
        now = datetime.datetime.now().strftime('%H:%M:%S')
        error_msg = '%s[%s] Error: Clinica does not handle 4D volumes for %s%s'\
                    % (Fore.RED, now, image_id.replace('_', ' | '), Fore.RESET)
        cprint(error_msg)
        raise NotImplementedError(error_msg)

    # Print begin message
    print_begin_image(image_id)

    return pet_nii


def create_binary_mask(tissues, threshold=0.3):
    """

    Args:
        tissues:
        threshold:

    Returns:

    """
    import nibabel as nib
    import numpy as np
    from os import getcwd
    from os.path import join, basename

    if len(tissues) == 0:
        raise RuntimeError('The length of the list of tissues must be greater than zero.')

    img_0 = nib.load(tissues[0])
    shape = list(img_0.get_data().shape)

    data = np.zeros(shape=shape)

    for image in tissues:
        data = data + nib.load(image).get_data()

    data = (data > threshold) * 1.0
    out_mask = join(getcwd(), basename(tissues[0]) + '_brainmask.nii')

    mask = nib.Nifti1Image(data, img_0.affine, header=img_0.header)
    nib.save(mask, out_mask)
    return out_mask


def apply_binary_mask(image, binary_mask):
    import nibabel as nib
    from os import getcwd
    from os.path import join, basename

    original_image = nib.load(image)
    mask = nib.load(binary_mask)

    data = original_image.get_data() * mask.get_data()

    masked_image_path = join(getcwd(), 'masked_' + basename(image))
    masked_image = nib.Nifti1Image(data, original_image.affine, header=original_image.header)
    nib.save(masked_image, masked_image_path)
    return masked_image_path


def create_pvc_mask(tissues):

    import nibabel as nib
    import numpy as np
    from os import getcwd
    from os.path import join

    if len(tissues) == 0:
        raise RuntimeError('The length of the list of tissues must be greater than zero.')

    img_0 = nib.load(tissues[0])
    shape = img_0.get_data().shape
    background = np.zeros(shape=shape)

    shape += tuple([len(tissues) + 1])
    data = np.empty(shape=shape, dtype=np.float64)

    for i in range(len(tissues)):
        image = nib.load(tissues[i])
        data[..., i] = np.array(image.get_data())
        background = background + image.get_data()

    background = 1.0 - background
    data[..., len(tissues)] = np.array(background)

    out_mask = join(getcwd(), 'pvc_mask.nii')
    mask = nib.Nifti1Image(data, img_0.affine, header=img_0.header)
    nib.save(mask, out_mask)
    return out_mask


def pet_pvc_name(pet_image, pvc_method):
    from os.path import basename
    pet_pvc_path = 'pvc-' + pvc_method.lower() + '_' + basename(pet_image)
    return pet_pvc_path


def normalize_to_reference(pet_image, region_mask):
    import nibabel as nib
    import numpy as np
    from os import getcwd
    from os.path import basename, join

    pet = nib.load(pet_image)
    ref = nib.load(region_mask)

    region = np.multiply(pet.get_data(), ref.get_data())
    region_mean = np.nanmean(np.where(region != 0, region, np.nan))

    data = pet.get_data() / region_mean

    suvr_pet_path = join(getcwd(), 'suvr_' + basename(pet_image))

    suvr_pet = nib.Nifti1Image(data, pet.affine, header=pet.header)
    nib.save(suvr_pet, suvr_pet_path)

    return suvr_pet_path


def atlas_statistics(in_image, in_atlas_list):
    """
    For each atlas name provided it calculates for the input image the mean
    for each region in the atlas and saves it to a TSV file.

    Args:
        in_image: A Nifti image
        in_atlas_list: List of names of atlas to be applied

    Returns:
        List of paths to tsv files
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
                out_atlas_statistics = abspath(join(getcwd(), base + '_space-' + atlas + '_statistics.tsv'))
                statistics_on_atlas(in_image, atlas_class(), out_atlas_statistics)
                atlas_statistics_list.append(out_atlas_statistics)
                break

    return atlas_statistics_list


def pet_container_from_filename(pet_filename):
    import re
    from os.path import join
    m = re.search(r'(sub-[a-zA-Z0-9]+)_(ses-[a-zA-Z0-9]+)_', pet_filename)

    if m is None:
        raise ValueError('Input filename is not in a BIDS or CAPS compliant format. It does not contain the subject' +
                         ' and session information.')

    subject = m.group(1)
    session = m.group(2)

    return join('subjects', subject, session, 'pet/preprocessing')


def get_from_list(in_list, index):
    return in_list[index]
