import nipype.pipeline.engine as pe
import nipype.interfaces.io as nio


def create_datagrabber(name, subjects, sessions, base_directory, template):

    dg = pe.Node(nio.DataGrabber(infields=['subject_id', 'session', 'subject_repeat', 'session_repeat'],
                                 outfields=['out_files']), name=name)
    dg.inputs.base_directory = base_directory
    dg.inputs.template = template
    dg.inputs.subject_id = subjects
    dg.inputs.session = sessions
    dg.inputs.subject_repeat = subjects
    dg.inputs.session_repeat = sessions
    dg.inputs.sort_filelist = False
    return dg


def create_tissues_datagrabber(name, subject, session, tissues, base_directory, template):

    dg = pe.Node(nio.DataGrabber(infields=['subject_id', 'session', 'tissues', 'subject_repeat', 'session_repeat'],
                                 outfields=['out_files']), name=name)
    dg.inputs.base_directory = base_directory
    dg.inputs.template = template
    dg.inputs.subject_id = [subject] * len(tissues)
    dg.inputs.session = [session] * len(tissues)
    dg.inputs.tissues = [str(x) for x in tissues]
    dg.inputs.subject_repeat = [subject] * len(tissues)
    dg.inputs.session_repeat = [session] * len(tissues)
    dg.inputs.sort_filelist = False
    return dg


def create_binary_mask(tissues, threshold=0.3):
    import nibabel as nib
    import numpy as np
    from os import getcwd
    from os.path import join

    if len(tissues) == 0:
        raise RuntimeError('The length of the list of tissues must be greater than zero.')

    img_0 = nib.load(tissues[0])
    shape = list(img_0.get_data().shape)

    data = np.zeros(shape=shape)

    for image in tissues:
        data = data + nib.load(image).get_data()

    data = (data > threshold) * 1.0
    out_mask = join(getcwd(), 'binary_mask.nii')

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

    shape.append(len(tissues) + 1)
    data = np.empty(shape, dtype=np.float64)

    for i in range(len(tissues)):
        image = tissues[i]
        data[..., i] = np.array(image.get_data())
        background = background + nib.load(image).get_data()

    background = 1.0 - background
    data[..., len(tissues)] = np.array(background)

    out_mask = join(getcwd(), 'pvc_mask.nii')
    mask = nib.Nifti1Image(data, img_0.affine, header=img_0.header)
    nib.save(mask, out_mask)
    return out_mask




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
