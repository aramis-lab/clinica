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


def create_mask(c1, c2, c3):
    import nibabel as nib
    import numpy as np
    from os import getcwd
    from os.path import join

    img_c1 = nib.load(c1)
    img_c2 = nib.load(c2)
    img_c3 = nib.load(c3)
    tissues = img_c1.get_data() + img_c2.get_data() + img_c3.get_data()
    background = 1.0 - tissues

    dshape = list(img_c1.shape)
    dshape.append(4)
    data = np.empty(dshape, dtype=np.float64)
    data[...,0] = np.array(img_c1.get_data())
    data[...,1] = np.array(img_c2.get_data())
    data[...,2] = np.array(img_c3.get_data())
    data[...,3] = np.array(background)

    out_mask = join(getcwd(), 'pvc_mask.nii')
    mask = nib.Nifti1Image(data, img_c1.affine, header=img_c1.header)
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
    region_mean = np.nanmean(np.where(region!=0, region, np.nan))

    data = pet.get_data() / region_mean

    suvr_pet_path = join(getcwd(), 'suvr_' + basename(pet_image))

    suvr_pet = nib.Nifti1Image(data, pet.affine, header=pet.header)
    nib.save(suvr_pet, suvr_pet_path)

    return suvr_pet_path
