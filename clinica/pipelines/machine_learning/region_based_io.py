import nibabel as nib
import numpy as np
import pandas as pd


def load_data(image_list, subjects):
    """

    Args:
        image_list:
        subjects:

    Returns:

    """

    subj_average = []
    all_vector = np.array([])
    read_file = pd.read_csv(image_list[0], sep="\t", usecols=[2], header=0)
    read_file = read_file.mean_scalar
    data = np.zeros((len(subjects), len(read_file)))

    for i in range(len(image_list)):
        tsv_file = pd.read_csv(image_list[i], sep="\t", usecols=[2], header=0)
        subj_average = tsv_file.mean_scalar
        all_vector = np.append(all_vector, subj_average)
    data_temp = np.split(all_vector, len(image_list))

    for i in range(len(image_list)):
        for j in range(len(subj_average)):
            data[i][j] = data_temp[i][j]
    return data


def features_weights(image_list, dual_coefficients, sv_indices, scaler=None):
    """

    Args:
        image_list:
        dual_coefficients:
        sv_indices:
        scaler:

    Returns:

    """

    if len(sv_indices) != len(dual_coefficients):
        raise ValueError(
            f"The number of support vectors indices and the number of coefficients must be the same.\n"
            f"- Number of dual coefficients: {len(dual_coefficients)}\n"
            f"- Number of indices:: {len(sv_indices)}\n"
        )

    if len(image_list) == 0:
        raise ValueError("The number of images must be greater than 0.")

    sv_images = [image_list[i] for i in sv_indices]

    shape = pd.read_csv(sv_images[0], sep="\t", usecols=[2], header=0)
    weights = np.zeros(len(shape))

    for i in range(len(sv_images)):
        subj = pd.read_csv(sv_images[i], sep="\t", usecols=[2], header=0)
        subj_data = subj.mean_scalar
        weights += dual_coefficients[i] * subj_data

    return weights


def weights_to_nifti(weights, atlas, output_filename):
    """

    Args:
        atlas:
        weights:
        output_filename:

    Returns:

    """
    from clinica.utils.atlas import AtlasAbstract

    atlas_path = None
    atlas_classes = AtlasAbstract.__subclasses__()
    for atlas_class in atlas_classes:
        if atlas_class.get_name_atlas() == atlas:
            atlas_path = atlas_class.get_atlas_labels()

    if not atlas_path:
        raise ValueError("Atlas path not found for atlas name " + atlas)

    atlas_image = nib.load(atlas_path)
    atlas_data = atlas_image.get_fdata(dtype="float32")
    labels = list(set(atlas_data.ravel()))
    output_image_weights = np.array(atlas_data, dtype="f")
    for i, n in enumerate(labels):
        index = np.array(np.where(atlas_data == n))
        output_image_weights[index[0, :], index[1, :], index[2, :]] = weights[i]

    affine = atlas_image.get_affine()
    header = atlas_image.get_header()
    output_image = nib.Nifti1Image(output_image_weights, affine, header)
    nib.save(output_image, output_filename)
