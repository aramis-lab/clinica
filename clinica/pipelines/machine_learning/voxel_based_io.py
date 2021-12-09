import nibabel as nib
import numpy as np


def load_data(image_list, mask=True):
    """

    Args:
        image_list:
        mask:

    Returns:

    """
    data = None
    shape = None
    data_mask = None
    first = True

    for i in range(len(image_list)):
        subj = nib.load(image_list[i])
        subj_data = np.nan_to_num(subj.get_fdata(dtype="float32").flatten())

        # Memory allocation for ndarray containing all data to avoid copying the array for each new subject
        if first:
            data = np.ndarray(
                shape=(len(image_list), subj_data.shape[0]), dtype=float, order="C"
            )
            shape = subj.get_fdata(dtype="float32").shape
            first = False

        data[i, :] = subj_data

    if mask:
        data_mask = (data != 0).sum(axis=0) != 0
        data = data[:, data_mask]

    return data, shape, data_mask


def revert_mask(weights, mask, shape):
    """

    Args:
        weights:
        mask:
        shape:

    Returns:

    """

    z = np.zeros(np.prod(shape))
    z[mask] = weights

    new_weights = np.reshape(z, shape)

    return new_weights


def features_weights(image_list, dual_coefficients, sv_indices, scaler=None, mask=None):

    if len(sv_indices) != len(dual_coefficients):
        print("Length dual coefficients: " + str(len(dual_coefficients)))
        print("Length indices: " + str(len(sv_indices)))
        raise ValueError(
            "The number of support vectors indices and the number of coefficients must be the same."
        )

    if len(image_list) == 0:
        raise ValueError("The number of images must be greater than 0.")

    sv_images = [image_list[i] for i in sv_indices]

    shape = nib.load(sv_images[0]).get_fdata(dtype="float32").shape
    weights = np.zeros(shape)

    for i in range(len(sv_images)):
        subj = nib.load(sv_images[i])
        subj_data = np.nan_to_num(subj.get_fdata(dtype="float32"))

        if scaler is not None and mask is not None:
            subj_data = subj_data.flatten()[mask]
            subj_data = scaler.transform(subj_data)
            subj_data = revert_mask(subj_data, mask, shape)
        weights += dual_coefficients[i] * subj_data

    return weights


def weights_to_nifti(weights, image, output_filename):
    """

    Args:
        weights:
        image:
        output_filename:

    Returns:

    """

    # Normalize with 2-norm
    # features = 2 * weights / np.power(norm(weights.flatten(), 2), 2)

    # Normalize inf-norm
    features = weights / abs(weights).max()

    img = nib.load(image)
    canonical_img = nib.as_closest_canonical(img)
    hd = canonical_img.header
    qform = np.zeros((4, 4))

    for i in range(1, 4):
        qform[i - 1, i - 1] = hd["pixdim"][i]
        qform[i - 1, 3] = -1.0 * hd["pixdim"][i] * hd["dim"][i] / 2.0

    output_image = nib.Nifti1Image(features, qform)
    nib.save(output_image, output_filename)
