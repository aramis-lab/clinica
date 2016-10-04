
import numpy as np
import nibabel as nib
from numpy.linalg import norm


def load_data(image_list, mask=True):
    data = None
    shape = None
    data_mask = None
    first = True

    for i in range(len(image_list)):
        subj = nib.load(image_list[i])
        subj_data = subj.get_data().flatten()

        # Memory allocation for ndarray containing all data to avoid copying the array for each new subject
        if first:
            data = np.ndarray(shape=(len(image_list), subj_data.shape[0]), dtype=float, order='C')
            shape = subj.get_data().shape
            first = False

        data[i, :] = subj_data

    if mask:
        data_mask = (data != 0).sum(axis=0) != 0
        data = data[:, data_mask]

    return data, shape, data_mask


def revert_mask(weights, mask, shape):

    z = np.zeros(np.prod(shape))
    z[mask] = weights

    new_weights = np.reshape(z, shape)

    return new_weights


def weights_to_nifti(weights, template, output_filename):

    # Normalize with 2-norm
    # comparable_features = 2 * weights.flatten() / np.power(norm(weights.flatten(), 2), 2)

    # Normalize inf-norm
    comparable_features = weights.flatten() / max(abs(weights.flatten()))

    comparable_features = comparable_features.reshape(weights.shape)

    img = nib.load(template)
    img.get_data()[:] = comparable_features
    nib.save(img, output_filename)

    # Recommended way, resulting images not working with Anatomist on Linux.
    # affine = img.get_affine()
    # new_img = nib.Nifti1Image(comparable_features, affine)
    # nib.save(new_img, output_filename)


