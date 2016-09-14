
import numpy as np
import nibabel as nib
from os import path
from scipy.linalg import norm


def load_data(subjects, data_dir='', mask=False):

    first = True
    for i in range(len(subjects)):
        subj_path = path.join(data_dir, subjects[i])
        subj = nib.load(subj_path)
        subj_data = subj.get_data().flatten()

        # Memory allocation for ndarray containing all data to avoid copying the array for each new subject
        if first:
            data = np.ndarray(shape=(len(subjects), subj_data.shape[0]), dtype=float, order='C')
            shape = subj.get_data().shape
            first = False

        data[i, :] = subj_data

    if mask:
        data_mask = (data != 0).sum(axis=0) != 0
        data = data[:, data_mask]

    return data, data_mask, shape


def revert_mask(weights, mask, shape):

    z = np.zeros(np.prod(shape))
    z[mask] = weights

    new_weights = np.reshape(z, shape)

    return new_weights


def weights_to_nifti(weights, template, output_filename):

    comparable_features = 2 * weights.flatten() / np.power(norm(weights.flatten()), 2)
    comparable_features = comparable_features.reshape(weights.shape)

    img = nib.load(template)
    affine = img.get_affine()

    new_img = nib.Nifti1Image(comparable_features, affine)
    nib.save(new_img, output_filename)

