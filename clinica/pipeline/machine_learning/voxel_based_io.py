
import numpy as np
import pandas as pd
import nibabel as nib
from numpy.linalg import norm
import csv
from ntpath import basename, splitext
from scipy.spatial.distance import squareform


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

    # Not recommendeed way
    img.get_data()[:] = comparable_features
    nib.save(img, output_filename)

    # Recommended way, resulting images not working with Anatomist on Linux.
    # affine = img.get_affine()
    # new_img = nib.Nifti1Image(comparable_features, affine)
    # nib.save(new_img, output_filename)


def save_subjects_prediction(subjects, diagnosis, y, y_hat, output_file):

    with open(output_file, 'w') as csvfile:

        fieldnames = ['Subject', 'Diagnose', 'Class_label', 'Predicted_label', 'Correct']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for i in range(len(subjects)):
            s = splitext(basename(subjects[i]))[0]
            dx = diagnosis[i]
            writer.writerow({'Subject': s,
                             'Diagnose': dx,
                             'Class_label': y[i],
                             'Predicted_label': int(y_hat[i]),
                             'Correct': int(y[i] == y_hat[i])
                             })


def results_to_csv(results, diagnose_list, output_file):

    with open(output_file, 'w') as output:
        s = 'Balanced Accuracy\n'
        balanced_accuracy_list = [round(res[1]['balanced_accuracy'], 2) for res in sorted(results.items())]
        df = pd.DataFrame(squareform(balanced_accuracy_list), index=diagnose_list, columns=diagnose_list)
        print df
        s += df.to_csv()

        s += '\nAccuracy\n'
        accuracy_list = [round(res[1]['accuracy'], 2) for res in sorted(results.items())]
        df = pd.DataFrame(squareform(accuracy_list), index=diagnose_list, columns=diagnose_list)
        print df
        s += df.to_csv()

        s += '\nSensitivity\n'
        sensitivity_list = [round(res[1]['sensitivity'], 2) for res in sorted(results.items())]
        df = pd.DataFrame(squareform(sensitivity_list), index=diagnose_list, columns=diagnose_list)
        print df
        s += df.to_csv()

        s += '\nSpecificity\n'
        specificity_list = [round(res[1]['specificity'], 2) for res in sorted(results.items())]
        df = pd.DataFrame(squareform(specificity_list), index=diagnose_list, columns=diagnose_list)
        print df
        s += df.to_csv()

        output.write(s)
