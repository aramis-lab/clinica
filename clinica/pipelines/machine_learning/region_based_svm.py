# coding: utf8


from clinica.pipelines.machine_learning.region_based_io import get_caps_t1_list, get_caps_pet_list,load_data, features_weights, weights_to_nifti
from clinica.pipelines.machine_learning.svm_utils import evaluate_prediction, gram_matrix_linear, save_subjects_prediction, results_to_tsv
import numpy as np
from sklearn.preprocessing import scale
from os.path import join
from cv_svm import cv_svm
import sharedmem

__author__ = "Jorge Samper Gonzalez"
__copyright__ = "Copyright 2016-2018, The Aramis Lab Team"
__credits__ = ["Jorge Samper Gonzalez", "Simona Bottani"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


def svm_binary_classification(input_image_atlas,
                              subjects_visits_tsv,
                              image_list, diagnosis_list,
                              output_directory,
                              kernel_function=None,
                              existing_gram_matrix=None,
                              mask_zeros=True,
                              scale_data=False,
                              balanced=False,
                              outer_folds=10,
                              inner_folds=10,
                              n_threads=10,
                              c_range=np.logspace(-10, 2, 1000),
                              save_gram_matrix=False,
                              save_subject_classification=False,
                              save_dual_coefficients=False,
                              scaler=None,
                              data_mask=None,
                              save_original_weights=False,
                              save_features_image=True):

    if (kernel_function is None and existing_gram_matrix is None) | (
            kernel_function is not None and existing_gram_matrix is not None):
        raise ValueError('Kernel_function and existing_gram_matrix are mutually exclusive parameters.')

    results = dict()
    dx_filter = np.unique(diagnosis_list)

    print 'Loading ' + str(len(image_list)) + ' subjects'
    x0 = load_data(image_list,subjects_visits_tsv)
    print 'Subjects loaded'
    if scale_data:
        x_all = scale(x0)
    else:
        x_all = x0

    if existing_gram_matrix is None:
        if kernel_function is not None:
            print 'Calculating Gram matrix'
            gram_matrix = kernel_function(x_all)
            print 'Gram matrix calculated'
        else:
            raise ValueError(
                'If a Gram matrix is not provided a function to calculate it (kernel_function) is a required input.')
    else:
        gram_matrix = existing_gram_matrix
        if (gram_matrix.shape[0] != gram_matrix.shape[1]) | (gram_matrix.shape[0] != len(image_list)):
            raise ValueError(
                'The existing Gram matrix must be a square matrix with number of rows and columns equal to the number of images.')

    if save_gram_matrix:
        np.savetxt(join(output_directory, 'gram_matrix.txt'), gram_matrix)

    shared_x = sharedmem.copy(x_all)
    x_all = None
    gc.collect()

    for i in range(len(dx_filter)):
        for j in range(i + 1, len(dx_filter)):
            print j
            dx1 = dx_filter[i]
            dx2 = dx_filter[j]

            ind1 = []
            ind2 = []
            for k in range(len(diagnosis_list)):
                if diagnosis_list[k] == dx1:
                    ind1.append(k)
                if diagnosis_list[k] == dx2:
                    ind2.append(k)

            indices = ind1 + ind2

            current_subjects = [image_list[k] for k in indices]
            current_diagnosis = [diagnosis_list[k] for k in indices]

            y = np.array([0] * len(ind1) + [1] * len(ind2))
            gm = gram_matrix[indices, :][:, indices]

            classification_str = dx1 + '_vs_' + dx2 + ('_balanced' if balanced else '_not_balanced')
            print 'Running ' + dx1 + ' vs ' + dx2 + ' classification'

            y_hat, dual_coefficients, sv_indices, intersect, c, auc = cv_svm(gm, shared_x, np.array(indices), y,
                                                                             c_range, balanced=balanced,
                                                                             outer_folds=outer_folds,
                                                                             inner_folds=inner_folds,
                                                                             n_threads=n_threads)

            evaluation = evaluate_prediction(y, y_hat)
            evaluation['auc'] = auc

            print '\nTrue positive %0.2f' % len(evaluation['predictions'][0])
            print 'True negative %0.2f' % len(evaluation['predictions'][1])
            print 'False positive %0.2f' % len(evaluation['predictions'][2])
            print 'False negative %0.2f' % len(evaluation['predictions'][3])

            print 'AUC %0.2f' % auc
            print 'Accuracy %0.2f' % evaluation['accuracy']
            print 'Balanced accuracy %0.2f' % evaluation['balanced_accuracy']
            print 'Sensitivity %0.2f' % evaluation['sensitivity']
            print 'Specificity %0.2f' % evaluation['specificity']
            print 'Positive predictive value %0.2f' % evaluation['ppv']
            print 'Negative predictive value %0.2f \n' % evaluation['npv']

            if save_dual_coefficients:
                np.save(join(output_directory, classification_str + '__dual_coefficients'), dual_coefficients[0])
                np.save(join(output_directory, classification_str + '__sv_indices'), sv_indices)
                np.save(join(output_directory, classification_str + '__intersect'), intersect)

            if save_original_weights or save_features_image:
                weights_orig = features_weights(current_subjects, dual_coefficients[0], sv_indices, scaler, data_mask)

            if save_original_weights:
                np.save(join(output_directory, classification_str + '__weights'), weights_orig)

            if save_features_image:
                output_image=weights_to_nifti(input_image_atlas,weights_orig)
                output_image.to_filename(join(output_directory,classification_str+'__weights.nii'))

            if save_subject_classification:
                save_subjects_prediction(current_subjects, current_diagnosis, y, y_hat,
                                             join(output_directory, classification_str + '__subjects.tsv'))

            results[(dx1, dx2)] = evaluation  # evaluate_prediction(y, y_hat)

    results_to_tsv(results, dx_filter,
                   join(output_directory, 'resume' + ('_balanced' if balanced else '_not_balanced') + '.tsv'))
    shared_x = None
    gc.collect()
