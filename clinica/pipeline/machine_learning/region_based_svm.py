from clinica.pipeline.machine_learning.region_based_io import get_caps_t1_list, get_caps_pet_list,load_data, features_weights, weights_to_nifti
from clinica.pipeline.machine_learning.voxel_based_utils import evaluate_prediction, gram_matrix_linear
from multiprocessing.pool import ThreadPool
import numpy as np
from os.path import join, isdir, dirname
from os import makedirs

def svm_binary_classification(input_image_atlas,image_list, diagnosis_list, output_directory, kernel_function=None, existing_gram_matrix=None, mask_zeros=True, scale_data=False, balanced=False, outer_folds=10, inner_folds=10, n_threads=10, c_range=np.logspace(-10, 2, 1000), save_gram_matrix=False, save_subject_classification=False, save_dual_coefficients=False, scaler=None, data_mask=None, save_original_weights=False, save_features_image=True):

    if (kernel_function is None and existing_gram_matrix is None) | (kernel_function is not None and existing_gram_matrix is not None):
        raise ValueError('Kernel_function and existing_gram_matrix are mutually exclusive parameters.')

    results = dict()
    dx_filter = np.unique(diagnosis_list)

    if kernel_function is not None:
        print 'Loading ' + str(len(image_list)) + ' subjects'
        x0, orig_shape, data_mask = load_data(image_list, mask=mask_zeros)
        print 'Subjects loaded'
        print 'Calculating Gram matrix'

        if scale_data:
            x_all = scale(x0)
        else:
            x_all = x0

        gram_matrix = kernel_function(x_all)
        print 'Gram matrix calculated'

    if existing_gram_matrix is not None:
        gram_matrix = existing_gram_matrix
        if (gram_matrix.shape[0] != gram_matrix.shape[1]) | (gram_matrix.shape[0] != len(image_list)):
            raise ValueError('The existing Gram matrix must be a square matrix with number of rows and columns equal to the number of images.')

    if save_gram_matrix:
        np.savetxt(join(output_directory, 'gram_matrix.txt'), gram_matrix)

    for i in range(len(dx_filter)):
        for j in range(i + 1, len(dx_filter)):
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

            y_hat, dual_coefficients, sv_indices, intersect, c = nested_folds(gm, y, c_range, balanced=balanced, outer_folds=outer_folds, inner_folds=inner_folds, n_threads=n_threads)

            evaluation = evaluate_prediction(y, y_hat)

            print '\nTrue positive %0.2f' % len(evaluation['predictions'][0])
            print 'True negative %0.2f' % len(evaluation['predictions'][1])
            print 'False positive %0.2f' % len(evaluation['predictions'][2])
            print 'False negative %0.2f' % len(evaluation['predictions'][3])

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
                weights_orig = features_weights(current_subjects, dual_coefficients[0], sv_indices)

            if save_original_weights:
                np.save(join(output_directory, classification_str + '__weights'), weights_orig)

            if save_features_image:
                output_image=weights_to_nifti(input_image_atlas,weights_orig)
                output_image.to_filename(join(output_directory,classification_str+'__weights.nii'))

            if save_subject_classification:
                save_subjects_prediction(current_subjects, current_diagnosis, y, y_hat, join(output_directory, classification_str + '__subjects.csv'))

            results[(dx1, dx2)] = evaluate_prediction(y, y_hat)

    results_to_csv(results, dx_filter, join(output_directory, 'resume' + ('_balanced' if balanced else '_not_balanced') + '.csv'))