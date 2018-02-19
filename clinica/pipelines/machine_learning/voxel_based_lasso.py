# coding: utf8

from multiprocessing.pool import ThreadPool
import numpy as np
from os.path import join, isdir, dirname
from os import makedirs
import errno
from sklearn.linear_model import Lasso
from sklearn.cross_validation import StratifiedKFold
from sklearn.preprocessing import scale
from scipy.stats import mode
from clinica.pipelines.machine_learning.voxel_based_io import get_caps_t1_list_OLD, load_data, revert_mask, weights_to_nifti, save_subjects_prediction, results_to_csv
from clinica.pipelines.machine_learning.voxel_based_utils import evaluate_prediction, gram_matrix_linear
# from multiprocessing.managers import BaseManager
import sharedmem
import gc

__author__ = "Jorge Samper Gonzalez"
__copyright__ = "Copyright 2016-2018, The Aramis Lab Team"
__credits__ = ["Jorge Samper Gonzalez"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


def launch_lasso(x_train, x_test, y_train, alpha, positive=False):

    lasso = Lasso(alpha=alpha, fit_intercept=True, precompute=False, positive=positive)
    lasso.fit(x_train, y_train)
    y_hat = lasso.predict(x_test)
    return y_hat


def inner_grid_search(shared_x, indices, train_index, inner_train_index, inner_test_index, y_train, y_test, alpha, positive=False):

    x_train = shared_x[indices[train_index[inner_train_index]], :]
    x_test = shared_x[indices[train_index[inner_test_index]], :]

    y_hat = launch_lasso(x_train, x_test, y_train, alpha, positive)
    res = evaluate_prediction(y_test, y_hat)
    return res['accuracy']


def select_best_alpha(async_res):
    best_alpha = 0
    best_acc = 0

    k = async_res.keys()
    k.sort()

    for alpha in k:
        mean_acc = np.mean([res.get() for res in async_res[alpha].values()])
        if mean_acc > best_acc:
            best_acc = mean_acc
            best_alpha = alpha

    return best_alpha, best_acc


def outer_cross_validation(shared_x, indices, train_index, test_index, y_train, y_test, async_res, positive=False):

    best_alpha, best_acc = select_best_alpha(async_res)

    x_train = shared_x[indices[train_index], :]
    x_test = shared_x[indices[test_index], :]

    y_hat = launch_lasso(x_train, x_test, y_train, best_alpha, positive)

    result = dict()
    result['best_alpha'] = best_alpha
    result['best_acc'] = best_acc
    result['evaluation'] = evaluate_prediction(y_test, y_hat)
    result['y_hat'] = y_hat

    return result


def nested_folds(shared_x, indices, y, alphas, positive=False, outer_folds=10, inner_folds=10, n_threads=15):

    y_hat = np.zeros(len(y))

    async_result = {}
    for i in range(outer_folds):
        async_result[i] = {}
        for alpha in alphas:
            async_result[i][alpha] = {}

    inner_pool = ThreadPool(n_threads)

    skf = StratifiedKFold(y, n_folds=outer_folds, shuffle=True)
    outer_cv = list(skf)

    #print 'Launching inner threads'
    for i in range(len(outer_cv)):
        train_index, test_index = outer_cv[i]

        # outer_gram_matrix = gram_matrix[train_index, :][:, train_index]
        # x_train = x[train_index, :]
        y_train = y[train_index]

        skf = StratifiedKFold(y_train, n_folds=inner_folds, shuffle=True)
        inner_cv = list(skf)

        for j in range(len(inner_cv)):
            inner_train_index, inner_test_index = inner_cv[j]

            # inner_gram_matrix = outer_gram_matrix[inner_train_index, :][:, inner_train_index]

            # x_train_inner, x_test_inner = x_train[inner_train_index, :], x_train[inner_test_index, :]
            y_train_inner, y_test_inner = y_train[inner_train_index], y_train[inner_test_index]

            for alpha in alphas:
                # print 'Launched %d, %d, %0.5f' %(i, j, c)
                async_result[i][alpha][j] = inner_pool.apply_async(inner_grid_search, (shared_x, indices, train_index, inner_train_index, inner_test_index, y_train_inner, y_test_inner, alpha, positive))

    #print 'All inner threads launched'
    inner_pool.close()
    inner_pool.join()
    #print 'All inner threads finished'

    outer_pool = ThreadPool(n_threads)
    outer_async_result = {}

    for i in range(outer_folds):
        train_index, test_index = outer_cv[i]

        # outer_gram_matrix = gram_matrix[train_index, :][:, train_index]
        # x_test = gram_matrix[test_index, :][:, train_index]
        y_train, y_test = y[train_index], y[test_index]
        outer_async_result[i] = outer_pool.apply_async(outer_cross_validation, (shared_x, indices, train_index, test_index, y_train, y_test, async_result[i], positive))

    #print 'All outer threads launched'
    outer_pool.close()
    outer_pool.join()
    #print 'All outer threads finished'

    keys = outer_async_result.keys()
    keys.sort()

    best_alpha_list = []
    for i in keys:
        train_index, test_index = outer_cv[i]
        res = outer_async_result[i].get()
        y_hat[test_index] = res['y_hat']
        best_alpha_list.append(res['best_alpha'])

    # Mode of best alpha for each iteration is selected
    best_alpha = mode(best_alpha_list)[0][0]

    lasso = Lasso(alpha=best_alpha, fit_intercept=True, precompute=False, positive=positive)
    lasso.fit(shared_x[indices, :], y)

    return y_hat, lasso.coef_, lasso.intercept_, best_alpha


def lasso_binary_classification(image_list, diagnosis_list, output_directory,
                                existing_gram_matrix=None, mask_zeros=True, scale_data=False, positive=False,
                                outer_folds=10, inner_folds=10, n_threads=10, alphas=np.arange(0.1, 1.1, 0.1),
                                save_gram_matrix=False, save_subject_classification=False,
                                save_weights=True, save_features_image=True):

    results = dict()
    dx_filter = np.unique(diagnosis_list)

    print 'Loading ' + str(len(image_list)) + ' subjects'
    x0, orig_shape, data_mask = load_data(image_list, mask=mask_zeros)
    print 'Subjects loaded'

    if scale_data:
        x_all = scale(x0)
    else:
        x_all = x0

    # if existing_gram_matrix is not None:
    #     gram_matrix = existing_gram_matrix
    #     if (gram_matrix.shape[0] != gram_matrix.shape[1]) | (gram_matrix.shape[0] != len(image_list)):
    #         raise ValueError('The existing Gram matrix must be a square matrix with number of rows and columns equal to the number of images.')
    # else:
    #     print 'Calculating Gram matrix'
    #     gram_matrix = gram_matrix_linear(x_all)
    #     print 'Gram matrix calculated'

    # if save_gram_matrix:
    #     np.savetxt(join(output_directory, 'gram_matrix.txt'), gram_matrix)

    # BaseManager.register('ndarray', type(x_all))
    # manager = BaseManager()
    # manager.start()
    # shared_x = manager.ndarray(x_all.shape, buffer=x_all)

    shared_x = sharedmem.copy(x_all)
    x_all = None
    gc.collect()

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
            indices = np.array(indices)

            current_subjects = [image_list[k] for k in indices]
            current_diagnosis = [diagnosis_list[k] for k in indices]

            # x = x_all[indices, :]
            y = np.array([0] * len(ind1) + [1] * len(ind2))
            # gm = gram_matrix[indices, :][:, indices]

            classification_str = dx1 + '_vs_' + dx2 + ('_positive' if positive else '')
            print 'Running ' + dx1 + ' vs ' + dx2 + ' classification'

            y_hat, coefficients, intersect, alpha = nested_folds(shared_x, indices, y, alphas, positive=positive, outer_folds=outer_folds, inner_folds=inner_folds, n_threads=n_threads)
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

            if save_weights or save_features_image:
                weights_orig = revert_mask(coefficients, data_mask, orig_shape)

            if save_weights:
                np.save(join(output_directory, classification_str + '__intersect'), intersect)
                np.save(join(output_directory, classification_str + '__weights'), weights_orig)

            if save_features_image:
                weights_to_nifti(weights_orig, image_list[0], join(output_directory, classification_str + '__features_image.nii'))

            if save_subject_classification:
                save_subjects_prediction(current_subjects, current_diagnosis, y, y_hat, join(output_directory, classification_str + '__subjects.csv'))

            results[(dx1, dx2)] = evaluate_prediction(y, y_hat)

    results_to_csv(results, dx_filter, join(output_directory, 'resume' + ('_positive' if positive else '') + '.csv'))


def lasso_binary_classification_caps(caps_directory,
                                     subjects_visits_tsv,
                                     group_id,
                                     diagnosis_list,
                                     prefix,
                                     tissue,
                                     mask_zeros=True,
                                     scale_data=False,
                                     positive=False,
                                     outer_folds=10,
                                     inner_folds=10,
                                     n_threads=10,
                                     alphas=np.arange(0.05, 1.05,0.05),
                                     save_gram_matrix=False,
                                     save_subject_classification=False,
                                     save_weights=True,
                                     save_features_image=True):

    output_directory = join(caps_directory, '/group-' + group_id + '/machine_learning/voxel_based_lasso/')
    try:
        makedirs(output_directory)

    except OSError as exc:
        if exc.errno == errno.EEXIST and isdir(dirname(output_directory)):
            pass
        else:
            raise

    image_list = get_caps_t1_list_OLD(caps_directory, subjects_visits_tsv, group_id, prefix, tissue)

    lasso_binary_classification(image_list, diagnosis_list, output_directory, mask_zeros=mask_zeros,
                                scale_data=scale_data, positive=positive, outer_folds=outer_folds,
                                inner_folds=inner_folds, n_threads=n_threads, alphas=alphas,
                                save_gram_matrix=save_gram_matrix,
                                save_subject_classification=save_subject_classification,
                                save_weights=save_weights,
                                save_features_image=save_features_image)
