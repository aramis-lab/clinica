from multiprocessing.pool import ThreadPool
import numpy as np
from os.path import join, isdir, dirname
from os import makedirs
import errno
from sklearn.svm import SVC
from sklearn.cross_validation import StratifiedKFold
from sklearn.preprocessing import scale
from scipy.stats import mode
from clinica.pipeline.machine_learning.voxel_based_io import get_caps_image_list, load_data, revert_mask, weights_to_nifti, save_subjects_prediction, results_to_csv
from clinica.pipeline.machine_learning.voxel_based_utils import evaluate_prediction, gram_matrix_linear

__author__ = "Jorge Samper Gonzalez"
__copyright__ = "Copyright 2016, The Aramis Lab Team"
__credits__ = ["Jorge Samper Gonzalez"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"

def launch_svc(kernel_train, x_test, y_train, c, balanced=False):
    if balanced:
        svc = SVC(C=c, kernel='precomputed', tol=1e-6, class_weight='balanced')
    else:
        svc = SVC(C=c, kernel='precomputed', tol=1e-6)

    svc.fit(kernel_train, y_train)
    y_hat = svc.predict(x_test)
    return y_hat



def inner_grid_search(kernel_train, x_test, y_train, y_test, c, balanced=False):

    y_hat = launch_svc(kernel_train, x_test, y_train, c, balanced)
    res = evaluate_prediction(y_test, y_hat)

    return res['accuracy']


def select_best_c(async_res):
    best_c = 0
    best_acc = 0

    k = async_res.keys()
    k.sort()

    for c in k:
        mean_acc = np.mean([res.get() for res in async_res[c].values()])
        if mean_acc > best_acc:
            best_acc = mean_acc
            best_c = c

    return best_c, best_acc


def outer_cross_validation(kernel_train, x_test, y_train, y_test, async_res, balanced=False):

    best_c, best_acc = select_best_c(async_res)
    y_hat = launch_svc(kernel_train, x_test, y_train, best_c, balanced)

    result = dict()
    result['best_c'] = best_c
    result['best_acc'] = best_acc
    result['evaluation'] = evaluate_prediction(y_test, y_hat)
    result['y_hat'] = y_hat

    return result


def nested_folds(gram_matrix, x, y, c_range, balanced=False, outer_folds=10, inner_folds=10, n_threads=15):

    y_hat = np.zeros(len(y))

    async_result = {}
    for i in range(outer_folds):
        async_result[i] = {}
        for c in c_range:
            async_result[i][c] = {}

    inner_pool = ThreadPool(n_threads)

    skf = StratifiedKFold(y, n_folds=outer_folds, shuffle=True)
    outer_cv = list(skf)

    #print 'Launching inner threads'
    for i in range(len(outer_cv)):
        train_index, test_index = outer_cv[i]

        outer_gram_matrix = gram_matrix[train_index, :][:, train_index]
        y_train = y[train_index]

        skf = StratifiedKFold(y_train, n_folds=inner_folds, shuffle=True)
        inner_cv = list(skf)

        for j in range(len(inner_cv)):
            inner_train_index, inner_test_index = inner_cv[j]

            inner_gram_matrix = outer_gram_matrix[inner_train_index,:][:,inner_train_index]
            x_test_inner = outer_gram_matrix[inner_test_index,:][:,inner_train_index]
            y_train_inner, y_test_inner = y_train[inner_train_index], y_train[inner_test_index]

            for c in c_range:
                # print 'Launched %d, %d, %0.5f' %(i, j, c)
                async_result[i][c][j] = inner_pool.apply_async(inner_grid_search, (inner_gram_matrix, x_test_inner, y_train_inner, y_test_inner, c, balanced))

    #print 'All inner threads launched'
    inner_pool.close()
    inner_pool.join()
    #print 'All inner threads finished'

    outer_pool = ThreadPool(n_threads)
    outer_async_result = {}

    for i in range(outer_folds):
        train_index, test_index = outer_cv[i]

        outer_gram_matrix = gram_matrix[train_index, :][:, train_index]
        x_test = gram_matrix[test_index, :][:, train_index]
        y_train, y_test = y[train_index], y[test_index]
        outer_async_result[i] = outer_pool.apply_async(outer_cross_validation, (outer_gram_matrix, x_test, y_train, y_test, async_result[i], balanced))

    #print 'All outer threads launched'
    outer_pool.close()
    outer_pool.join()
    #print 'All outer threads finished'

    keys = outer_async_result.keys()
    keys.sort()

    best_c_list = []
    for i in keys:
        train_index, test_index = outer_cv[i]
        res = outer_async_result[i].get()
        y_hat[test_index] = res['y_hat']
        best_c_list.append(res['best_c'])

    # Mode of best c for each iteration is selected
    best_c = mode(best_c_list)[0][0]
    if balanced:
        svc = SVC(C=best_c, kernel='linear', tol=1e-6, class_weight='balanced')
    else:
        svc = SVC(C=best_c, kernel='linear', tol=1e-6)

    svc.fit(x, y)

    return y_hat, svc.coef_, best_c


def linear_svm_binary_classification(image_list, diagnose_list, output_directory, mask_zeros=True, balanced=False, outer_folds=10, inner_folds=10, n_threads=10, c_range=np.logspace(-6, 2, 17), save_gram_matrix=False, save_subject_classification=False, save_original_weights=False, save_features_image=True):

    results = dict()
    dx_filter = np.unique(diagnose_list)

    print 'Loading ' + str(len(image_list)) + ' subjects'

    x0, orig_shape, data_mask = load_data(image_list, mask=mask_zeros)

    print 'Subjects loaded'
    print 'Calculating Gram matrix'

    x_all = scale(np.nan_to_num(x0))
    gram_matrix = gram_matrix_linear(x_all)

    print 'Gram matrix calculated'

    if save_gram_matrix:
        np.savetxt(join(output_directory, 'gram_matrix.txt'), gram_matrix)

    # Allow loading precalculated gram_matrix?
    # gram_matrix = np.loadtxt(input_gram_matrix)

    for i in range(len(dx_filter)):
        for j in range(i + 1, len(dx_filter)):
            dx1 = dx_filter[i]
            dx2 = dx_filter[j]

            ind1 = []
            ind2 = []
            for k in range(len(diagnose_list)):
                if diagnose_list[k] == dx1:
                    ind1.append(k)
                if diagnose_list[k] == dx2:
                    ind2.append(k)

            indices = ind1 + ind2

            current_subjects = [image_list[k] for k in indices]
            current_diagnosis = [diagnose_list[k] for k in indices]

            x = [x_all[k] for k in indices]
            y = np.array([0] * len(ind1) + [1] * len(ind2))
            gm = gram_matrix[indices, :][:, indices]

            classification_str = dx1 + '_vs_' + dx2 + ('_balanced' if balanced else '_not_balanced')
            print 'Running ' + dx1 + ' vs ' + dx2 + ' classification'

            y_hat, w, c = nested_folds(gm, x, y, c_range, balanced=balanced, outer_folds=outer_folds, inner_folds=inner_folds, n_threads=n_threads)

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

            weights_orig = revert_mask(w, data_mask, orig_shape)

            if save_original_weights:
                np.save(join(output_directory, classification_str + '__weights'), weights_orig)

            if save_features_image:
                weights_to_nifti(weights_orig, image_list[0], join(output_directory, classification_str + '__features_image.nii'))

            if save_subject_classification:
                save_subjects_prediction(current_subjects, current_diagnosis, y, y_hat, join(output_directory, classification_str + '__subjects.csv'))

            results[(dx1, dx2)] = evaluate_prediction(y, y_hat)

    results_to_csv(results, dx_filter, join(output_directory, 'resume.csv'))


def linear_svm_binary_classification_caps(caps_directory,
                                          subjects_visits_tsv,
                                          group_id,
                                          diagnoses_list,
                                          prefix='smwc1',
                                          mask_zeros=True,
                                          balanced=False,
                                          outer_folds=10,
                                          inner_folds=10,
                                          n_threads=10,
                                          c_range=np.logspace(-6, 2, 17),
                                          save_gram_matrix=False,
                                          save_subject_classification=False,
                                          save_original_weights=False,
                                          save_features_image=True):

    output_directory = join(caps_directory, 'group-' + group_id + '/machine_learning/voxel_based_svm/')
    try:
        makedirs(output_directory)

    except OSError as exc:
        if exc.errno == errno.EEXIST and isdir(dirname(output_directory)):
            pass
        else:
            raise

    image_list = get_caps_image_list(caps_directory, subjects_visits_tsv, group_id, prefix)

    linear_svm_binary_classification(image_list, diagnoses_list, output_directory, mask_zeros, balanced, outer_folds,
                                     inner_folds, n_threads, c_range, save_gram_matrix, save_subject_classification,
                                     save_original_weights, save_features_image)
