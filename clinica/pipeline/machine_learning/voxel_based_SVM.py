from multiprocessing.pool import ThreadPool
import numpy as np
from voxel_based_utils import evaluate_prediction
from sklearn.svm import SVC
from sklearn.cross_validation import StratifiedKFold
#
#
# def class_weights(y):
#     ones = sum(y)
#     zeros = len(y) - sum(y)
#
#     min_class = 0 if ones > zeros else 1
#     proportion = (float(ones) / float(zeros)) if ones > zeros else (float(zeros) / float(ones))
#
#     cw = {int(min_class): proportion, int(((min_class - 1) * -1)): 1.0}
#     # print cw
#
#     return cw


def launch_svc(kernel_train, x_test, y_train, y_test, c, balanced=False):

    if balanced:
        svc = SVC(C=c, kernel='precomputed', tol=1e-6, class_weight='balanced')
    else:
        svc = SVC(C=c, kernel='precomputed', tol=1e-6)

    svc.fit(kernel_train, y_train)
    y_hat = svc.predict(x_test)

    res = evaluate_prediction(y_test, y_hat)

    return res['accuracy']


def select_best(kernel_train, x_test, y_train, y_test, async_res, balanced=False):
    best_c = 0
    best_acc = 0

    k = async_res.keys()
    k.sort()

    for c in k:
        mean_acc = np.mean([res.get() for res in async_res[c].values()])
        if mean_acc > best_acc:
            best_acc = mean_acc
            best_c = c

    if balanced:
        svc = SVC(C=best_c, kernel='precomputed', tol=1e-6, class_weight='balanced')
    else:
        svc = SVC(C=best_c, kernel='precomputed', tol=1e-6)

    svc.fit(kernel_train, y_train)
    y_hat = svc.predict(x_test)

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
                async_result[i][c][j] = inner_pool.apply_async(launch_svc, (inner_gram_matrix, x_test_inner, y_train_inner, y_test_inner, c, balanced))

    #print 'All inner threads launched'
    inner_pool.close()
    inner_pool.join()
    #print 'All inner threads finished'

    outer_pool = ThreadPool(n_threads)
    outer_async_result = {}

    for i in range(outer_folds):
        train_index, test_index = outer_cv[i]

        outer_gram_matrix = gram_matrix[train_index,:][:,train_index]
        x_test = gram_matrix[test_index,:][:,train_index]
        y_train, y_test = y[train_index], y[test_index]
        outer_async_result[i] = outer_pool.apply_async(select_best, (outer_gram_matrix, x_test, y_train, y_test, async_result[i], balanced))

    #print 'All outer threads launched'
    outer_pool.close()
    outer_pool.join()
    #print 'All outer threads finished'

    keys = outer_async_result.keys()
    keys.sort()
    for i in keys:
        train_index, test_index = outer_cv[i]
        res = outer_async_result[i].get()
        y_hat[test_index] = res['y_hat']

        print '\nPartition: %d' %i
        print 'Best c: ' + str(res['best_c'])
        print 'Best accuracy for c: ' + str(res['best_acc'])
        #print res['evaluation']

    #print '\nPredicted labels: '
    #print y_hat

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
    print 'Negative predictive value %0.2f' % evaluation['npv']
    print

    if balanced:
        svc = SVC(C=c, kernel='linear', tol=1e-6, class_weight='balanced')
    else:
        svc = SVC(C=c, kernel='linear', tol=1e-6)

    svc.fit(x, y)

    return y_hat, svc.coef_


def SVM_binary_classification(subject_list, diagnose_list, data_directory):

    results = dict()

    dx_filter = ['Depr', 'EOAD', 'LOAD', 'DCB', 'DCL', 'DFT', 'APPl', 'APPs']
    subjects, labels = data_morin_all.get_all_subj_dx(diagnose_file, data_directory + smooth, dx_filter)

    print str(len(subjects)) + ' Subjects'
    # print data_directory + smooth

    x0, data_mask, orig_shape = load_data(subjects, data_directory + smooth, True)

    print '\nX0.shape: '
    print x0.shape
    print

    x_all = preprocessing.scale(np.nan_to_num(x0))
    gram_matrix = gram_matrix_linear(x_all)
    np.savetxt(output_directory + 'gram_matrix_smooth_' + smooth + '.txt', gram_matrix)

    # gram_matrix = np.loadtxt(output_directory + 'gram_matrix_smooth_' + smooth + '.txt')

    c_range = np.logspace(-6, 2, 17)

    for i in range(len(dx_filter)):
        for j in range(i + 1, len(dx_filter)):
            dx1 = dx_filter[i]
            dx2 = dx_filter[j]

            ind1 = labels[dx1]
            ind2 = labels[dx2]

            indices = ind1 + ind2

            subj = [subjects[k] for k in indices]
            x = [x_all[k] for k in indices]
            y = np.array([0] * len(ind1) + [1] * len(ind2))
            gm = gram_matrix[indices, :][:, indices]

            # Balanced
            class_str = dx1 + '_vs_' + dx2 + '_smooth_' + smooth + '_balanced'
            print class_str
            y_hat, w = nested_folds(gm, x, y, c_range, balanced=True, n_threads=25)

            weights_orig = revert_mask(w, data_mask, orig_shape)
            np.save(output_directory + class_str + '__weights', weights_orig)
            weights_to_nifti(weights_orig, template, output_directory + class_str + '__features_image.nii')
            data_morin_all.save_subjects(subj, y, y_hat, diagnose_file, output_directory + class_str + '__subjects.csv')

            results[class_str] = evaluate_prediction(y, y_hat)

            # Not Balanced
            class_str = dx1 + '_vs_' + dx2 + '_smooth_' + smooth + '_not_balanced'
            print class_str
            y_hat, w = nested_folds(gm, x, y, c_range, balanced=False, n_threads=25)

            weights_orig = revert_mask(w, data_mask, orig_shape)
            np.save(output_directory + class_str + '__weights', weights_orig)
            weights_to_nifti(weights_orig, template, output_directory + class_str + '__features_image.nii')
            data_morin_all.save_subjects(subj, y, y_hat, diagnose_file, output_directory + class_str + '__subjects.csv')

            results[class_str] = evaluate_prediction(y, y_hat)

            # print str(datetime.now())

    data_morin_all.save_all_results(results, output_directory + 'resume.csv')

    print str(datetime.now())


