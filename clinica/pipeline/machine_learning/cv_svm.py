from multiprocessing.pool import ThreadPool
import numpy as np
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold
from scipy.stats import mode
from clinica.pipeline.machine_learning.svm_utils import evaluate_prediction, calculate_auc


def launch_svc(kernel_train, x_test, y_train, y_test, c, balanced=False, shared_x=None, train_indices=None, test_indices=None):
    if balanced:
        svc = SVC(C=c, kernel='precomputed', tol=1e-6, class_weight='balanced')
    else:
        svc = SVC(C=c, kernel='precomputed', tol=1e-6)

    svc.fit(kernel_train, y_train)
    y_hat = svc.predict(x_test)
    auc = None

    if shared_x is not None and train_indices is not None and test_indices is not None:
        auc = calculate_auc(svc, shared_x, train_indices, test_indices, y_test)

    return y_hat, auc


def inner_grid_search(kernel_train, x_test, y_train, y_test, c, balanced=False):

    y_hat, _ = launch_svc(kernel_train, x_test, y_train, y_test, c, balanced)
    res = evaluate_prediction(y_test, y_hat)

    return res['balanced_accuracy']


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


def outer_cross_validation(kernel_train, shared_x, train_indices, test_indices, x_test, y_train, y_test, async_res, balanced=False):

    best_c, best_acc = select_best_c(async_res)
    y_hat, auc = launch_svc(kernel_train, x_test, y_train, y_test, best_c, balanced, shared_x, train_indices, test_indices)

    result = dict()
    result['best_c'] = best_c
    result['best_acc'] = best_acc
    result['evaluation'] = evaluate_prediction(y_test, y_hat)
    result['y_hat'] = y_hat
    result['auc'] = auc

    return result


def cv_svm(gram_matrix, shared_x, indices, y, c_range, balanced=False, outer_folds=10, inner_folds=10, n_threads=15):

    y_hat = np.zeros(len(y))

    async_result = {}
    for i in range(outer_folds):
        async_result[i] = {}
        for c in c_range:
            async_result[i][c] = {}

    inner_pool = ThreadPool(n_threads)

    skf = StratifiedKFold(n_splits=outer_folds, shuffle=True)
    outer_cv = list(skf.split(np.zeros(len(y)), y))

    #print 'Launching inner threads'
    for i in range(len(outer_cv)):
        train_index, test_index = outer_cv[i]

        outer_gram_matrix = gram_matrix[train_index, :][:, train_index]
        y_train = y[train_index]

        skf = StratifiedKFold(n_splits=inner_folds, shuffle=True)
        inner_cv = list(skf.split(np.zeros(len(y_train)), y_train))

        for j in range(len(inner_cv)):
            inner_train_index, inner_test_index = inner_cv[j]

            inner_gram_matrix = outer_gram_matrix[inner_train_index, :][:, inner_train_index]
            x_test_inner = outer_gram_matrix[inner_test_index, :][:, inner_train_index]
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
        train_indices, test_indices = indices[train_index], indices[test_index]

        outer_async_result[i] = outer_pool.apply_async(outer_cross_validation, (outer_gram_matrix, shared_x, train_indices, test_indices, x_test, y_train, y_test, async_result[i], balanced))

    #print 'All outer threads launched'
    outer_pool.close()
    outer_pool.join()
    #print 'All outer threads finished'

    keys = outer_async_result.keys()
    keys.sort()

    best_c_list = []
    auc_list = []
    for i in keys:
        train_index, test_index = outer_cv[i]
        res = outer_async_result[i].get()
        y_hat[test_index] = res['y_hat']
        best_c_list.append(res['best_c'])
        auc_list.append(res['auc'])

    # Mode of best c for each iteration is selected
    best_c = mode(best_c_list)[0][0]
    # Mean AUC
    auc = np.mean(auc_list)

    if balanced:
        svc = SVC(C=best_c, kernel='precomputed', tol=1e-6, class_weight='balanced')
    else:
        svc = SVC(C=best_c, kernel='precomputed', tol=1e-6)

    svc.fit(gram_matrix, y)

    return y_hat, svc.dual_coef_, svc.support_, svc.intercept_, best_c, auc
