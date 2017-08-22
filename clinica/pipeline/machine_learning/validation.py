
from clinica.pipeline.machine_learning import base

from multiprocessing.pool import ThreadPool
import numpy as np
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold


class KFoldCV(base.MLValidation):

    def __init__(self, ml_algorithm):
        self._ml_algorithm = ml_algorithm
        self._best_parameters = []

    def outer_cross_validation(kernel_train, shared_x, train_indices, test_indices, x_test, y_train, y_test,
                               best_c, best_acc, balanced=False):

        y_hat, auc = launch_svc(kernel_train, x_test, y_train, y_test, best_c, balanced, shared_x, train_indices,
                                test_indices)

        result = dict()
        result['best_c'] = best_c
        result['best_acc'] = best_acc
        result['evaluation'] = evaluate_prediction(y_test, y_hat)
        result['y_hat'] = y_hat
        result['auc'] = auc

        return result

    def cv_svm(self, gram_matrix, shared_x, indices, y, balanced=False, outer_folds=10, inner_folds=10,
               n_threads=15):

        y_hat = np.zeros(len(y))

        skf = StratifiedKFold(n_splits=outer_folds, shuffle=True)
        outer_cv = list(skf.split(np.zeros(len(y)), y))
        outer_pool = ThreadPool(n_threads)
        outer_async_result = {}

        for i in range(outer_folds):

            train_index, test_index = outer_cv[i]

            self._best_parameters.append(self._ml_algorithm.parameter_estimation(train_index))

            outer_kernel = gram_matrix[train_index, :][:, train_index]
            x_test = gram_matrix[test_index, :][:, train_index]
            y_train, y_test = y[train_index], y[test_index]
            train_indices, test_indices = indices[train_index], indices[test_index]

            outer_async_result[i] = outer_pool.apply_async(outer_cross_validation, (
            outer_kernel, shared_x, train_indices, test_indices, x_test, y_train, y_test, async_result[i],
            balanced))

        # print 'All outer threads launched'
        outer_pool.close()
        outer_pool.join()
        # print 'All outer threads finished'

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

        # # Mean of best c of each fold is selected
        # best_c = mode(best_c_list)[0][0]
        best_c = np.power(10, np.mean(np.log10(best_c_list)))

        # Mean AUC
        auc = np.mean(auc_list)

        if balanced:
            svc = SVC(C=best_c, kernel='precomputed', tol=1e-6, class_weight='balanced')
        else:
            svc = SVC(C=best_c, kernel='precomputed', tol=1e-6)

        svc.fit(gram_matrix, y)

        return y_hat, svc.dual_coef_, svc.support_, svc.intercept_, best_c, auc
