
from clinica.pipeline.machine_learning import base
from multiprocessing.pool import ThreadPool
import numpy as np
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold
from scipy.stats import mode
from clinica.pipeline.machine_learning.svm_utils import evaluate_prediction, calculate_auc


class KernelSVMAlgorithm(base.MLAlgorithm):

    def __init__(self, kernel, y, balanced=False, grid_search_folds=10, c_range=np.logspace(-6, 2, 17), n_threads=15):
        self._kernel = kernel
        self._y = y
        self._balanced = balanced
        self._grid_search_folds = grid_search_folds
        self._c_range = c_range
        self._n_threads = n_threads
        self._async_result = {}

    def launch_svc(self, kernel_train, x_test, y_train, y_test, c, shared_x=None, train_indices=None,
                   test_indices=None):

        if self._balanced:
            svc = SVC(C=c, kernel='precomputed', tol=1e-6, class_weight='balanced')
        else:
            svc = SVC(C=c, kernel='precomputed', tol=1e-6)

        svc.fit(kernel_train, y_train)
        y_hat = svc.predict(x_test)
        auc = None

        if shared_x is not None and train_indices is not None and test_indices is not None:
            auc = calculate_auc(svc, shared_x, train_indices, test_indices, y_test)

        return y_hat, auc

    def inner_grid_search(self, kernel_train, x_test, y_train, y_test, c):

        y_hat, _ = self.launch_svc(kernel_train, x_test, y_train, y_test, c)
        res = evaluate_prediction(y_test, y_hat)

        return res['balanced_accuracy']

    def _select_best_parameter(self):

        c_values = []
        accuracies = []
        for fold in self._async_result.keys():
            best_c = 0
            best_acc = 0
            for c, acc in self._async_result[fold].iteritems():
                if acc > best_acc:
                    best_c = c
                    best_acc = acc
            c_values.append(best_c)
            accuracies.append(best_acc)

        best_acc = np.mean(accuracies)
        best_c = np.power(10, np.mean(np.log10(c_values)))

        return {'c': best_c, 'balanced_accuracy': best_acc}

    def parameter_estimation(self, train_index):

        inner_pool = ThreadPool(self._n_threads)

        for i in range(self._grid_search_folds):
            self._async_result[i] = {}

        outer_kernel = self._kernel[train_index, :][:, train_index]
        y_train = self._y[train_index]

        skf = StratifiedKFold(n_splits=self._grid_search_folds, shuffle=True)
        inner_cv = list(skf.split(np.zeros(len(y_train)), y_train))

        for i in range(len(inner_cv)):
            inner_train_index, inner_test_index = inner_cv[j]

            inner_kernel = outer_kernel[inner_train_index, :][:, inner_train_index]
            x_test_inner = outer_kernel[inner_test_index, :][:, inner_train_index]
            y_train_inner, y_test_inner = y_train[inner_train_index], y_train[inner_test_index]

            for c in self._c_range:
                # print 'Launched %d, %d, %0.5f' %(i, j, c)
                self._async_result[i][c] = inner_pool.apply_async(self.inner_grid_search,
                                                                  (inner_kernel, x_test_inner,
                                                                   y_train_inner, y_test_inner, c))

        # print 'All inner threads launched'
        inner_pool.close()
        inner_pool.join()
        # print 'All inner threads finished'
        return self._select_best_parameter()
