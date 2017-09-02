
from os import path
import json
from multiprocessing.pool import ThreadPool

import numpy as np
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import roc_auc_score

from clinica.pipeline.machine_learning import base
import clinica.pipeline.machine_learning.svm_utils as utils


class DualSVMAlgorithm(base.MLAlgorithm):

    def __init__(self, kernel, y, balanced=True, grid_search_folds=10, c_range=np.logspace(-6, 2, 17), n_threads=15):
        self._kernel = kernel
        self._y = y
        self._balanced = balanced
        self._grid_search_folds = grid_search_folds
        self._c_range = c_range
        self._n_threads = n_threads

    def _launch_svc(self, kernel_train, x_test, y_train, y_test, c):

        if self._balanced:
            svc = SVC(C=c, kernel='precomputed', probability=True, tol=1e-6, class_weight='balanced')
        else:
            svc = SVC(C=c, kernel='precomputed', probability=True, tol=1e-6)

        svc.fit(kernel_train, y_train)
        y_hat = svc.predict(x_test)
        proba_test = svc.predict_proba(x_test)[:, 1]
        auc = roc_auc_score(y_test, proba_test)

        return svc, y_hat, auc

    def _grid_search(self, kernel_train, x_test, y_train, y_test, c):

        _, y_hat, _ = self._launch_svc(kernel_train, x_test, y_train, y_test, c)
        res = utils.evaluate_prediction(y_test, y_hat)

        return res['balanced_accuracy']

    def _select_best_parameter(self, async_result):

        c_values = []
        accuracies = []
        for fold in async_result.keys():
            best_c = -1
            best_acc = -1

            for c, async_acc in async_result[fold].iteritems():

                acc = async_acc.get()
                if acc > best_acc:
                    best_c = c
                    best_acc = acc
            c_values.append(best_c)
            accuracies.append(best_acc)

        best_acc = np.mean(accuracies)
        best_c = np.power(10, np.mean(np.log10(c_values)))

        return {'c': best_c, 'balanced_accuracy': best_acc}

    def evaluate(self, train_index, test_index):

        inner_pool = ThreadPool(self._n_threads)
        async_result = {}
        for i in range(self._grid_search_folds):
            async_result[i] = {}

        outer_kernel = self._kernel[train_index, :][:, train_index]
        y_train = self._y[train_index]

        skf = StratifiedKFold(n_splits=self._grid_search_folds, shuffle=True)
        inner_cv = list(skf.split(np.zeros(len(y_train)), y_train))

        for i in range(len(inner_cv)):
            inner_train_index, inner_test_index = inner_cv[i]

            inner_kernel = outer_kernel[inner_train_index, :][:, inner_train_index]
            x_test_inner = outer_kernel[inner_test_index, :][:, inner_train_index]
            y_train_inner, y_test_inner = y_train[inner_train_index], y_train[inner_test_index]

            for c in self._c_range:
                async_result[i][c] = inner_pool.apply_async(self._grid_search,
                                                            (inner_kernel, x_test_inner,
                                                             y_train_inner, y_test_inner, c))
                #print i, c, async_result[i][c]
        inner_pool.close()
        inner_pool.join()

        best_parameter = self._select_best_parameter(async_result)
        x_test = self._kernel[test_index, :][:, train_index]
        y_train, y_test = self._y[train_index], self._y[test_index]

        _, y_hat, auc = self._launch_svc(outer_kernel, x_test, y_train, y_test, best_parameter['c'])

        result = dict()
        result['best_parameter'] = best_parameter
        result['evaluation'] = utils.evaluate_prediction(y_test, y_hat)
        result['y_hat'] = y_hat
        result['y'] = y_test
        result['y_index'] = test_index
        result['auc'] = auc

        return result

    def apply_best_parameters(self, results_list):

        best_c_list = []
        bal_acc_list = []

        for result in results_list:
            best_c_list.append(result['best_parameter']['c'])
            bal_acc_list.append(result['best_parameter']['balanced_accuracy'])

        # 10^(mean of log10 of best Cs of each fold) is selected
        best_c = np.power(10, np.mean(np.log10(best_c_list)))
        # Mean balanced accuracy
        mean_bal_acc = np.mean(bal_acc_list)

        if self._balanced:
            svc = SVC(C=best_c, kernel='precomputed', probability=True, tol=1e-6, class_weight='balanced')
        else:
            svc = SVC(C=best_c, kernel='precomputed', probability=True, tol=1e-6)

        svc.fit(self._kernel, self._y)

        return svc, {'c': best_c, 'balanced_accuracy': mean_bal_acc}

    def save_classifier(self, classifier, output_dir):

        np.savetxt(path.join(output_dir, 'dual_coefficients.txt'), classifier.dual_coef_)
        np.savetxt(path.join(output_dir, 'support_vectors_indices.txt'), classifier.support_)
        np.savetxt(path.join(output_dir, 'intersect.txt'), classifier.intercept_)

    def save_weights(self, classifier, x, output_dir):

        dual_coefficients = classifier.dual_coef_
        sv_indices = classifier.support_

        weighted_sv = dual_coefficients.transpose() * x[sv_indices]
        weights = np.sum(weighted_sv, 0)

        np.savetxt(path.join(output_dir, 'weights.txt'), weights)

        return weights

    def save_parameters(self, parameters_dict, output_dir):
        with open(path.join(output_dir, 'best_parameters.json'), 'w') as f:
            json.dump(parameters_dict, f)


class LogisticReg(base.MLAlgorithm):
    
    def __init__(self, x, y, penalty='l2', balanced=False, grid_search_folds=10, c_range=np.logspace(-6, 2, 17), n_threads=15):
        """ penalty can either be 'l2' or 'l1'"""
        self._penalty = penalty
        self._x = x
        self._y = y
        self._balanced = balanced
        self._grid_search_folds = grid_search_folds
        self._c_range = c_range
        self._n_threads = n_threads
    
    def _launch_logistic_reg(self, x_train, x_test, y_train, y_test, c, shared_x=None, train_indices=None,
                             test_indices=None):
        
        # x_train_, mean_x, std_x = centered_normalised_data(x_train)
        # x_test_ = (x_test - mean_x)/std_x
        
        if self._balanced:
            classifier = LogisticRegression(penalty=self._penalty, tol=1e-6, C=c, class_weight='balanced')
        else:
            classifier = LogisticRegression(penalty=self._penalty, tol=1e-6, C=c)
        
        classifier.fit(x_train, y_train)
        y_hat = classifier.predict(x_test)
        proba_test = classifier.predict_proba(x_test)[:, 1]
        auc = roc_auc_score(y_test, proba_test)

        return classifier, y_hat, auc

    def _grid_search(self, x_train, x_test, y_train, y_test, c):
        
        _, y_hat, _ = self._launch_logistic_reg(x_train, x_test, y_train, y_test, c)
        res = utils.evaluate_prediction(y_test, y_hat)
        
        return res['balanced_accuracy']

    def _select_best_parameter(self, async_result):
        
        c_values = []
        accuracies = []
        for fold in async_result.keys():
            best_c = -1
            best_acc = -1
            
            for c, async_acc in async_result[fold].iteritems():
                
                acc = async_acc.get()
                if acc > best_acc:
                    best_c = c
                    best_acc = acc
            c_values.append(best_c)
            accuracies.append(best_acc)
        
        best_acc = np.mean(accuracies)
        best_c = np.power(10, np.mean(np.log10(c_values)))
        
        return {'c': best_c, 'balanced_accuracy': best_acc}
    
    def evaluate(self, train_index, test_index):
        
        inner_pool = ThreadPool(self._n_threads)
        async_result = {}
        for i in range(self._grid_search_folds):
            async_result[i] = {}

        x_train = self._x[train_index]
        y_train = self._y[train_index]
        
        skf = StratifiedKFold(n_splits=self._grid_search_folds, shuffle=True)
        inner_cv = list(skf.split(np.zeros(len(y_train)), y_train))
        
        for i in range(len(inner_cv)):
            inner_train_index, inner_test_index = inner_cv[i]
            
            x_train_inner = x_train[inner_train_index]
            x_test_inner = x_train[inner_test_index]
            y_train_inner = y_train[inner_train_index]
            y_test_inner = y_train[inner_test_index]

            for c in self._c_range:
                async_result[i][c] = inner_pool.apply_async(self._grid_search,
                                                            (x_train_inner, x_test_inner,
                                                             y_train_inner, y_test_inner, c))
        inner_pool.close()
        inner_pool.join()
        
        best_parameter = self._select_best_parameter(async_result)
        x_test = self._x[test_index]
        y_test = self._y[test_index]
            
        _, y_hat, auc = self._launch_logistic_reg(x_train, x_test, y_train, y_test, best_parameter['c'])
            
        result = dict()
        result['best_parameter'] = best_parameter
        result['evaluation'] = utils.evaluate_prediction(y_test, y_hat)
        result['y_hat'] = y_hat
        result['y'] = y_test
        result['y_index'] = test_index
        result['auc'] = auc
        
        return result

    def apply_best_parameters(self, results_list):
    
        best_c_list = []
        bal_acc_list = []
        
        for result in results_list:
            best_c_list.append(result['best_parameter']['c'])
            bal_acc_list.append(result['best_parameter']['balanced_accuracy'])

        # 10^(mean of log10 of best Cs of each fold) is selected
        best_c = np.power(10, np.mean(np.log10(best_c_list)))
        # Mean balanced accuracy
        mean_bal_acc = np.mean(bal_acc_list)
        
        if self._balanced:
            classifier = LogisticRegression(C=best_c, penalty=self._penalty, tol=1e-6, class_weight='balanced')
        else:
            classifier = LogisticRegression(C=best_c, penalty=self._penalty, tol=1e-6)
        
        classifier.fit(self._x, self._y)
        
        return classifier, {'c': best_c, 'balanced_accuracy': mean_bal_acc}

    def save_classifier(self, classifier, output_dir):
        
        np.savetxt(path.join(output_dir, 'weights.txt'), classifier.coef_)
        np.savetxt(path.join(output_dir, 'intercept.txt'), classifier.intercept_)

    def save_weights(self, classifier, output_dir):
    
        np.savetxt(path.join(output_dir, 'weights.txt'), classifier.coef_)
        return classifier.coef_
    
    def save_parameters(self, parameters_dict, output_dir):
        
        with open(path.join(output_dir, 'best_parameters.json'), 'w') as f:
            json.dump(parameters_dict, f)

    def _centered_normalised_data(features):
        std = np.std(features, axis=0)
        std[np.where(std == 0)[0]] = 1.
        mean = np.mean(features, axis=0)
        features_bis = (features - mean)/std
        return features_bis, mean, std




from sklearn.ensemble import RandomForestClassifier

class RandomForest(base.MLAlgorithm):
    
    def __init__(self, x, y, balanced=False, grid_search_folds=10, n_estimators_range=range(5, 15, 1), n_threads=15):
        self._x = x
        self._y = y
        self._balanced = balanced
        self._grid_search_folds = grid_search_folds
        self._n_estimators_range = n_estimators_range
        self._n_threads = n_threads
    
    def _launch_random_forest(self, x_train, x_test, y_train, y_test, n_estimators):
        
        if self._balanced:
            classifier = RandomForestClassifier(n_estimators=n_estimators, class_weight='balanced')
        else:
            classifier = RandomForestClassifier(n_estimators=n_estimators)
        
        classifier.fit(x_train, y_train)
        y_hat = classifier.predict(x_test)
        proba_test = classifier.predict_proba(x_test)[:, 1]
        auc = roc_auc_score(y_test, proba_test)
        
        return classifier, y_hat, auc

    def _grid_search(self, x_train, x_test, y_train, y_test, n_estimators):
    
        _, y_hat, _ = self._launch_random_forest(x_train, x_test, y_train, y_test, n_estimators)
        res = utils.evaluate_prediction(y_test, y_hat)
        
        return res['balanced_accuracy']
    
    def _select_best_parameter(self, async_result):
        
        n_estimators_list = []
        accuracies = []
        for fold in async_result.keys():
            best_n_estimators = -1
            best_acc = -1
            
            for n_estimators, async_acc in async_result[fold].iteritems():
                
                acc = async_acc.get()
                if acc > best_acc:
                    best_n_estimators = n_estimators
                    best_acc = acc
            n_estimators_list.append(best_n_estimators)
            accuracies.append(best_acc)
        
        best_acc = np.mean(accuracies)
        best_n_estimators = int(round(np.mean(n_estimators_list)))
        
        return {'n_estimators': best_n_estimators, 'balanced_accuracy': best_acc}
    
    def evaluate(self, train_index, test_index):
        
        inner_pool = ThreadPool(self._n_threads)
        async_result = {}
        for i in range(self._grid_search_folds):
            async_result[i] = {}
        
        x_train = self._x[train_index]
        y_train = self._y[train_index]
        
        skf = StratifiedKFold(n_splits=self._grid_search_folds, shuffle=True)
        inner_cv = list(skf.split(np.zeros(len(y_train)), y_train))
        
        for i in range(len(inner_cv)):
            inner_train_index, inner_test_index = inner_cv[i]
            
            x_train_inner = x_train[inner_train_index]
            x_test_inner = x_train[inner_test_index]
            y_train_inner = y_train[inner_train_index]
            y_test_inner = y_train[inner_test_index]
            
            for n_estimators in self._n_estimators_range:
                async_result[i][n_estimators] = inner_pool.apply_async(self._grid_search,
                                                                       (x_train_inner, x_test_inner,
                                                                        y_train_inner, y_test_inner, n_estimators))
        inner_pool.close()
        inner_pool.join()
        
        best_parameter = self._select_best_parameter(async_result)
        x_test = self._x[test_index]
        y_test = self._y[test_index]
        
        _, y_hat, auc = self._launch_random_forest(x_train, x_test, y_train, y_test, best_parameter['n_estimators'])
        
        result = dict()
        result['best_parameter'] = best_parameter
        result['evaluation'] = utils.evaluate_prediction(y_test, y_hat)
        result['y_hat'] = y_hat
        result['y'] = y_test
        result['y_index'] = test_index
        result['auc'] = auc
        
        return result

    def apply_best_parameters(self, results_list):
    
        best_n_estimators_list = []
        bal_acc_list = []
        
        for result in results_list:
            best_n_estimators_list.append(result['best_parameter']['n_estimators'])
            bal_acc_list.append(result['best_parameter']['balanced_accuracy'])
    
        # best n_estimators is the average of all n_estimators
        best_n_estimators = int(round(np.mean(best_n_estimators_list)))
        # Mean balanced accuracy
        mean_bal_acc = np.mean(bal_acc_list)
        
        if self._balanced:
            classifier = RandomForestClassifier(n_estimators=best_n_estimators, class_weight='balanced')
        else:
            classifier = RandomForestClassifier(n_estimators=best_n_estimators)
        
        classifier.fit(self._x, self._y)
        
        return classifier, {'n_estimators': best_n_estimators, 'balanced_accuracy': mean_bal_acc}

    def save_classifier(self, classifier, output_dir):
        np.savetxt(path.join(output_dir, 'feature_importances.txt'), classifier.feature_importances_)
        #print classifier.estimators_
        #np.savetxt(path.join(output_dir, 'estimators.txt'), str(classifier.estimators_))

    def save_weights(self, classifier, output_dir):

        np.savetxt(path.join(output_dir, 'weights.txt'), classifier.feature_importances_)
        return classifier.coef_
    
    def save_parameters(self, parameters_dict, output_dir):
        
        with open(path.join(output_dir, 'best_parameters.json'), 'w') as f:
            json.dump(parameters_dict, f)
