# coding: utf8

import os
from os import path
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold, StratifiedShuffleSplit
from multiprocessing.pool import ThreadPool

from clinica.pipelines.machine_learning import base

__author__ = "Jorge Samper-Gonzalez"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Jorge Samper-Gonzalez"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper-Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


class KFoldCV(base.MLValidation):

    def __init__(self, ml_algorithm):
        self._ml_algorithm = ml_algorithm
        self._fold_results = []
        self._classifier = None
        self._best_params = None
        self._cv = None

    def validate(self, y, n_folds=10, splits_indices=None, n_threads=15):

        if splits_indices is None:
            skf = StratifiedKFold(n_splits=n_folds, shuffle=True)
            self._cv = list(skf.split(np.zeros(len(y)), y))
        else:
            self._cv = splits_indices

        async_pool = ThreadPool(n_threads)
        async_result = {}

        for i in range(n_folds):

            train_index, test_index = self._cv[i]
            async_result[i] = async_pool.apply_async(self._ml_algorithm.evaluate, (train_index, test_index))

        async_pool.close()
        async_pool.join()

        for i in range(n_folds):
            self._fold_results.append(async_result[i].get())

        self._classifier, self._best_params = self._ml_algorithm.apply_best_parameters(self._fold_results)

        return self._classifier, self._best_params, self._fold_results

    def save_results(self, output_dir):
        if self._fold_results is None:
            raise Exception("No results to save. Method validate() must be run before save_results().")

        subjects_folds = []
        results_folds = []
        container_dir = path.join(output_dir, 'folds')

        if not path.exists(container_dir):
            os.makedirs(container_dir)

        for i in range(len(self._fold_results)):
            subjects_df = pd.DataFrame({'y': self._fold_results[i]['y'],
                                        'y_hat': self._fold_results[i]['y_hat'],
                                        'y_index': self._fold_results[i]['y_index']})
            subjects_df.to_csv(path.join(container_dir, 'subjects_fold-' + str(i) + '.tsv'),
                               index=False, sep='\t', encoding='utf-8')
            subjects_folds.append(subjects_df)

            results_df = pd.DataFrame({'balanced_accuracy': self._fold_results[i]['evaluation']['balanced_accuracy'],
                                       'auc': self._fold_results[i]['auc'],
                                       'accuracy': self._fold_results[i]['evaluation']['accuracy'],
                                       'sensitivity': self._fold_results[i]['evaluation']['sensitivity'],
                                       'specificity': self._fold_results[i]['evaluation']['specificity'],
                                       'ppv': self._fold_results[i]['evaluation']['ppv'],
                                       'npv': self._fold_results[i]['evaluation']['npv'],
                                       'train_balanced_accuracy': self._fold_results[i]['evaluation_train']['balanced_accuracy'],
                                       'train_accuracy': self._fold_results[i]['evaluation_train']['accuracy'],
                                       'train_sensitivity': self._fold_results[i]['evaluation_train']['sensitivity'],
                                       'train_specificity': self._fold_results[i]['evaluation_train']['specificity'],
                                       'train_ppv': self._fold_results[i]['evaluation_train']['ppv'],
                                       'train_npv': self._fold_results[i]['evaluation_train']['npv']
                                       }, index=['i', ])

            results_df.to_csv(path.join(container_dir, 'results_fold-' + str(i) + '.tsv'),
                              index=False, sep='\t', encoding='utf-8')
            results_folds.append(results_df)

        all_subjects = pd.concat(subjects_folds)
        all_subjects.to_csv(path.join(output_dir, 'subjects.tsv'),
                            index=False, sep='\t', encoding='utf-8')

        all_results = pd.concat(results_folds)
        all_results.to_csv(path.join(output_dir, 'results.tsv'),
                           index=False, sep='\t', encoding='utf-8')

        mean_results = pd.DataFrame(all_results.apply(np.nanmean).to_dict(), columns=all_results.columns, index=[0, ])
        mean_results.to_csv(path.join(output_dir, 'mean_results.tsv'),
                            index=False, sep='\t', encoding='utf-8')

        print("Mean results of the classification:")
        print("Balanced accuracy: %s" % (mean_results['balanced_accuracy'].to_string(index=False)))
        print("specificity: %s" % (mean_results['specificity'].to_string(index=False)))
        print("sensitivity: %s" % (mean_results['sensitivity'].to_string(index=False)))
        print("auc: %s" % (mean_results['auc'].to_string(index=False)))


class RepeatedKFoldCV(base.MLValidation):

    def __init__(self, ml_algorithm):
        self._ml_algorithm = ml_algorithm
        self._repeated_fold_results = []
        self._classifier = None
        self._best_params = None
        self._cv = None

    def validate(self, y, n_iterations=100, n_folds=10, n_threads=15):

        async_pool = ThreadPool(n_threads)
        async_result = {}
        self._cv = []

        for r in range(n_iterations):
            skf = StratifiedKFold(n_splits=n_folds, shuffle=True)
            self._cv.append(list(skf.split(np.zeros(len(y)), y)))
            async_result[r] = {}
            self._repeated_fold_results.append([])

            for i in range(n_folds):

                train_index, test_index = self._cv[r][i]
                async_result[r][i] = async_pool.apply_async(self._ml_algorithm.evaluate, (train_index, test_index))

        async_pool.close()
        async_pool.join()
        for r in range(n_iterations):
            for i in range(n_folds):
                self._repeated_fold_results[r].append(async_result[r][i].get())

        # TODO Find a better way to estimate best parameter
        flat_results = [result for fold in self._repeated_fold_results for result in fold]
        self._classifier, self._best_params = self._ml_algorithm.apply_best_parameters(flat_results)

        return self._classifier, self._best_params, self._repeated_fold_results

    def save_results(self, output_dir):
        if self._repeated_fold_results is None:
            raise Exception("No results to save. Method validate() must be run before save_results().")

        all_results_list = []
        all_subjects_list = []

        for iteration in range(len(self._repeated_fold_results)):

            iteration_dir = path.join(output_dir, 'iteration-' + str(iteration))
            if not path.exists(iteration_dir):
                os.makedirs(iteration_dir)

            iteration_subjects_list = []
            iteration_results_list = []
            folds_dir = path.join(iteration_dir, 'folds')

            if not path.exists(folds_dir):
                os.makedirs(folds_dir)

            for i in range(len(self._repeated_fold_results[iteration])):
                subjects_df = pd.DataFrame({'y': self._repeated_fold_results[iteration][i]['y'],
                                            'y_hat': self._repeated_fold_results[iteration][i]['y_hat'],
                                            'y_index': self._repeated_fold_results[iteration][i]['y_index']})
                subjects_df.to_csv(path.join(folds_dir, 'subjects_fold-' + str(i) + '.tsv'),
                                   index=False, sep='\t', encoding='utf-8')
                iteration_subjects_list.append(subjects_df)

                results_df = pd.DataFrame(
                    {'balanced_accuracy': self._repeated_fold_results[iteration][i]['evaluation']['balanced_accuracy'],
                     'auc': self._repeated_fold_results[iteration][i]['auc'],
                     'accuracy': self._repeated_fold_results[iteration][i]['evaluation']['accuracy'],
                     'sensitivity': self._repeated_fold_results[iteration][i]['evaluation']['sensitivity'],
                     'specificity': self._repeated_fold_results[iteration][i]['evaluation']['specificity'],
                     'ppv': self._repeated_fold_results[iteration][i]['evaluation']['ppv'],
                     'npv': self._repeated_fold_results[iteration][i]['evaluation']['npv'],
                     'train_balanced_accuracy': self._repeated_fold_results[iteration][i]['evaluation_train']['balanced_accuracy'],
                     'train_accuracy': self._repeated_fold_results[iteration][i]['evaluation_train']['accuracy'],
                     'train_sensitivity': self._repeated_fold_results[iteration][i]['evaluation_train']['sensitivity'],
                     'train_specificity': self._repeated_fold_results[iteration][i]['evaluation_train']['specificity'],
                     'train_ppv': self._repeated_fold_results[iteration][i]['evaluation_train']['ppv'],
                     'train_npv': self._repeated_fold_results[iteration][i]['evaluation_train']['npv']
                     }, index=['i', ])
                results_df.to_csv(path.join(folds_dir, 'results_fold-' + str(i) + '.tsv'),
                                  index=False, sep='\t', encoding='utf-8')
                iteration_results_list.append(results_df)

            iteration_subjects_df = pd.concat(iteration_subjects_list)
            iteration_subjects_df.to_csv(path.join(iteration_dir, 'subjects.tsv'),
                                         index=False, sep='\t', encoding='utf-8')
            all_subjects_list.append(iteration_subjects_df)

            iteration_results_df = pd.concat(iteration_results_list)
            iteration_results_df.to_csv(path.join(iteration_dir, 'results.tsv'),
                                        index=False, sep='\t', encoding='utf-8')

            mean_results_df = pd.DataFrame(iteration_results_df.apply(np.nanmean).to_dict(),
                                           columns=iteration_results_df.columns, index=[0, ])
            mean_results_df.to_csv(path.join(iteration_dir, 'mean_results.tsv'),
                                   index=False, sep='\t', encoding='utf-8')
            all_results_list.append(mean_results_df)

        all_subjects_df = pd.concat(all_subjects_list)
        all_subjects_df.to_csv(path.join(output_dir, 'subjects.tsv'),
                               index=False, sep='\t', encoding='utf-8')

        all_results_df = pd.concat(all_results_list)
        all_results_df.to_csv(path.join(output_dir, 'results.tsv'),
                              index=False, sep='\t', encoding='utf-8')

        mean_results_df = pd.DataFrame(all_results_df.apply(np.nanmean).to_dict(),
                                       columns=all_results_df.columns, index=[0, ])
        mean_results_df.to_csv(path.join(output_dir, 'mean_results.tsv'),
                               index=False, sep='\t', encoding='utf-8')

        print("Mean results of the classification:")
        print("Balanced accuracy: %s" % (mean_results_df['balanced_accuracy'].to_string(index=False)))
        print("specificity: %s" % (mean_results_df['specificity'].to_string(index=False)))
        print("sensitivity: %s" % (mean_results_df['sensitivity'].to_string(index=False)))
        print("auc: %s" % (mean_results_df['auc'].to_string(index=False)))


class RepeatedHoldOut(base.MLValidation):

    def __init__(self, ml_algorithm, n_iterations=100, test_size=0.3):
        self._ml_algorithm = ml_algorithm
        self._split_results = []
        self._classifier = None
        self._best_params = None
        self._cv = None
        self._n_iterations = n_iterations
        self._test_size = test_size
        self._error_resampled_t = None
        self._error_corrected_resampled_t = None
        self._bal_accuracy_resampled_t = None
        self._bal_accuracy_corrected_resampled_t = None

    def validate(self, y, n_threads=15, splits_indices=None, inner_cv=True):

        if splits_indices is None:
            splits = StratifiedShuffleSplit(n_splits=self._n_iterations, test_size=self._test_size)
            self._cv = list(splits.split(np.zeros(len(y)), y))
        else:
            self._cv = splits_indices
        async_pool = ThreadPool(n_threads)
        async_result = {}

        for i in range(self._n_iterations):

            train_index, test_index = self._cv[i]
            if inner_cv:
                async_result[i] = async_pool.apply_async(self._ml_algorithm.evaluate, (train_index, test_index))
            else:
                async_result[i] = async_pool.apply_async(self._ml_algorithm.evaluate_no_cv, (train_index, test_index))

        async_pool.close()
        async_pool.join()

        for i in range(self._n_iterations):
            self._split_results.append(async_result[i].get())

        self._classifier, self._best_params = self._ml_algorithm.apply_best_parameters(self._split_results)
        return self._classifier, self._best_params, self._split_results

    def save_results(self, output_dir):
        if self._split_results is None:
            raise Exception("No results to save. Method validate() must be run before save_results().")

        all_results_list = []
        all_train_subjects_list = []
        all_test_subjects_list = []

        for iteration in range(len(self._split_results)):

            iteration_dir = path.join(output_dir, 'iteration-' + str(iteration))
            if not path.exists(iteration_dir):
                os.makedirs(iteration_dir)
            iteration_train_subjects_df = pd.DataFrame({'iteration': iteration,
                                                        'y': self._split_results[iteration]['y_train'],
                                                        'y_hat': self._split_results[iteration]['y_hat_train'],
                                                        'subject_index': self._split_results[iteration]['x_index']})
            iteration_train_subjects_df.to_csv(path.join(iteration_dir, 'train_subjects.tsv'),
                                               index=False, sep='\t', encoding='utf-8')
            all_train_subjects_list.append(iteration_train_subjects_df)

            iteration_test_subjects_df = pd.DataFrame({'iteration': iteration,
                                                       'y': self._split_results[iteration]['y'],
                                                       'y_hat': self._split_results[iteration]['y_hat'],
                                                       'subject_index': self._split_results[iteration]['y_index']})
            iteration_test_subjects_df.to_csv(path.join(iteration_dir, 'test_subjects.tsv'),
                                              index=False, sep='\t', encoding='utf-8')
            all_test_subjects_list.append(iteration_test_subjects_df)

            iteration_results_df = pd.DataFrame(
                    {'balanced_accuracy': self._split_results[iteration]['evaluation']['balanced_accuracy'],
                     'auc': self._split_results[iteration]['auc'],
                     'accuracy': self._split_results[iteration]['evaluation']['accuracy'],
                     'sensitivity': self._split_results[iteration]['evaluation']['sensitivity'],
                     'specificity': self._split_results[iteration]['evaluation']['specificity'],
                     'ppv': self._split_results[iteration]['evaluation']['ppv'],
                     'npv': self._split_results[iteration]['evaluation']['npv'],
                     'train_balanced_accuracy': self._split_results[iteration]['evaluation_train']['balanced_accuracy'],
                     'train_accuracy': self._split_results[iteration]['evaluation_train']['accuracy'],
                     'train_sensitivity': self._split_results[iteration]['evaluation_train']['sensitivity'],
                     'train_specificity': self._split_results[iteration]['evaluation_train']['specificity'],
                     'train_ppv': self._split_results[iteration]['evaluation_train']['ppv'],
                     'train_npv': self._split_results[iteration]['evaluation_train']['npv']
                     }, index=['i', ])
            iteration_results_df.to_csv(path.join(iteration_dir, 'results.tsv'),
                                        index=False, sep='\t', encoding='utf-8')

            # mean_results_df = pd.DataFrame(iteration_results_df.apply(np.nanmean).to_dict(),
            #                                columns=iteration_results_df.columns, index=[0, ])
            # mean_results_df.to_csv(path.join(iteration_dir, 'mean_results.tsv'),
            #                        index=False, sep='\t', encoding='utf-8')
            all_results_list.append(iteration_results_df)

        all_train_subjects_df = pd.concat(all_train_subjects_list)
        all_train_subjects_df.to_csv(path.join(output_dir, 'train_subjects.tsv'),
                                     index=False, sep='\t', encoding='utf-8')

        all_test_subjects_df = pd.concat(all_test_subjects_list)
        all_test_subjects_df.to_csv(path.join(output_dir, 'test_subjects.tsv'),
                                    index=False, sep='\t', encoding='utf-8')

        all_results_df = pd.concat(all_results_list)
        all_results_df.to_csv(path.join(output_dir, 'results.tsv'),
                              index=False, sep='\t', encoding='utf-8')

        mean_results_df = pd.DataFrame(all_results_df.apply(np.nanmean).to_dict(),
                                       columns=all_results_df.columns, index=[0, ])
        mean_results_df.to_csv(path.join(output_dir, 'mean_results.tsv'),
                               index=False, sep='\t', encoding='utf-8')

        print("Mean results of the classification:")
        print("Balanced accuracy: %s" % (mean_results_df['balanced_accuracy'].to_string(index=False)))
        print("specificity: %s" % (mean_results_df['specificity'].to_string(index=False)))
        print("sensitivity: %s" % (mean_results_df['sensitivity'].to_string(index=False)))
        print("auc: %s" % (mean_results_df['auc'].to_string(index=False)))

        self.compute_error_variance()
        self.compute_accuracy_variance()

        variance_df = pd.DataFrame({'bal_accuracy_resampled_t': self._bal_accuracy_resampled_t,
                                    'bal_accuracy_corrected_resampled_t': self._bal_accuracy_corrected_resampled_t,
                                    'error_resampled_t': self._error_resampled_t,
                                    'error_corrected_resampled_t': self._error_corrected_resampled_t}, index=[0, ])

        variance_df.to_csv(path.join(output_dir, 'variance.tsv'),
                           index=False, sep='\t', encoding='utf-8')

    def _compute_variance(self, test_error_split):

        # compute average test error
        num_split = len(self._split_results)  # J in the paper

        # compute mu_{n_1}^{n_2}
        average_test_error = np.mean(test_error_split)

        approx_variance = np.sum((test_error_split - average_test_error)**2)/(num_split - 1)

        # compute variance (point 2 and 6 of Nadeau's paper)
        resampled_t = approx_variance / num_split
        corrected_resampled_t = (1/num_split + self._test_size/(1 - self._test_size)) * approx_variance

        return resampled_t, corrected_resampled_t

    def compute_error_variance(self):
        num_split = len(self._split_results)
        test_error_split = np.zeros((num_split, 1))  # this list will contain the list of mu_j hat for j = 1 to J
        for i in range(num_split):
            test_error_split[i] = self._compute_average_test_error(self._split_results[i]['y'],
                                                                   self._split_results[i]['y_hat'])

        self._error_resampled_t, self._error_corrected_resampled_t = self._compute_variance(test_error_split)

        return self._error_resampled_t, self._error_corrected_resampled_t

    def _compute_average_test_error(self, y_list, yhat_list):
        # return the average test error (denoted mu_j hat)
        return float(len(np.where(y_list != yhat_list)[0]))/float(len(y_list))

    def compute_accuracy_variance(self):
        num_split = len(self._split_results)
        test_accuracy_split = np.zeros((num_split, 1))  # this list will contain the list of mu_j hat for j = 1 to J
        for i in range(num_split):
            test_accuracy_split[i] = self._compute_average_test_accuracy(self._split_results[i]['y'],
                                                                         self._split_results[i]['y_hat'])

        self._bal_accuracy_resampled_t, self._bal_accuracy_corrected_resampled_t = self._compute_variance(test_accuracy_split)

        return self._bal_accuracy_resampled_t, self._bal_accuracy_corrected_resampled_t

    def _compute_average_test_accuracy(self, y_list, yhat_list):

        from clinica.pipelines.machine_learning.ml_utils import evaluate_prediction

        return evaluate_prediction(y_list, yhat_list)['balanced_accuracy']


class LearningCurveRepeatedHoldOut(base.MLValidation):

    def __init__(self, ml_algorithm, n_iterations=100, test_size=0.3, n_learning_points=10):
        self._ml_algorithm = ml_algorithm
        self._split_results = []
        self._classifier = None
        self._best_params = None
        self._cv = None
        self._n_iterations = n_iterations
        self._test_size = test_size
        self._n_learning_points = n_learning_points
        self._error_resampled_t = None
        self._error_corrected_resampled_t = None
        self._bal_accuracy_resampled_t = None
        self._bal_accuracy_corrected_resampled_t = None

    def validate(self, y, n_threads=15):

        splits = StratifiedShuffleSplit(n_splits=self._n_iterations, test_size=self._test_size)
        self._cv = list(splits.split(np.zeros(len(y)), y))
        async_pool = ThreadPool(n_threads)
        async_result = {}

        for i in range(self._n_iterations):
            train_index, test_index = self._cv[i]
            async_result[i] = {}

            skf = StratifiedKFold(n_splits=self._n_learning_points, shuffle=False)
            inner_cv = list(skf.split(np.zeros(len(y[train_index])), y[train_index]))

            for j in range(self._n_learning_points):
                inner_train_index = np.concatenate([indexes[1] for indexes in inner_cv[:j + 1]]).ravel()
                async_result[i][j] = async_pool.apply_async(self._ml_algorithm.evaluate, (train_index[inner_train_index], test_index))

        async_pool.close()
        async_pool.join()

        for j in range(self._n_learning_points):
            learning_point_results = []
            for i in range(self._n_iterations):
                learning_point_results.append(async_result[i][j].get())

            self._split_results.append(learning_point_results)

        self._classifier = []
        self._best_params = []
        for j in range(self._n_learning_points):
            classifier, best_params = self._ml_algorithm.apply_best_parameters(self._split_results[j])
            self._classifier.append(classifier)
            self._best_params.append(best_params)

        return self._classifier, self._best_params, self._split_results

    def save_results(self, output_dir):
        if self._split_results is None:
            raise Exception("No results to save. Method validate() must be run before save_results().")

        for learning_point in range(self._n_learning_points):

            all_results_list = []
            all_subjects_list = []

            learning_point_dir = path.join(output_dir, 'learning_split-' + str(learning_point))

            for iteration in range(self._n_iterations):

                iteration_dir = path.join(learning_point_dir, 'iteration-' + str(iteration))
                if not path.exists(iteration_dir):
                    os.makedirs(iteration_dir)
                iteration_subjects_df = pd.DataFrame({'y': self._split_results[learning_point][iteration]['y'],
                                                      'y_hat': self._split_results[learning_point][iteration]['y_hat'],
                                                      'y_index': self._split_results[learning_point][iteration]['y_index']})
                iteration_subjects_df.to_csv(path.join(iteration_dir, 'subjects.tsv'),
                                             index=False, sep='\t', encoding='utf-8')
                all_subjects_list.append(iteration_subjects_df)

                iteration_results_df = pd.DataFrame(
                        {'balanced_accuracy': self._split_results[learning_point][iteration]['evaluation']['balanced_accuracy'],
                         'auc': self._split_results[learning_point][iteration]['auc'],
                         'accuracy': self._split_results[learning_point][iteration]['evaluation']['accuracy'],
                         'sensitivity': self._split_results[learning_point][iteration]['evaluation']['sensitivity'],
                         'specificity': self._split_results[learning_point][iteration]['evaluation']['specificity'],
                         'ppv': self._split_results[learning_point][iteration]['evaluation']['ppv'],
                         'npv': self._split_results[learning_point][iteration]['evaluation']['npv'],
                         'train_balanced_accuracy': self._split_results[learning_point][iteration]['evaluation_train']['balanced_accuracy'],
                         'train_accuracy': self._split_results[learning_point][iteration]['evaluation_train']['accuracy'],
                         'train_sensitivity': self._split_results[learning_point][iteration]['evaluation_train']['sensitivity'],
                         'train_specificity': self._split_results[learning_point][iteration]['evaluation_train']['specificity'],
                         'train_ppv': self._split_results[learning_point][iteration]['evaluation_train']['ppv'],
                         'train_npv': self._split_results[learning_point][iteration]['evaluation_train']['npv']}, index=['i', ])

                iteration_results_df.to_csv(path.join(iteration_dir, 'results.tsv'),
                                            index=False, sep='\t', encoding='utf-8')

                mean_results_df = pd.DataFrame(iteration_results_df.apply(np.nanmean).to_dict(),
                                               columns=iteration_results_df.columns, index=[0, ])
                mean_results_df.to_csv(path.join(iteration_dir, 'mean_results.tsv'),
                                       index=False, sep='\t', encoding='utf-8')
                all_results_list.append(mean_results_df)

            all_subjects_df = pd.concat(all_subjects_list)
            all_subjects_df.to_csv(path.join(learning_point_dir, 'subjects.tsv'),
                                   index=False, sep='\t', encoding='utf-8')

            all_results_df = pd.concat(all_results_list)
            all_results_df.to_csv(path.join(learning_point_dir, 'results.tsv'),
                                  index=False, sep='\t', encoding='utf-8')

            mean_results_df = pd.DataFrame(all_results_df.apply(np.nanmean).to_dict(),
                                           columns=all_results_df.columns, index=[0, ])
            mean_results_df.to_csv(path.join(learning_point_dir, 'mean_results.tsv'),
                                   index=False, sep='\t', encoding='utf-8')

            self.compute_error_variance(learning_point)
            self.compute_accuracy_variance(learning_point)

            variance_df = pd.DataFrame({'bal_accuracy_resampled_t': self._bal_accuracy_resampled_t,
                                        'bal_accuracy_corrected_resampled_t': self._bal_accuracy_corrected_resampled_t,
                                        'error_resampled_t': self._error_resampled_t,
                                        'error_corrected_resampled_t': self._error_corrected_resampled_t}, index=[0, ])

            variance_df.to_csv(path.join(learning_point_dir, 'variance.tsv'),
                               index=False, sep='\t', encoding='utf-8')

    def _compute_variance(self, test_error_split):

        # compute average test error
        num_split = self._n_iterations  # J in the paper

        # compute mu_{n_1}^{n_2}
        average_test_error = np.mean(test_error_split)

        approx_variance = np.sum((test_error_split - average_test_error)**2)/(num_split - 1)

        # compute variance (point 2 and 6 of Nadeau's paper)
        resampled_t = approx_variance / num_split
        corrected_resampled_t = (1/num_split + self._test_size/(1 - self._test_size)) * approx_variance

        return resampled_t, corrected_resampled_t

    def compute_error_variance(self, learning_point):
        num_split = self._n_iterations
        test_error_split = np.zeros((num_split, 1))  # this list will contain the list of mu_j hat for j = 1 to J
        for i in range(num_split):
            test_error_split[i] = self._compute_average_test_error(self._split_results[learning_point][i]['y'],
                                                                   self._split_results[learning_point][i]['y_hat'])

            self._error_resampled_t, self._error_corrected_resampled_t = self._compute_variance(test_error_split)

        return self._error_resampled_t, self._error_corrected_resampled_t

    def _compute_average_test_error(self, y_list, yhat_list):
        # return the average test error (denoted mu_j hat)
        return float(len(np.where(y_list != yhat_list)[0]))/float(len(y_list))

    def compute_accuracy_variance(self, learning_point):
        num_split = self._n_iterations
        test_accuracy_split = np.zeros((num_split, 1))  # this list will contain the list of mu_j hat for j = 1 to J
        for i in range(num_split):
            test_accuracy_split[i] = self._compute_average_test_accuracy(self._split_results[learning_point][i]['y'],
                                                                         self._split_results[learning_point][i]['y_hat'])

        self._bal_accuracy_resampled_t, self._bal_accuracy_corrected_resampled_t = self._compute_variance(test_accuracy_split)

        return self._bal_accuracy_resampled_t, self._bal_accuracy_corrected_resampled_t

    def _compute_average_test_accuracy(self, y_list, yhat_list):

        from clinica.pipelines.machine_learning.ml_utils import evaluate_prediction

        return evaluate_prediction(y_list, yhat_list)['balanced_accuracy']


class RepeatedKFoldCV_Multiclass(base.MLValidation):

    def __init__(self, ml_algorithm):
        self._ml_algorithm = ml_algorithm
        self._repeated_fold_results = []
        self._classifier = None
        self._best_params = None
        self._cv = None

    def validate(self, y, n_iterations=100, n_folds=10, n_threads=15):

        async_pool = ThreadPool(n_threads)
        async_result = {}
        self._cv = []

        for r in range(n_iterations):
            skf = StratifiedKFold(n_splits=n_folds, shuffle=True)
            self._cv.append(list(skf.split(np.zeros(len(y)), y)))
            async_result[r] = {}
            self._repeated_fold_results.append([])

            for i in range(n_folds):

                train_index, test_index = self._cv[r][i]
                async_result[r][i] = async_pool.apply_async(self._ml_algorithm.evaluate, (train_index, test_index))

        async_pool.close()
        async_pool.join()
        for r in range(n_iterations):
            for i in range(n_folds):
                self._repeated_fold_results[r].append(async_result[r][i].get())

        # TODO Find a better way to estimate best parameter
        flat_results = [result for fold in self._repeated_fold_results for result in fold]
        self._classifier, self._best_params = self._ml_algorithm.apply_best_parameters(flat_results)

        return self._classifier, self._best_params, self._repeated_fold_results

    def save_results(self, output_dir):
        if self._repeated_fold_results is None:
            raise Exception("No results to save. Method validate() must be run before save_results().")

        all_results_list = []
        all_subjects_list = []

        for iteration in range(len(self._repeated_fold_results)):

            iteration_dir = path.join(output_dir, 'iteration-' + str(iteration))
            if not path.exists(iteration_dir):
                os.makedirs(iteration_dir)

            iteration_subjects_list = []
            iteration_results_list = []
            folds_dir = path.join(iteration_dir, 'folds')

            if not path.exists(folds_dir):
                os.makedirs(folds_dir)

            for i in range(len(self._repeated_fold_results[iteration])):
                subjects_df = pd.DataFrame({'y': self._repeated_fold_results[iteration][i]['y'],
                                            'y_hat': self._repeated_fold_results[iteration][i]['y_hat'],
                                            'y_index': self._repeated_fold_results[iteration][i]['y_index']})
                subjects_df.to_csv(path.join(folds_dir, 'subjects_fold-' + str(i) + '.tsv'),
                                   index=False, sep='\t', encoding='utf-8')
                iteration_subjects_list.append(subjects_df)

                results_df = pd.DataFrame(
                    {'balanced_accuracy': self._repeated_fold_results[iteration][i]['evaluation']['balanced_accuracy'],
                     'accuracy': self._repeated_fold_results[iteration][i]['evaluation']['accuracy'],
                     'train_balanced_accuracy': self._repeated_fold_results[iteration][i]['evaluation_train']['balanced_accuracy'],
                     'train_accuracy': self._repeated_fold_results[iteration][i]['evaluation_train']['accuracy']
                     }, index=['i', ])
                results_df.to_csv(path.join(folds_dir, 'results_fold-' + str(i) + '.tsv'),
                                  index=False, sep='\t', encoding='utf-8')
                iteration_results_list.append(results_df)

            iteration_subjects_df = pd.concat(iteration_subjects_list)
            iteration_subjects_df.to_csv(path.join(iteration_dir, 'subjects.tsv'),
                                         index=False, sep='\t', encoding='utf-8')
            all_subjects_list.append(iteration_subjects_df)

            iteration_results_df = pd.concat(iteration_results_list)
            iteration_results_df.to_csv(path.join(iteration_dir, 'results.tsv'),
                                        index=False, sep='\t', encoding='utf-8')

            mean_results_df = pd.DataFrame(iteration_results_df.apply(np.nanmean).to_dict(),
                                           columns=iteration_results_df.columns, index=[0, ])
            mean_results_df.to_csv(path.join(iteration_dir, 'mean_results.tsv'),
                                   index=False, sep='\t', encoding='utf-8')
            all_results_list.append(mean_results_df)

        all_subjects_df = pd.concat(all_subjects_list)
        all_subjects_df.to_csv(path.join(output_dir, 'subjects.tsv'),
                               index=False, sep='\t', encoding='utf-8')

        all_results_df = pd.concat(all_results_list)
        all_results_df.to_csv(path.join(output_dir, 'results.tsv'),
                              index=False, sep='\t', encoding='utf-8')

        mean_results_df = pd.DataFrame(all_results_df.apply(np.nanmean).to_dict(),
                                       columns=all_results_df.columns, index=[0, ])
        mean_results_df.to_csv(path.join(output_dir, 'mean_results.tsv'),
                               index=False, sep='\t', encoding='utf-8')

        print("Mean results of the classification:")
        print("Balanced accuracy: %s" % (mean_results_df['balanced_accuracy'].to_string(index=False)))
