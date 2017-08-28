
import os
from os import path
import json
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold, StratifiedShuffleSplit
from multiprocessing.pool import ThreadPool

from clinica.pipeline.machine_learning import base


class KFoldCV(base.MLValidation):

    def __init__(self, ml_algorithm):
        self._ml_algorithm = ml_algorithm
        self._fold_results = []
        self._classifier = None
        self._best_params = None
        self._cv = None

    def validate(self, y, n_folds=10, n_threads=15):

        skf = StratifiedKFold(n_splits=n_folds, shuffle=True)
        self._cv = list(skf.split(np.zeros(len(y)), y))
        async_pool = ThreadPool(n_threads)
        async_result = {}

        for i in range(n_folds):

            train_index, test_index = self._cv[i]

            # self._fold_results.append(self._ml_algorithm.parameter_estimation(train_index))

            async_result[i] = async_pool.apply_async(self._ml_algorithm.evaluate, (train_index, test_index))

        # print "lanzando threads"
        async_pool.close()
        async_pool.join()
        # print "terminaron threads"

        for i in range(n_folds):
            self._fold_results.append(async_result[i].get())

        self._classifier, self._best_params = self._ml_algorithm.apply_best_parameters(self._fold_results)

        return self._classifier, self._best_params, self._fold_results

    def save_results(self, output_dir):
        if self._fold_results is None:
            raise Exception("No results to save. Method validate() must be run before save_results().")

        results_dict = {'balanced_accuracy': np.nanmean([r['evaluation']['balanced_accuracy'] for r in self._fold_results]),
                        'auc': np.nanmean([r['auc'] for r in self._fold_results]),
                        'accuracy': np.nanmean([r['evaluation']['accuracy'] for r in self._fold_results]),
                        'sensitivity': np.nanmean([r['evaluation']['sensitivity'] for r in self._fold_results]),
                        'specificity': np.nanmean([r['evaluation']['specificity'] for r in self._fold_results]),
                        'ppv': np.nanmean([r['evaluation']['ppv'] for r in self._fold_results]),
                        'npv': np.nanmean([r['evaluation']['npv'] for r in self._fold_results])}

        # t1_df = pd.DataFrame(columns=t1_col_df)
        # t1_df = t1_df.append(row_to_append, ignore_index=True)

        results_df = pd.DataFrame(results_dict, index=['i', ])
        results_df.to_csv(path.join(output_dir, 'results.tsv'),
                          index=False, sep='\t', encoding='utf-8')

        subjects_folds = []
        container_dir = path.join(output_dir, 'folds')
        os.mkdir(container_dir)
        for i in range(len(self._fold_results)):
            subjects_df = pd.DataFrame({'y': self._fold_results[i]['y'],
                                        'y_hat': self._fold_results[i]['y_hat'],
                                        'y_index': self._fold_results[i]['y_index']})

            subjects_df.to_csv(path.join(container_dir, 'subjects_fold-' + str(i) + '.tsv'),
                               index=False, sep='\t', encoding='utf-8')

            subjects_folds.append(subjects_df)

        all_subjects = pd.concat(subjects_folds)
        all_subjects.to_csv(path.join(output_dir, 'subjects.tsv'),
                            index=False, sep='\t', encoding='utf-8')


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

        all_iterations_results = []
        all_fold_subjects = []
        for iteration in range(len(self._repeated_fold_results)):

            iteration_dir = path.join(output_dir, 'iteration-' + str(iteration))
            os.mkdir(iteration_dir)

            iteration_results_dict = {'balanced_accuracy': np.nanmean([r['evaluation']['balanced_accuracy'] for r in self._repeated_fold_results[iteration]]),
                                      'auc': np.nanmean([r['auc'] for r in self._repeated_fold_results[iteration]]),
                                      'accuracy': np.nanmean([r['evaluation']['accuracy'] for r in self._repeated_fold_results[iteration]]),
                                      'sensitivity': np.nanmean([r['evaluation']['sensitivity'] for r in self._repeated_fold_results[iteration]]),
                                      'specificity': np.nanmean([r['evaluation']['specificity'] for r in self._repeated_fold_results[iteration]]),
                                      'ppv': np.nanmean([r['evaluation']['ppv'] for r in self._repeated_fold_results[iteration]]),
                                      'npv': np.nanmean([r['evaluation']['npv'] for r in self._repeated_fold_results[iteration]])}

            iteration_results_df = pd.DataFrame(iteration_results_dict, index=['i', ])
            iteration_results_df.to_csv(path.join(iteration_dir, 'results_iteration-' + str(iteration) + '.tsv'),
                                        index=False, sep='\t', encoding='utf-8')

            all_iterations_results.append(iteration_results_df)

            fold_subjects = []
            container_dir = path.join(iteration_dir, 'folds')
            os.mkdir(container_dir)
            for i in range(len(self._repeated_fold_results[iteration])):
                subjects_df = pd.DataFrame({'y': self._repeated_fold_results[iteration][i]['y'],
                                            'y_hat': self._repeated_fold_results[iteration][i]['y_hat'],
                                            'y_index': self._repeated_fold_results[iteration][i]['y_index']})

                subjects_df.to_csv(path.join(container_dir, 'subjects_fold-' + str(i) + '.tsv'),
                                   index=False, sep='\t', encoding='utf-8')

                fold_subjects.append(subjects_df)
                all_fold_subjects.append(subjects_df)

            iteration_subjects = pd.concat(fold_subjects)
            iteration_subjects.to_csv(path.join(iteration_dir, 'subjects_iteration-' + str(iteration) + '.tsv'),
                                      index=False, sep='\t', encoding='utf-8')

        all_results = pd.concat(all_iterations_results)
        all_results.to_csv(path.join(output_dir, 'results.tsv'),
                           index=False, sep='\t', encoding='utf-8')

        all_subjects = pd.concat(all_fold_subjects)
        all_subjects.to_csv(path.join(output_dir, 'subjects.tsv'),
                            index=False, sep='\t', encoding='utf-8')


class RepeatedSplit(base.MLValidation):

    def __init__(self, ml_algorithm):
        self._ml_algorithm = ml_algorithm
        self._split_results = []
        self._classifier = None
        self._best_params = None
        self._cv = None

    def validate(self, y, n_iterations=100, test_size=0.3, n_threads=15):

        splits = StratifiedShuffleSplit(n_splits=n_iterations, test_size=test_size)
        self._cv = list(splits.split(np.zeros(len(y)), y))
        async_pool = ThreadPool(n_threads)
        async_result = {}

        for i in range(n_iterations):

            train_index, test_index = self._cv[i]
            async_result[i] = async_pool.apply_async(self._ml_algorithm.evaluate, (train_index, test_index))

        async_pool.close()
        async_pool.join()

        for i in range(n_iterations):
            self._split_results.append(async_result[i].get())

        self._classifier, self._best_params = self._ml_algorithm.apply_best_parameters(self._split_results)

        return self._classifier, self._best_params, self._split_results

    def save_results(self, output_dir):
        if self._split_results is None:
            raise Exception("No results to save. Method validate() must be run before save_results().")

        results_dict = {'balanced_accuracy': np.nanmean([r['evaluation']['balanced_accuracy'] for r in self._split_results]),
                        'auc': np.nanmean([r['auc'] for r in self._split_results]),
                        'accuracy': np.nanmean([r['evaluation']['accuracy'] for r in self._split_results]),
                        'sensitivity': np.nanmean([r['evaluation']['sensitivity'] for r in self._split_results]),
                        'specificity': np.nanmean([r['evaluation']['specificity'] for r in self._split_results]),
                        'ppv': np.nanmean([r['evaluation']['ppv'] for r in self._split_results]),
                        'npv': np.nanmean([r['evaluation']['npv'] for r in self._split_results])}

        # t1_df = pd.DataFrame(columns=t1_col_df)
        # t1_df = t1_df.append(row_to_append, ignore_index=True)

        results_df = pd.DataFrame(results_dict, index=['i', ])
        results_df.to_csv(path.join(output_dir, 'results.tsv'),
                          index=False, sep='\t', encoding='utf-8')

        subjects_folds = []
        container_dir = path.join(output_dir, 'iterations')
        os.mkdir(container_dir)
        for i in range(len(self._split_results)):
            subjects_df = pd.DataFrame({'y': self._split_results[i]['y'],
                                        'y_hat': self._split_results[i]['y_hat'],
                                        'y_index': self._split_results[i]['y_index']})

            subjects_df.to_csv(path.join(container_dir, 'subjects_iteration-' + str(i) + '.tsv'),
                               index=False, sep='\t', encoding='utf-8')

            subjects_folds.append(subjects_df)

        all_subjects = pd.concat(subjects_folds)
        all_subjects.to_csv(path.join(output_dir, 'subjects.tsv'),
                            index=False, sep='\t', encoding='utf-8')

