
from os import path
import json
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
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
        for i in range(len(self._fold_results)):
            subjects_df = pd.DataFrame({'y': self._fold_results[i]['y'],
                                        'y_hat': self._fold_results[i]['y_hat'],
                                        'y_index': self._fold_results[i]['y_index']})

            subjects_df.to_csv(path.join(output_dir, 'subjects_fold-' + str(i) + '.tsv'),
                               index=False, sep='\t', encoding='utf-8')

            subjects_folds.append(subjects_df)

        all_subjects = pd.concat(subjects_folds)
        all_subjects.to_csv(path.join(output_dir, 'subjects.tsv'),
                            index=False, sep='\t', encoding='utf-8')
