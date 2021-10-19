import os
from multiprocessing.pool import ThreadPool
from os import path

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold, StratifiedShuffleSplit

from clinica.pipelines.machine_learning import base


class KFoldCV(base.MLValidation):
    def validate(self, y):

        if not self._validation_params["splits_indices"]:
            skf = StratifiedKFold(
                n_splits=self._validation_params["n_folds"], shuffle=True
            )
            self._validation_params["splits_indices"] = list(
                skf.split(np.zeros(len(y)), y)
            )

        async_pool = ThreadPool(self._validation_params["n_threads"])
        async_result = {}

        for i in range(self._validation_params["n_folds"]):

            train_index, test_index = self._validation_params["splits_indices"][i]
            async_result[i] = async_pool.apply_async(
                self._ml_algorithm.evaluate, (train_index, test_index)
            )

        async_pool.close()
        async_pool.join()

        for i in range(self._validation_params["n_folds"]):
            self._validation_results.append(async_result[i].get())

        self._classifier, self._best_params = self._ml_algorithm.apply_best_parameters(
            self._validation_results
        )

        return self._classifier, self._best_params, self._validation_results

    def save_results(self, output_dir: str):
        if not self._validation_results:
            raise Exception(
                "No results to save. Method validate() must be run before save_results()."
            )

        subjects_folds = []
        results_folds = []
        container_dir = path.join(output_dir, "folds")

        os.makedirs(container_dir, exist_ok=True)

        for i in range(len(self._validation_results)):
            subjects_df = pd.DataFrame(
                {
                    "y": self._validation_results[i]["y"],
                    "y_hat": self._validation_results[i]["y_hat"],
                    "y_index": self._validation_results[i]["y_index"],
                }
            )
            subjects_df.to_csv(
                path.join(container_dir, "subjects_fold-" + str(i) + ".tsv"),
                index=False,
                sep="\t",
                encoding="utf-8",
            )
            subjects_folds.append(subjects_df)

            # fmt: off
            results_df = pd.DataFrame(
                {
                    "balanced_accuracy": self._validation_results[i]["evaluation"]["balanced_accuracy"],
                    "auc": self._validation_results[i]["auc"],
                    "accuracy": self._validation_results[i]["evaluation"]["accuracy"],
                    "sensitivity": self._validation_results[i]["evaluation"]["sensitivity"],
                    "specificity": self._validation_results[i]["evaluation"]["specificity"],
                    "ppv": self._validation_results[i]["evaluation"]["ppv"],
                    "npv": self._validation_results[i]["evaluation"]["npv"],
                    "train_balanced_accuracy": self._validation_results[i]["evaluation_train"]["balanced_accuracy"],
                    "train_accuracy": self._validation_results[i]["evaluation_train"]["accuracy"],
                    "train_sensitivity": self._validation_results[i]["evaluation_train"]["sensitivity"],
                    "train_specificity": self._validation_results[i]["evaluation_train"]["specificity"],
                    "train_ppv": self._validation_results[i]["evaluation_train"]["ppv"],
                    "train_npv": self._validation_results[i]["evaluation_train"]["npv"],
                },
                index=["i"],
            )
            # fmt: on

            results_df.to_csv(
                path.join(container_dir, "results_fold-" + str(i) + ".tsv"),
                index=False,
                sep="\t",
                encoding="utf-8",
            )
            results_folds.append(results_df)

        all_subjects = pd.concat(subjects_folds)
        all_subjects.to_csv(
            path.join(output_dir, "subjects.tsv"),
            index=False,
            sep="\t",
            encoding="utf-8",
        )

        all_results = pd.concat(results_folds)
        all_results.to_csv(
            path.join(output_dir, "results.tsv"),
            index=False,
            sep="\t",
            encoding="utf-8",
        )

        mean_results = pd.DataFrame(
            all_results.apply(np.nanmean).to_dict(),
            columns=all_results.columns,
            index=[
                0,
            ],
        )
        mean_results.to_csv(
            path.join(output_dir, "mean_results.tsv"),
            index=False,
            sep="\t",
            encoding="utf-8",
        )

        print("Mean results of the classification:")
        print(
            "Balanced accuracy: %s"
            % (mean_results["balanced_accuracy"].to_string(index=False))
        )
        print("specificity: %s" % (mean_results["specificity"].to_string(index=False)))
        print("sensitivity: %s" % (mean_results["sensitivity"].to_string(index=False)))
        print("auc: %s" % (mean_results["auc"].to_string(index=False)))

    @staticmethod
    def get_default_parameters():

        parameters_dict = {
            "n_folds": 10,
            "n_threads": 15,
            "splits_indices": None,
            "inner_cv": True,
        }

        return parameters_dict


class RepeatedKFoldCV(base.MLValidation):
    def validate(self, y):

        if not self._validation_params["splits_indices"]:
            self._validation_params["splits_indices"] = []

            for i in range(self._validation_params["n_iterations"]):
                skf = StratifiedKFold(
                    n_splits=self._validation_params["n_folds"], shuffle=True
                )
                self._validation_params["splits_indices"].append(
                    list(skf.split(np.zeros(len(y)), y))
                )

        async_pool = ThreadPool(self._validation_params["n_threads"])
        async_result = {}

        for i in range(self._validation_params["n_iterations"]):

            train_index, test_index = self._validation_params["splits_indices"][i]
            async_result[i] = async_pool.apply_async(
                self._ml_algorithm.evaluate, (train_index, test_index)
            )

        for r in range(self._validation_params["n_iterations"]):

            async_result[r] = {}
            self._validation_results.append([])

            for i in range(self._validation_params["n_folds"]):

                train_index, test_index = self._validation_params["splits_indices"][r][
                    i
                ]
                async_result[r][i] = async_pool.apply_async(
                    self._ml_algorithm.evaluate, (train_index, test_index)
                )

        async_pool.close()
        async_pool.join()
        for r in range(self._validation_params["n_iterations"]):
            for i in range(self._validation_params["n_folds"]):
                self._validation_results[r].append(async_result[r][i].get())

        # TODO Find a better way to estimate best parameter
        flat_results = [result for fold in self._validation_results for result in fold]
        self._classifier, self._best_params = self._ml_algorithm.apply_best_parameters(
            flat_results
        )

        return self._classifier, self._best_params, self._validation_results

    def save_results(self, output_dir: str):
        if not self._validation_results:
            raise Exception(
                "No results to save. Method validate() must be run before save_results()."
            )

        all_results_list = []
        all_subjects_list = []

        for iteration in range(len(self._validation_results)):

            iteration_dir = path.join(output_dir, "iteration-" + str(iteration))
            os.makedirs(iteration_dir, exist_ok=True)

            iteration_subjects_list = []
            iteration_results_list = []
            folds_dir = path.join(iteration_dir, "folds")

            os.makedirs(folds_dir, exist_ok=True)

            for i in range(len(self._validation_results[iteration])):
                subjects_df = pd.DataFrame(
                    {
                        "y": self._validation_results[iteration][i]["y"],
                        "y_hat": self._validation_results[iteration][i]["y_hat"],
                        "y_index": self._validation_results[iteration][i]["y_index"],
                    }
                )
                subjects_df.to_csv(
                    path.join(folds_dir, "subjects_fold-" + str(i) + ".tsv"),
                    index=False,
                    sep="\t",
                    encoding="utf-8",
                )
                iteration_subjects_list.append(subjects_df)

                # fmt: off
                results_df = pd.DataFrame(
                    {
                        "balanced_accuracy": self._validation_results[iteration][i]["evaluation"]["balanced_accuracy"],
                        "auc": self._validation_results[iteration][i]["auc"],
                        "accuracy": self._validation_results[iteration][i]["evaluation"]["accuracy"],
                        "sensitivity": self._validation_results[iteration][i]["evaluation"]["sensitivity"],
                        "specificity": self._validation_results[iteration][i]["evaluation"]["specificity"],
                        "ppv": self._validation_results[iteration][i]["evaluation"]["ppv"],
                        "npv": self._validation_results[iteration][i]["evaluation"]["npv"],
                        "train_balanced_accuracy": self._validation_results[iteration][i]["evaluation_train"]["balanced_accuracy"],
                        "train_accuracy": self._validation_results[iteration][i]["evaluation_train"]["accuracy"],
                        "train_sensitivity": self._validation_results[iteration][i]["evaluation_train"]["sensitivity"],
                        "train_specificity": self._validation_results[iteration][i]["evaluation_train"]["specificity"],
                        "train_ppv": self._validation_results[iteration][i]["evaluation_train"]["ppv"],
                        "train_npv": self._validation_results[iteration][i]["evaluation_train"]["npv"],
                    },
                    index=["i"],
                )
                # fmt: on

                results_df.to_csv(
                    path.join(folds_dir, "results_fold-" + str(i) + ".tsv"),
                    index=False,
                    sep="\t",
                    encoding="utf-8",
                )
                iteration_results_list.append(results_df)

            iteration_subjects_df = pd.concat(iteration_subjects_list)
            iteration_subjects_df.to_csv(
                path.join(iteration_dir, "subjects.tsv"),
                index=False,
                sep="\t",
                encoding="utf-8",
            )
            all_subjects_list.append(iteration_subjects_df)

            iteration_results_df = pd.concat(iteration_results_list)
            iteration_results_df.to_csv(
                path.join(iteration_dir, "results.tsv"),
                index=False,
                sep="\t",
                encoding="utf-8",
            )

            mean_results_df = pd.DataFrame(
                iteration_results_df.apply(np.nanmean).to_dict(),
                columns=iteration_results_df.columns,
                index=[
                    0,
                ],
            )
            mean_results_df.to_csv(
                path.join(iteration_dir, "mean_results.tsv"),
                index=False,
                sep="\t",
                encoding="utf-8",
            )
            all_results_list.append(mean_results_df)

        all_subjects_df = pd.concat(all_subjects_list)
        all_subjects_df.to_csv(
            path.join(output_dir, "subjects.tsv"),
            index=False,
            sep="\t",
            encoding="utf-8",
        )

        all_results_df = pd.concat(all_results_list)
        all_results_df.to_csv(
            path.join(output_dir, "results.tsv"),
            index=False,
            sep="\t",
            encoding="utf-8",
        )

        mean_results_df = pd.DataFrame(
            all_results_df.apply(np.nanmean).to_dict(),
            columns=all_results_df.columns,
            index=[
                0,
            ],
        )
        mean_results_df.to_csv(
            path.join(output_dir, "mean_results.tsv"),
            index=False,
            sep="\t",
            encoding="utf-8",
        )

        print("Mean results of the classification:")
        print(
            "Balanced accuracy: %s"
            % (mean_results_df["balanced_accuracy"].to_string(index=False))
        )
        print(
            "specificity: %s" % (mean_results_df["specificity"].to_string(index=False))
        )
        print(
            "sensitivity: %s" % (mean_results_df["sensitivity"].to_string(index=False))
        )
        print("auc: %s" % (mean_results_df["auc"].to_string(index=False)))

    @staticmethod
    def get_default_parameters():

        parameters_dict = {
            "n_iterations": 100,
            "n_folds": 10,
            "n_threads": 15,
            "splits_indices": None,
            "inner_cv": True,
        }

        return parameters_dict


class RepeatedHoldOut(base.MLValidation):
    def validate(self, y):

        if not self._validation_params["splits_indices"]:
            splits = StratifiedShuffleSplit(
                n_splits=self._validation_params["n_iterations"],
                test_size=self._validation_params["test_size"],
            )
            self._validation_params["splits_indices"] = list(
                splits.split(np.zeros(len(y)), y)
            )

        async_pool = ThreadPool(self._validation_params["n_threads"])
        async_result = {}

        for i in range(self._validation_params["n_iterations"]):

            train_index, test_index = self._validation_params["splits_indices"][i]
            if self._validation_params["inner_cv"]:
                async_result[i] = async_pool.apply_async(
                    self._ml_algorithm.evaluate, (train_index, test_index)
                )
            else:
                async_result[i] = async_pool.apply_async(
                    self._ml_algorithm.evaluate_no_cv, (train_index, test_index)
                )

        async_pool.close()
        async_pool.join()

        for i in range(self._validation_params["n_iterations"]):
            self._validation_results.append(async_result[i].get())

        self._classifier, self._best_params = self._ml_algorithm.apply_best_parameters(
            self._validation_results
        )
        return self._classifier, self._best_params, self._validation_results

    def save_results(self, output_dir: str):
        if not self._validation_results:
            raise Exception(
                "No results to save. Method validate() must be run before save_results()."
            )

        all_results_list = []
        all_train_subjects_list = []
        all_test_subjects_list = []

        for iteration in range(len(self._validation_results)):

            iteration_dir = path.join(output_dir, "iteration-" + str(iteration))
            os.makedirs(iteration_dir, exist_ok=True)
            iteration_train_subjects_df = pd.DataFrame(
                {
                    "iteration": iteration,
                    "y": self._validation_results[iteration]["y_train"],
                    "y_hat": self._validation_results[iteration]["y_hat_train"],
                    "subject_index": self._validation_results[iteration]["x_index"],
                }
            )
            iteration_train_subjects_df.to_csv(
                path.join(iteration_dir, "train_subjects.tsv"),
                index=False,
                sep="\t",
                encoding="utf-8",
            )
            all_train_subjects_list.append(iteration_train_subjects_df)

            iteration_test_subjects_df = pd.DataFrame(
                {
                    "iteration": iteration,
                    "y": self._validation_results[iteration]["y"],
                    "y_hat": self._validation_results[iteration]["y_hat"],
                    "subject_index": self._validation_results[iteration]["y_index"],
                }
            )
            iteration_test_subjects_df.to_csv(
                path.join(iteration_dir, "test_subjects.tsv"),
                index=False,
                sep="\t",
                encoding="utf-8",
            )
            all_test_subjects_list.append(iteration_test_subjects_df)

            # fmt: off
            iteration_results_df = pd.DataFrame(
                {
                    "balanced_accuracy": self._validation_results[iteration]["evaluation"]["balanced_accuracy"],
                    "auc": self._validation_results[iteration]["auc"],
                    "accuracy": self._validation_results[iteration]["evaluation"]["accuracy"],
                    "sensitivity": self._validation_results[iteration]["evaluation"]["sensitivity"],
                    "specificity": self._validation_results[iteration]["evaluation"]["specificity"],
                    "ppv": self._validation_results[iteration]["evaluation"]["ppv"],
                    "npv": self._validation_results[iteration]["evaluation"]["npv"],
                    "train_balanced_accuracy": self._validation_results[iteration]["evaluation_train"]["balanced_accuracy"],
                    "train_accuracy": self._validation_results[iteration]["evaluation_train"]["accuracy"],
                    "train_sensitivity": self._validation_results[iteration]["evaluation_train"]["sensitivity"],
                    "train_specificity": self._validation_results[iteration]["evaluation_train"]["specificity"],
                    "train_ppv": self._validation_results[iteration]["evaluation_train"]["ppv"],
                    "train_npv": self._validation_results[iteration]["evaluation_train"]["npv"],
                },
                index=["i"],
            )
            # fmt: on

            iteration_results_df.to_csv(
                path.join(iteration_dir, "results.tsv"),
                index=False,
                sep="\t",
                encoding="utf-8",
            )

            # mean_results_df = pd.DataFrame(iteration_results_df.apply(np.nanmean).to_dict(),
            #                                columns=iteration_results_df.columns, index=[0, ])
            # mean_results_df.to_csv(path.join(iteration_dir, 'mean_results.tsv'),
            #                        index=False, sep='\t', encoding='utf-8')
            all_results_list.append(iteration_results_df)

        all_train_subjects_df = pd.concat(all_train_subjects_list)
        all_train_subjects_df.to_csv(
            path.join(output_dir, "train_subjects.tsv"),
            index=False,
            sep="\t",
            encoding="utf-8",
        )

        all_test_subjects_df = pd.concat(all_test_subjects_list)
        all_test_subjects_df.to_csv(
            path.join(output_dir, "test_subjects.tsv"),
            index=False,
            sep="\t",
            encoding="utf-8",
        )

        all_results_df = pd.concat(all_results_list)
        all_results_df.to_csv(
            path.join(output_dir, "results.tsv"),
            index=False,
            sep="\t",
            encoding="utf-8",
        )

        mean_results_df = pd.DataFrame(
            all_results_df.apply(np.nanmean).to_dict(),
            columns=all_results_df.columns,
            index=[
                0,
            ],
        )
        mean_results_df.to_csv(
            path.join(output_dir, "mean_results.tsv"),
            index=False,
            sep="\t",
            encoding="utf-8",
        )

        print("Mean results of the classification:")
        print(
            "Balanced accuracy: %s"
            % (mean_results_df["balanced_accuracy"].to_string(index=False))
        )
        print(
            "specificity: %s" % (mean_results_df["specificity"].to_string(index=False))
        )
        print(
            "sensitivity: %s" % (mean_results_df["sensitivity"].to_string(index=False))
        )
        print("auc: %s" % (mean_results_df["auc"].to_string(index=False)))

    @staticmethod
    def get_default_parameters():

        parameters_dict = {
            "n_iterations": 100,
            "test_size": 0.2,
            "n_threads": 15,
            "splits_indices": None,
            "inner_cv": True,
        }

        return parameters_dict


class LearningCurveRepeatedHoldOut(base.MLValidation):
    def validate(self, y):

        if not self._validation_params["splits_indices"]:
            splits = StratifiedShuffleSplit(
                n_splits=self._validation_params["n_iterations"],
                test_size=self._validation_params["test_size"],
            )
            self._validation_params["splits_indices"] = list(
                splits.split(np.zeros(len(y)), y)
            )

        async_pool = ThreadPool(self._validation_params["n_threads"])
        async_result = {}

        for i in range(self._validation_params["n_iterations"]):
            train_index, test_index = self._validation_params["splits_indices"][i]
            async_result[i] = {}

            skf = StratifiedKFold(
                n_splits=self._validation_params["n_learning_points"], shuffle=False
            )
            inner_cv_splits = list(
                skf.split(np.zeros(len(y[train_index])), y[train_index])
            )

            for j in range(self._validation_params["n_learning_points"]):
                inner_train_index = np.concatenate(
                    [indexes[1] for indexes in inner_cv_splits[: j + 1]]
                ).ravel()
                async_result[i][j] = async_pool.apply_async(
                    self._ml_algorithm.evaluate,
                    (train_index[inner_train_index], test_index),
                )

        async_pool.close()
        async_pool.join()

        for j in range(self._validation_params["n_learning_points"]):
            learning_point_results = []
            for i in range(self._validation_params["n_iterations"]):
                learning_point_results.append(async_result[i][j].get())

            self._validation_results.append(learning_point_results)

        self._classifier = []
        self._best_params = []
        for j in range(self._validation_params["n_learning_points"]):
            classifier, best_params = self._ml_algorithm.apply_best_parameters(
                self._validation_results[j]
            )
            self._classifier.append(classifier)
            self._best_params.append(best_params)

        return self._classifier, self._best_params, self._validation_results

    def save_results(self, output_dir: str):
        if not self._validation_results:
            raise Exception(
                "No results to save. Method validate() must be run before save_results()."
            )

        for learning_point in range(self._validation_params["n_learning_points"]):

            all_results_list = []
            all_subjects_list = []

            learning_point_dir = path.join(
                output_dir, "learning_split-" + str(learning_point)
            )

            for iteration in range(self._validation_params["n_iterations"]):

                iteration_dir = path.join(
                    learning_point_dir, "iteration-" + str(iteration)
                )
                os.makedirs(iteration_dir, exist_ok=True)
                iteration_subjects_df = pd.DataFrame(
                    {
                        "y": self._validation_results[learning_point][iteration]["y"],
                        "y_hat": self._validation_results[learning_point][iteration][
                            "y_hat"
                        ],
                        "y_index": self._validation_results[learning_point][iteration][
                            "y_index"
                        ],
                    }
                )
                iteration_subjects_df.to_csv(
                    path.join(iteration_dir, "subjects.tsv"),
                    index=False,
                    sep="\t",
                    encoding="utf-8",
                )
                all_subjects_list.append(iteration_subjects_df)

                # fmt: off
                iteration_results_df = pd.DataFrame(
                    {
                        "balanced_accuracy": self._validation_results[learning_point][iteration]["evaluation"]["balanced_accuracy"],
                        "auc": self._validation_results[learning_point][iteration]["auc"],
                        "accuracy": self._validation_results[learning_point][iteration]["evaluation"]["accuracy"],
                        "sensitivity": self._validation_results[learning_point][iteration]["evaluation"]["sensitivity"],
                        "specificity": self._validation_results[learning_point][iteration]["evaluation"]["specificity"],
                        "ppv": self._validation_results[learning_point][iteration]["evaluation"]["ppv"],
                        "npv": self._validation_results[learning_point][iteration]["evaluation"]["npv"],
                        "train_balanced_accuracy": self._validation_results[learning_point][iteration]["evaluation_train"]["balanced_accuracy"],
                        "train_accuracy": self._validation_results[learning_point][iteration]["evaluation_train"]["accuracy"],
                        "train_sensitivity": self._validation_results[learning_point][iteration]["evaluation_train"]["sensitivity"],
                        "train_specificity": self._validation_results[learning_point][iteration]["evaluation_train"]["specificity"],
                        "train_ppv": self._validation_results[learning_point][iteration]["evaluation_train"]["ppv"],
                        "train_npv": self._validation_results[learning_point][iteration]["evaluation_train"]["npv"],
                    },
                    index=["i"],
                )
                # fmt: on

                iteration_results_df.to_csv(
                    path.join(iteration_dir, "results.tsv"),
                    index=False,
                    sep="\t",
                    encoding="utf-8",
                )

                mean_results_df = pd.DataFrame(
                    iteration_results_df.apply(np.nanmean).to_dict(),
                    columns=iteration_results_df.columns,
                    index=[
                        0,
                    ],
                )
                mean_results_df.to_csv(
                    path.join(iteration_dir, "mean_results.tsv"),
                    index=False,
                    sep="\t",
                    encoding="utf-8",
                )
                all_results_list.append(mean_results_df)

            all_subjects_df = pd.concat(all_subjects_list)
            all_subjects_df.to_csv(
                path.join(learning_point_dir, "subjects.tsv"),
                index=False,
                sep="\t",
                encoding="utf-8",
            )

            all_results_df = pd.concat(all_results_list)
            all_results_df.to_csv(
                path.join(learning_point_dir, "results.tsv"),
                index=False,
                sep="\t",
                encoding="utf-8",
            )

            mean_results_df = pd.DataFrame(
                all_results_df.apply(np.nanmean).to_dict(),
                columns=all_results_df.columns,
                index=[
                    0,
                ],
            )
            mean_results_df.to_csv(
                path.join(learning_point_dir, "mean_results.tsv"),
                index=False,
                sep="\t",
                encoding="utf-8",
            )

    @staticmethod
    def get_default_parameters():

        parameters_dict = {
            "n_iterations": 100,
            "test_size": 0.2,
            "n_learning_points": 10,
            "n_threads": 15,
            "splits_indices": None,
            "inner_cv": True,
        }

        return parameters_dict


class RepeatedKFoldCV_Multiclass(base.MLValidation):
    def __init__(self, ml_algorithm):
        self._ml_algorithm = ml_algorithm
        self._repeated_validation_results = []
        self._classifier = None
        self._best_params = None
        self._cv = None

    def validate(self, y, n_iterations=100, n_folds=10, n_threads=15):

        async_pool = ThreadPool(self._validation_params["n_threads"])
        async_result = {}
        self._cv = []

        for r in range(n_iterations):
            skf = StratifiedKFold(n_splits=n_folds, shuffle=True)
            self._cv.append(list(skf.split(np.zeros(len(y)), y)))
            async_result[r] = {}
            self._repeated_validation_results.append([])

            for i in range(n_folds):

                train_index, test_index = self._cv[r][i]
                async_result[r][i] = async_pool.apply_async(
                    self._ml_algorithm.evaluate, (train_index, test_index)
                )

        async_pool.close()
        async_pool.join()
        for r in range(n_iterations):
            for i in range(n_folds):
                self._repeated_validation_results[r].append(async_result[r][i].get())

        # TODO Find a better way to estimate best parameter
        flat_results = [
            result for fold in self._repeated_validation_results for result in fold
        ]
        self._classifier, self._best_params = self._ml_algorithm.apply_best_parameters(
            flat_results
        )

        return self._classifier, self._best_params, self._repeated_validation_results

    def save_results(self, output_dir: str):
        if not self._repeated_validation_results:
            raise Exception(
                "No results to save. Method validate() must be run before save_results()."
            )

        all_results_list = []
        all_subjects_list = []

        for iteration in range(len(self._repeated_validation_results)):

            iteration_dir = path.join(output_dir, "iteration-" + str(iteration))
            os.makedirs(iteration_dir, exist_ok=True)

            iteration_subjects_list = []
            iteration_results_list = []
            folds_dir = path.join(iteration_dir, "folds")

            os.makedirs(folds_dir, exist_ok=True)

            for i in range(len(self._repeated_validation_results[iteration])):
                subjects_df = pd.DataFrame(
                    {
                        "y": self._repeated_validation_results[iteration][i]["y"],
                        "y_hat": self._repeated_validation_results[iteration][i][
                            "y_hat"
                        ],
                        "y_index": self._repeated_validation_results[iteration][i][
                            "y_index"
                        ],
                    }
                )
                subjects_df.to_csv(
                    path.join(folds_dir, "subjects_fold-" + str(i) + ".tsv"),
                    index=False,
                    sep="\t",
                    encoding="utf-8",
                )
                iteration_subjects_list.append(subjects_df)

                # fmt: off
                results_df = pd.DataFrame(
                    {
                        "balanced_accuracy": self._repeated_validation_results[iteration][i]["evaluation"]["balanced_accuracy"],
                        "accuracy": self._repeated_validation_results[iteration][i]["evaluation"]["accuracy"],
                        "train_balanced_accuracy": self._repeated_validation_results[iteration][i]["evaluation_train"]["balanced_accuracy"],
                        "train_accuracy": self._repeated_validation_results[iteration][i]["evaluation_train"]["accuracy"],
                    },
                    index=[
                        "i",
                    ],
                )
                # fmt: on

                results_df.to_csv(
                    path.join(folds_dir, "results_fold-" + str(i) + ".tsv"),
                    index=False,
                    sep="\t",
                    encoding="utf-8",
                )
                iteration_results_list.append(results_df)

            iteration_subjects_df = pd.concat(iteration_subjects_list)
            iteration_subjects_df.to_csv(
                path.join(iteration_dir, "subjects.tsv"),
                index=False,
                sep="\t",
                encoding="utf-8",
            )
            all_subjects_list.append(iteration_subjects_df)

            iteration_results_df = pd.concat(iteration_results_list)
            iteration_results_df.to_csv(
                path.join(iteration_dir, "results.tsv"),
                index=False,
                sep="\t",
                encoding="utf-8",
            )

            mean_results_df = pd.DataFrame(
                iteration_results_df.apply(np.nanmean).to_dict(),
                columns=iteration_results_df.columns,
                index=[
                    0,
                ],
            )
            mean_results_df.to_csv(
                path.join(iteration_dir, "mean_results.tsv"),
                index=False,
                sep="\t",
                encoding="utf-8",
            )
            all_results_list.append(mean_results_df)

        all_subjects_df = pd.concat(all_subjects_list)
        all_subjects_df.to_csv(
            path.join(output_dir, "subjects.tsv"),
            index=False,
            sep="\t",
            encoding="utf-8",
        )

        all_results_df = pd.concat(all_results_list)
        all_results_df.to_csv(
            path.join(output_dir, "results.tsv"),
            index=False,
            sep="\t",
            encoding="utf-8",
        )

        mean_results_df = pd.DataFrame(
            all_results_df.apply(np.nanmean).to_dict(),
            columns=all_results_df.columns,
            index=[
                0,
            ],
        )
        mean_results_df.to_csv(
            path.join(output_dir, "mean_results.tsv"),
            index=False,
            sep="\t",
            encoding="utf-8",
        )

        print("Mean results of the classification:")
        print(
            "Balanced accuracy: %s"
            % (mean_results_df["balanced_accuracy"].to_string(index=False))
        )
