import itertools
import json
from multiprocessing.pool import ThreadPool
from os import path

import numpy as np
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.multiclass import OneVsOneClassifier, OneVsRestClassifier
from sklearn.svm import SVC
from xgboost import XGBClassifier

import clinica.pipelines.machine_learning.ml_utils as utils
from clinica.pipelines.machine_learning import base


class DualSVMAlgorithm(base.MLAlgorithm):
    def _launch_svc(self, kernel_train, x_test, y_train, y_test, c):

        if self._algorithm_params["balanced"]:
            svc = SVC(
                C=c,
                kernel="precomputed",
                probability=True,
                tol=1e-6,
                class_weight="balanced",
            )
        else:
            svc = SVC(C=c, kernel="precomputed", probability=True, tol=1e-6)

        svc.fit(kernel_train, y_train)
        y_hat_train = svc.predict(kernel_train)
        y_hat = svc.predict(x_test)
        proba_test = svc.predict_proba(x_test)[:, 1]
        auc = roc_auc_score(y_test, proba_test)

        return svc, y_hat, auc, y_hat_train

    def _grid_search(self, kernel_train, x_test, y_train, y_test, c):

        _, y_hat, _, _ = self._launch_svc(kernel_train, x_test, y_train, y_test, c)
        res = utils.evaluate_prediction(y_test, y_hat)

        return res["balanced_accuracy"]

    def _select_best_parameter(self, async_result):

        c_values = []
        accuracies = []
        for fold in async_result.keys():
            best_c = -1
            best_acc = -1

            for c, async_acc in async_result[fold].items():

                acc = async_acc.get()
                if acc > best_acc:
                    best_c = c
                    best_acc = acc
            c_values.append(best_c)
            accuracies.append(best_acc)

        best_acc = np.mean(accuracies)
        best_c = np.power(10, np.mean(np.log10(c_values)))

        return {"c": best_c, "balanced_accuracy": best_acc}

    def evaluate(self, train_index, test_index):

        inner_pool = ThreadPool(self._algorithm_params["n_threads"])
        async_result = {}
        for i in range(self._algorithm_params["grid_search_folds"]):
            async_result[i] = {}

        outer_kernel = self._kernel[train_index, :][:, train_index]
        y_train = self._y[train_index]

        skf = StratifiedKFold(
            n_splits=self._algorithm_params["grid_search_folds"], shuffle=True
        )
        inner_cv = list(skf.split(np.zeros(len(y_train)), y_train))

        for i in range(len(inner_cv)):
            inner_train_index, inner_test_index = inner_cv[i]

            inner_kernel = outer_kernel[inner_train_index, :][:, inner_train_index]
            x_test_inner = outer_kernel[inner_test_index, :][:, inner_train_index]
            y_train_inner, y_test_inner = (
                y_train[inner_train_index],
                y_train[inner_test_index],
            )

            for c in self._algorithm_params["c_range"]:
                async_result[i][c] = inner_pool.apply_async(
                    self._grid_search,
                    (inner_kernel, x_test_inner, y_train_inner, y_test_inner, c),
                )
        inner_pool.close()
        inner_pool.join()

        best_parameter = self._select_best_parameter(async_result)
        x_test = self._kernel[test_index, :][:, train_index]
        y_train, y_test = self._y[train_index], self._y[test_index]

        _, y_hat, auc, y_hat_train = self._launch_svc(
            outer_kernel, x_test, y_train, y_test, best_parameter["c"]
        )

        result = dict()
        result["best_parameter"] = best_parameter
        result["evaluation"] = utils.evaluate_prediction(y_test, y_hat)
        result["evaluation_train"] = utils.evaluate_prediction(y_train, y_hat_train)
        result["y_hat"] = y_hat
        result["y_hat_train"] = y_hat_train
        result["y"] = y_test
        result["y_train"] = y_train
        result["y_index"] = test_index
        result["x_index"] = train_index
        result["auc"] = auc

        return result

    def apply_best_parameters(self, results_list):

        best_c_list = []
        bal_acc_list = []

        for result in results_list:
            best_c_list.append(result["best_parameter"]["c"])
            bal_acc_list.append(result["best_parameter"]["balanced_accuracy"])

        # 10^(mean of log10 of best Cs of each fold) is selected
        best_c = np.power(10, np.mean(np.log10(best_c_list)))
        # Mean balanced accuracy
        mean_bal_acc = np.mean(bal_acc_list)

        if self._algorithm_params["balanced"]:
            svc = SVC(
                C=best_c,
                kernel="precomputed",
                probability=True,
                tol=1e-6,
                class_weight="balanced",
            )
        else:
            svc = SVC(C=best_c, kernel="precomputed", probability=True, tol=1e-6)

        svc.fit(self._kernel, self._y)

        return svc, {"c": best_c, "balanced_accuracy": mean_bal_acc}

    def save_classifier(self, classifier, output_dir):

        np.savetxt(
            path.join(output_dir, "dual_coefficients.txt"), classifier.dual_coef_
        )
        np.savetxt(
            path.join(output_dir, "support_vectors_indices.txt"), classifier.support_
        )
        np.savetxt(path.join(output_dir, "intersect.txt"), classifier.intercept_)

    def save_weights(self, classifier, x, output_dir):

        dual_coefficients = classifier.dual_coef_
        sv_indices = classifier.support_

        weighted_sv = dual_coefficients.transpose() * x[sv_indices]
        weights = np.sum(weighted_sv, 0)

        np.savetxt(path.join(output_dir, "weights.txt"), weights)

        return weights

    def save_parameters(self, parameters_dict, output_dir):
        with open(path.join(output_dir, "best_parameters.json"), "w") as f:
            json.dump(parameters_dict, f)

    @staticmethod
    def uses_kernel():
        return True

    @staticmethod
    def get_default_parameters():

        parameters_dict = {
            "balanced": True,
            "grid_search_folds": 10,
            "c_range": np.logspace(-6, 2, 17),
            "n_threads": 15,
        }

        return parameters_dict


class LogisticReg(base.MLAlgorithm):
    def _launch_logistic_reg(self, x_train, x_test, y_train, y_test, c):

        if self._algorithm_params["balanced"]:
            classifier = LogisticRegression(
                penalty=self._algorithm_params["penalty"],
                tol=1e-6,
                C=c,
                class_weight="balanced",
            )
        else:
            classifier = LogisticRegression(
                penalty=self._algorithm_params["penalty"], tol=1e-6, C=c
            )

        classifier.fit(x_train, y_train)
        y_hat_train = classifier.predict(x_train)
        y_hat = classifier.predict(x_test)
        proba_test = classifier.predict_proba(x_test)[:, 1]
        auc = roc_auc_score(y_test, proba_test)

        return classifier, y_hat, auc, y_hat_train

    def _grid_search(self, x_train, x_test, y_train, y_test, c):

        _, y_hat, _, _ = self._launch_logistic_reg(x_train, x_test, y_train, y_test, c)
        res = utils.evaluate_prediction(y_test, y_hat)

        return res["balanced_accuracy"]

    def _select_best_parameter(self, async_result):

        c_values = []
        accuracies = []
        for fold in async_result.keys():
            best_c = -1
            best_acc = -1

            for c, async_acc in async_result[fold].items():

                acc = async_acc.get()
                if acc > best_acc:
                    best_c = c
                    best_acc = acc
            c_values.append(best_c)
            accuracies.append(best_acc)

        best_acc = np.mean(accuracies)
        best_c = np.power(10, np.mean(np.log10(c_values)))

        return {"c": best_c, "balanced_accuracy": best_acc}

    def evaluate(self, train_index, test_index):

        inner_pool = ThreadPool(self._algorithm_params["n_threads"])
        async_result = {}
        for i in range(self._algorithm_params["grid_search_folds"]):
            async_result[i] = {}

        x_train = self._x[train_index]
        y_train = self._y[train_index]

        skf = StratifiedKFold(
            n_splits=self._algorithm_params["grid_search_folds"], shuffle=True
        )
        inner_cv = list(skf.split(np.zeros(len(y_train)), y_train))

        for i in range(len(inner_cv)):
            inner_train_index, inner_test_index = inner_cv[i]

            x_train_inner = x_train[inner_train_index]
            x_test_inner = x_train[inner_test_index]
            y_train_inner = y_train[inner_train_index]
            y_test_inner = y_train[inner_test_index]

            for c in self._algorithm_params["c_range"]:
                async_result[i][c] = inner_pool.apply_async(
                    self._grid_search,
                    (x_train_inner, x_test_inner, y_train_inner, y_test_inner, c),
                )
        inner_pool.close()
        inner_pool.join()

        best_parameter = self._select_best_parameter(async_result)
        x_test = self._x[test_index]
        y_test = self._y[test_index]

        _, y_hat, auc, y_hat_train = self._launch_logistic_reg(
            x_train, x_test, y_train, y_test, best_parameter["c"]
        )

        result = dict()
        result["best_parameter"] = best_parameter
        result["evaluation"] = utils.evaluate_prediction(y_test, y_hat)
        result["evaluation_train"] = utils.evaluate_prediction(y_train, y_hat_train)
        result["y_hat"] = y_hat
        result["y_hat_train"] = y_hat_train
        result["y"] = y_test
        result["y_train"] = y_train
        result["y_index"] = test_index
        result["x_index"] = train_index
        result["auc"] = auc

        return result

    def apply_best_parameters(self, results_list):

        best_c_list = []
        bal_acc_list = []

        for result in results_list:
            best_c_list.append(result["best_parameter"]["c"])
            bal_acc_list.append(result["best_parameter"]["balanced_accuracy"])

        # 10^(mean of log10 of best Cs of each fold) is selected
        best_c = np.power(10, np.mean(np.log10(best_c_list)))
        # Mean balanced accuracy
        mean_bal_acc = np.mean(bal_acc_list)

        if self._algorithm_params["balanced"]:
            classifier = LogisticRegression(
                C=best_c,
                penalty=self._algorithm_params["penalty"],
                tol=1e-6,
                class_weight="balanced",
            )
        else:
            classifier = LogisticRegression(
                C=best_c, penalty=self._algorithm_params["penalty"], tol=1e-6
            )

        classifier.fit(self._x, self._y)

        return classifier, {"c": best_c, "balanced_accuracy": mean_bal_acc}

    def save_classifier(self, classifier, output_dir):

        np.savetxt(path.join(output_dir, "weights.txt"), classifier.coef_.transpose())
        np.savetxt(path.join(output_dir, "intercept.txt"), classifier.intercept_)

    def save_weights(self, classifier, x, output_dir):

        np.savetxt(path.join(output_dir, "weights.txt"), classifier.coef_.transpose())
        return classifier.coef_.transpose()

    def save_parameters(self, parameters_dict, output_dir):

        with open(path.join(output_dir, "best_parameters.json"), "w") as f:
            json.dump(parameters_dict, f)

    @staticmethod
    def _centered_normalised_data(features):
        std = np.std(features, axis=0)
        std[np.where(std == 0)[0]] = 1.0
        mean = np.mean(features, axis=0)
        features_bis = (features - mean) / std
        return features_bis, mean, std

    @staticmethod
    def uses_kernel():
        return False

    @staticmethod
    def get_default_parameters():
        parameters_dict = {
            "penalty": "l2",
            "balanced": False,
            "grid_search_folds": 10,
            "c_range": np.logspace(-6, 2, 17),
            "n_threads": 15,
        }

        return parameters_dict


class RandomForest(base.MLAlgorithm):
    def _launch_random_forest(
        self,
        x_train,
        x_test,
        y_train,
        y_test,
        n_estimators,
        max_depth,
        min_samples_split,
        max_features,
    ):

        if self._algorithm_params["balanced"]:
            classifier = RandomForestClassifier(
                n_estimators=n_estimators,
                max_depth=max_depth,
                min_samples_split=min_samples_split,
                max_features=max_features,
                class_weight="balanced",
                n_jobs=self._algorithm_params["n_threads"],
            )
        else:
            classifier = RandomForestClassifier(
                n_estimators=n_estimators,
                max_depth=max_depth,
                min_samples_split=min_samples_split,
                max_features=max_features,
                n_jobs=self._algorithm_params["n_threads"],
            )

        classifier.fit(x_train, y_train)
        y_hat_train = classifier.predict(x_train)
        y_hat = classifier.predict(x_test)
        proba_test = classifier.predict_proba(x_test)[:, 1]
        auc = roc_auc_score(y_test, proba_test)

        return classifier, y_hat, auc, y_hat_train

    def _grid_search(
        self,
        x_train,
        x_test,
        y_train,
        y_test,
        n_estimators,
        max_depth,
        min_samples_split,
        max_features,
    ):

        _, y_hat, _, _ = self._launch_random_forest(
            x_train,
            x_test,
            y_train,
            y_test,
            n_estimators,
            max_depth,
            min_samples_split,
            max_features,
        )
        res = utils.evaluate_prediction(y_test, y_hat)

        return res["balanced_accuracy"]

    def _select_best_parameter(self, async_result):

        params_list = []
        accuracies = []

        all_params_acc = []

        for fold in async_result.keys():
            best_params = None
            best_acc = -1

            for params, async_acc in async_result[fold].items():
                acc = async_acc.get()
                if acc > best_acc:
                    best_params = params
                    best_acc = acc

                all_params_acc.append(
                    pd.DataFrame(
                        {
                            "n_estimators": params[0],
                            "max_depth": params[1],
                            "min_samples_split": params[2],
                            "max_features": params[3],
                            "balanced_accuracy": acc,
                        },
                        index=[
                            "i",
                        ],
                    )
                )

            params_list.append(best_params)
            accuracies.append(best_acc)

        best_acc = np.mean(accuracies)
        best_n_estimators = int(round(np.mean([x[0] for x in params_list])))
        best_max_depth = int(
            round(np.mean([x[1] if x[1] is not None else 50 for x in params_list]))
        )
        best_min_samples_split = int(round(np.mean([x[2] for x in params_list])))

        def max_feature_to_float(m):
            if type(m) is float:
                return m
            if type(m) is int:
                return float(m) / float(self._x.shape[1])
            if m == "auto" or m == "sqrt":
                return np.sqrt(self._x.shape[1]) / float(self._x.shape[1])
            if m == "log2":
                return np.log2(self._x.shape[1]) / float(self._x.shape[1])
            raise ValueError("Not valid value for max_feature: %s" % m)

        float_max_feat = [max_feature_to_float(x[3]) for x in params_list]
        best_max_features = np.mean(float_max_feat)

        return {
            "n_estimators": best_n_estimators,
            "max_depth": best_max_depth,
            "min_samples_split": best_min_samples_split,
            "max_features": best_max_features,
            "balanced_accuracy": best_acc,
        }

    def evaluate(self, train_index, test_index):

        inner_pool = ThreadPool(self._algorithm_params["n_threads"])
        async_result = {}
        for i in range(self._algorithm_params["grid_search_folds"]):
            async_result[i] = {}

        x_train = self._x[train_index]
        y_train = self._y[train_index]

        skf = StratifiedKFold(
            n_splits=self._algorithm_params["grid_search_folds"], shuffle=True
        )
        inner_cv = list(skf.split(np.zeros(len(y_train)), y_train))

        parameters_combinations = list(
            itertools.product(
                self._algorithm_params["n_estimators_range"],
                self._algorithm_params["max_depth_range"],
                self._algorithm_params["min_samples_split_range"],
                self._algorithm_params["max_features_range"],
            )
        )

        for i in range(len(inner_cv)):
            inner_train_index, inner_test_index = inner_cv[i]

            x_train_inner = x_train[inner_train_index]
            x_test_inner = x_train[inner_test_index]
            y_train_inner = y_train[inner_train_index]
            y_test_inner = y_train[inner_test_index]

            for parameters in parameters_combinations:
                async_result[i][parameters] = inner_pool.apply_async(
                    self._grid_search,
                    (
                        x_train_inner,
                        x_test_inner,
                        y_train_inner,
                        y_test_inner,
                        parameters[0],
                        parameters[1],
                        parameters[2],
                        parameters[3],
                    ),
                )
        inner_pool.close()
        inner_pool.join()
        best_parameter = self._select_best_parameter(async_result)
        x_test = self._x[test_index]
        y_test = self._y[test_index]

        _, y_hat, auc, y_hat_train = self._launch_random_forest(
            x_train,
            x_test,
            y_train,
            y_test,
            best_parameter["n_estimators"],
            best_parameter["max_depth"],
            best_parameter["min_samples_split"],
            best_parameter["max_features"],
        )

        result = dict()
        result["best_parameter"] = best_parameter
        result["evaluation"] = utils.evaluate_prediction(y_test, y_hat)
        result["evaluation_train"] = utils.evaluate_prediction(y_train, y_hat_train)
        result["y_hat"] = y_hat
        result["y_hat_train"] = y_hat_train
        result["y"] = y_test
        result["y_train"] = y_train
        result["y_index"] = test_index
        result["x_index"] = train_index
        result["auc"] = auc

        return result

    def evaluate_no_cv(self, train_index, test_index):

        x_train = self._x[train_index]
        y_train = self._y[train_index]
        x_test = self._x[test_index]
        y_test = self._y[test_index]

        best_parameter = dict()
        best_parameter["n_estimators"] = self._algorithm_params["n_estimators_range"]
        best_parameter["max_depth"] = self._algorithm_params["max_depth_range"]
        best_parameter["min_samples_split"] = self._algorithm_params[
            "min_samples_split_range"
        ]
        best_parameter["max_features"] = self._algorithm_params["max_features_range"]

        _, y_hat, auc, y_hat_train = self._launch_random_forest(
            x_train,
            x_test,
            y_train,
            y_test,
            self._algorithm_params["n_estimators_range"],
            self._algorithm_params["max_depth_range"],
            self._algorithm_params["min_samples_split_range"],
            self._algorithm_params["max_features_range"],
        )
        result = dict()
        result["best_parameter"] = best_parameter
        result["evaluation"] = utils.evaluate_prediction(y_test, y_hat)
        best_parameter["balanced_accuracy"] = result["evaluation"]["balanced_accuracy"]
        result["evaluation_train"] = utils.evaluate_prediction(y_train, y_hat_train)
        result["y_hat"] = y_hat
        result["y_hat_train"] = y_hat_train
        result["y"] = y_test
        result["y_train"] = y_train
        result["y_index"] = test_index
        result["x_index"] = train_index
        result["auc"] = auc

        return result

    def apply_best_parameters(self, results_list):

        mean_bal_acc = np.mean(
            [result["best_parameter"]["balanced_accuracy"] for result in results_list]
        )
        best_n_estimators = int(
            round(
                np.mean(
                    [
                        result["best_parameter"]["n_estimators"]
                        for result in results_list
                    ]
                )
            )
        )
        best_max_depth = int(
            round(
                np.mean(
                    [
                        result["best_parameter"]["max_depth"]
                        if result["best_parameter"]["max_depth"] is not None
                        else 50
                        for result in results_list
                    ]
                )
            )
        )
        best_min_samples_split = int(
            round(
                np.mean(
                    [
                        result["best_parameter"]["min_samples_split"]
                        for result in results_list
                    ]
                )
            )
        )

        max_feat = []
        n_features = self._x.shape[1]
        for result in results_list:
            result_feat = result["best_parameter"]["max_features"]

            if not result_feat:
                max_features = 1.0
            elif result_feat in ["auto", "sqrt"]:
                max_features = np.sqrt(n_features) / n_features
            elif result_feat == "log2":
                max_features = np.log2(n_features) / n_features
            elif isinstance(result_feat, int):
                max_features = float(result_feat) / n_features
            elif isinstance(result_feat, float):
                max_features = result_feat
            else:
                raise ("Unknown max_features type")

            max_feat.append(max_features)
        best_max_features = np.mean(max_feat)

        if self._algorithm_params["balanced"]:
            classifier = RandomForestClassifier(
                n_estimators=best_n_estimators,
                max_depth=best_max_depth,
                min_samples_split=best_min_samples_split,
                max_features=best_max_features,
                class_weight="balanced",
                n_jobs=self._algorithm_params["n_threads"],
            )
        else:
            classifier = RandomForestClassifier(
                n_estimators=best_n_estimators,
                max_depth=best_max_depth,
                min_samples_split=best_min_samples_split,
                max_features=best_max_features,
                n_jobs=self._algorithm_params["n_threads"],
            )

        classifier.fit(self._x, self._y)

        return classifier, {
            "n_estimators": best_n_estimators,
            "max_depth": best_max_depth,
            "min_samples_split": best_min_samples_split,
            "max_features": best_max_features,
            "balanced_accuracy": mean_bal_acc,
        }

    def save_classifier(self, classifier, output_dir):
        np.savetxt(
            path.join(output_dir, "feature_importances.txt"),
            classifier.feature_importances_,
        )
        # print classifier.estimators_
        # np.savetxt(path.join(output_dir, 'estimators.txt'), str(classifier.estimators_))

    def save_weights(self, classifier, x, output_dir):

        np.savetxt(
            path.join(output_dir, "weights.txt"), classifier.feature_importances_
        )
        return classifier.feature_importances_

    def save_parameters(self, parameters_dict, output_dir):

        with open(path.join(output_dir, "best_parameters.json"), "w") as f:
            json.dump(parameters_dict, f)

    @staticmethod
    def uses_kernel():
        return False

    @staticmethod
    def get_default_parameters():

        parameters_dict = {
            "balanced": False,
            "grid_search_folds": 10,
            "n_estimators_range": (10, 25, 50, 100, 150, 200, 500),
            "max_depth_range": (None, 6, 8, 10, 12),
            "min_samples_split_range": (2, 4, 6, 8),
            "max_features_range": ("auto", 0.1, 0.2, 0.3, 0.4, 0.5),
            "n_threads": 15,
        }

        return parameters_dict


class XGBoost(base.MLAlgorithm):
    def _launch_xgboost(
        self,
        x_train,
        x_test,
        y_train,
        y_test,
        max_depth,
        learning_rate,
        n_estimators,
        colsample_bytree,
    ):

        if self._algorithm_params["balanced"]:
            # set scale_pos_weight
            # http://xgboost.readthedocs.io/en/latest//how_to/param_tuning.html
            scale_pos_weight = float(len(self._y - sum(self._y)) / sum(self._y))
            classifier = XGBClassifier(
                max_depth=max_depth,
                learning_rate=learning_rate,
                n_estimators=n_estimators,
                n_jobs=self._algorithm_params["n_threads"],
                colsample_bytree=colsample_bytree,
                reg_alpha=self._algorithm_params["reg_alpha"],
                reg_lambda=self._algorithm_params["reg_lambda"],
                scale_pos_weight=scale_pos_weight,
            )
        else:
            classifier = XGBClassifier(
                max_depth=max_depth,
                learning_rate=learning_rate,
                n_estimators=n_estimators,
                n_jobs=self._algorithm_params["n_threads"],
                colsample_bytree=colsample_bytree,
                reg_alpha=self._algorithm_params["reg_alpha"],
                reg_lambda=self._algorithm_params["reg_lambda"],
            )

        classifier.fit(x_train, y_train)
        y_hat_train = classifier.predict(x_train)
        y_hat = classifier.predict(x_test)
        proba_test = classifier.predict_proba(x_test)[:, 1]
        auc = roc_auc_score(y_test, proba_test)

        return classifier, y_hat, auc, y_hat_train

    def _grid_search(
        self,
        x_train,
        x_test,
        y_train,
        y_test,
        max_depth,
        learning_rate,
        n_estimators,
        colsample_bytree,
    ):

        _, y_hat, _, _ = self._launch_xgboost(
            x_train,
            x_test,
            y_train,
            y_test,
            max_depth,
            learning_rate,
            n_estimators,
            colsample_bytree,
        )
        res = utils.evaluate_prediction(y_test, y_hat)

        return res["balanced_accuracy"]

    def _select_best_parameter(self, async_result):

        params_list = []
        accuracies = []

        all_params_acc = []

        for fold in async_result.keys():
            best_params = None
            best_acc = -1

            for params, async_acc in async_result[fold].items():
                acc = async_acc.get()
                if acc > best_acc:
                    best_params = params
                    best_acc = acc

                all_params_acc.append(
                    pd.DataFrame(
                        {
                            "max_depth": params[0],
                            "learning_rate": params[1],
                            "n_estimators": params[2],
                            "colsample_bytree": params[3],
                            "balanced_accuracy": acc,
                        },
                        index=[
                            "i",
                        ],
                    )
                )

            params_list.append(best_params)
            accuracies.append(best_acc)

        best_acc = np.mean(accuracies)

        best_max_depth = int(round(np.mean([x[0] for x in params_list])))
        best_learning_rate = np.mean([x[1] for x in params_list])
        best_n_estimators = int(round(np.mean([x[2] for x in params_list])))
        best_colsample_bytree = np.mean([x[3] for x in params_list])

        return {
            "max_depth": best_max_depth,
            "learning_rate": best_learning_rate,
            "n_estimators": best_n_estimators,
            "colsample_bytree": best_colsample_bytree,
            "balanced_accuracy": best_acc,
        }

    def evaluate(self, train_index, test_index):

        inner_pool = ThreadPool(self._algorithm_params["n_threads"])
        async_result = {}
        for i in range(self._algorithm_params["grid_search_folds"]):
            async_result[i] = {}

        x_train = self._x[train_index]
        y_train = self._y[train_index]

        skf = StratifiedKFold(
            n_splits=self._algorithm_params["grid_search_folds"], shuffle=True
        )
        inner_cv = list(skf.split(np.zeros(len(y_train)), y_train))

        parameters_combinations = list(
            itertools.product(
                self._algorithm_params["max_depth_range"],
                self._algorithm_params["learning_rate_range"],
                self._algorithm_params["n_estimators_range"],
                self._algorithm_params["colsample_bytree_range"],
            )
        )

        for i in range(len(inner_cv)):
            inner_train_index, inner_test_index = inner_cv[i]

            x_train_inner = x_train[inner_train_index]
            x_test_inner = x_train[inner_test_index]
            y_train_inner = y_train[inner_train_index]
            y_test_inner = y_train[inner_test_index]

            for parameters in parameters_combinations:
                async_result[i][parameters] = inner_pool.apply_async(
                    self._grid_search,
                    (
                        x_train_inner,
                        x_test_inner,
                        y_train_inner,
                        y_test_inner,
                        parameters[0],
                        parameters[1],
                        parameters[2],
                        parameters[3],
                    ),
                )
        inner_pool.close()
        inner_pool.join()
        best_parameter = self._select_best_parameter(async_result)
        x_test = self._x[test_index]
        y_test = self._y[test_index]

        _, y_hat, auc, y_hat_train = self._launch_xgboost(
            x_train,
            x_test,
            y_train,
            y_test,
            best_parameter["max_depth"],
            best_parameter["learning_rate"],
            best_parameter["n_estimators"],
            best_parameter["colsample_bytree"],
        )

        result = dict()
        result["best_parameter"] = best_parameter
        result["evaluation"] = utils.evaluate_prediction(y_test, y_hat)
        result["evaluation_train"] = utils.evaluate_prediction(y_train, y_hat_train)
        result["y_hat"] = y_hat
        result["y_hat_train"] = y_hat_train
        result["y"] = y_test
        result["y_train"] = y_train
        result["y_index"] = test_index
        result["x_index"] = train_index
        result["auc"] = auc

        return result

    def evaluate_no_cv(self, train_index, test_index):

        x_train = self._x[train_index]
        y_train = self._y[train_index]
        x_test = self._x[test_index]
        y_test = self._y[test_index]

        best_parameter = dict()
        best_parameter["max_depth"] = self._algorithm_params["max_depth_range"]
        best_parameter["learning_rate"] = self._algorithm_params["learning_rate_range"]
        best_parameter["n_estimators"] = self._algorithm_params["n_estimators_range"]
        best_parameter["colsample_bytree"] = self._algorithm_params[
            "colsample_bytree_range"
        ]

        _, y_hat, auc, y_hat_train = self._launch_xgboost(
            x_train,
            x_test,
            y_train,
            y_test,
            self._algorithm_params["max_depth_range"],
            self._algorithm_params["learning_rate_range"],
            self._algorithm_params["n_estimators_range"],
            self._algorithm_params["colsample_bytree_range"],
        )
        result = dict()
        result["best_parameter"] = best_parameter
        result["evaluation"] = utils.evaluate_prediction(y_test, y_hat)
        best_parameter["balanced_accuracy"] = result["evaluation"]["balanced_accuracy"]
        result["evaluation_train"] = utils.evaluate_prediction(y_train, y_hat_train)
        result["y_hat"] = y_hat
        result["y_hat_train"] = y_hat_train
        result["y"] = y_test
        result["y_train"] = y_train
        result["y_index"] = test_index
        result["x_index"] = train_index
        result["auc"] = auc

        return result

    def apply_best_parameters(self, results_list):

        mean_bal_acc = np.mean(
            [result["best_parameter"]["balanced_accuracy"] for result in results_list]
        )

        best_max_depth = int(
            round(
                np.mean(
                    [result["best_parameter"]["max_depth"] for result in results_list]
                )
            )
        )
        best_learning_rate = np.mean(
            [result["best_parameter"]["learning_rate"] for result in results_list]
        )
        best_n_estimators = int(
            round(
                np.mean(
                    [
                        result["best_parameter"]["n_estimators"]
                        for result in results_list
                    ]
                )
            )
        )
        best_colsample_bytree = np.mean(
            [result["best_parameter"]["colsample_bytree"] for result in results_list]
        )

        if self._algorithm_params["balanced"]:
            scale_pos_weight = float(len(self._y - sum(self._y)) / sum(self._y))

            classifier = XGBClassifier(
                max_depth=best_max_depth,
                learning_rate=best_learning_rate,
                n_estimators=best_n_estimators,
                n_jobs=self._algorithm_params["n_threads"],
                colsample_bytree=best_colsample_bytree,
                reg_alpha=self._algorithm_params["reg_alpha"],
                reg_lambda=self._algorithm_params["reg_lambda"],
                scale_pos_weight=scale_pos_weight,
            )
        else:
            classifier = XGBClassifier(
                max_depth=best_max_depth,
                learning_rate=best_learning_rate,
                n_estimators=best_n_estimators,
                n_jobs=self._algorithm_params["n_threads"],
                colsample_bytree=best_colsample_bytree,
                reg_alpha=self._algorithm_params["reg_alpha"],
                reg_lambda=self._algorithm_params["reg_lambda"],
            )

        classifier.fit(self._x, self._y)

        return classifier, {
            "max_depth": best_max_depth,
            "learning_rate": best_learning_rate,
            "n_estimators": best_n_estimators,
            "colsample_bytree": best_colsample_bytree,
            "balanced_accuracy": mean_bal_acc,
        }

    def save_classifier(self, classifier, output_dir):
        np.savetxt(
            path.join(output_dir, "feature_importances.txt"),
            classifier.feature_importances_,
        )
        # print classifier.estimators_
        # np.savetxt(path.join(output_dir, 'estimators.txt'), str(classifier.estimators_))

    def save_weights(self, classifier, x, output_dir):

        np.savetxt(
            path.join(output_dir, "weights.txt"), classifier.feature_importances_
        )
        return classifier.feature_importances_

    def save_parameters(self, parameters_dict, output_dir):

        with open(path.join(output_dir, "best_parameters.json"), "w") as f:
            json.dump(parameters_dict, f)

    @staticmethod
    def uses_kernel():
        return False

    @staticmethod
    def get_default_parameters():
        parameters_dict = {
            "balanced": False,
            "grid_search_folds": 10,
            "max_depth_range": (0, 6),
            "learning_rate_range": (0.1, 0.3),
            "n_estimators_range": (100, 200),
            "colsample_bytree_range": (0.5, 1),
            "reg_alpha": 0,
            "reg_lambda": 1,
            "n_threads": 15,
        }

        return parameters_dict


class OneVsOneSVM(base.MLAlgorithm):
    def _launch_svc(self, kernel_train, x_test, y_train, y_test, c):

        if self._algorithm_params["balanced"]:
            svc = OneVsOneClassifier(
                SVC(
                    C=c,
                    kernel="precomputed",
                    probability=True,
                    tol=1e-6,
                    class_weight="balanced",
                )
            )
        else:
            svc = OneVsOneClassifier(
                SVC(C=c, kernel="precomputed", probability=True, tol=1e-6)
            )

        svc.fit(kernel_train, y_train)
        y_hat_train = svc.predict(kernel_train)
        y_hat = svc.predict(x_test)
        proba_test = svc.predict_proba(x_test)[:, 1]

        return svc, y_hat, y_hat_train

    def _grid_search(self, kernel_train, x_test, y_train, y_test, c):
        # y_hat is the value predicted

        _, y_hat, _ = self._launch_svc(kernel_train, x_test, y_train, y_test, c)
        res = utils.evaluate_prediction_multiclass(y_test, y_hat)

        return res["balanced_accuracy"]

    def _select_best_parameter(self, async_result):

        c_values = []
        accuracies = []
        for fold in async_result.keys():
            best_c = -1
            best_acc = -1

            for c, async_acc in async_result[fold].items():

                acc = async_acc.get()
                if acc > best_acc:
                    best_c = c
                    best_acc = acc
            c_values.append(best_c)
            accuracies.append(best_acc)

        best_acc = np.mean(accuracies)
        best_c = np.power(10, np.mean(np.log10(c_values)))

        return {"c": best_c, "balanced_accuracy": best_acc}

    def evaluate(self, train_index, test_index):

        inner_pool = ThreadPool(self._algorithm_params["n_threads"])
        async_result = {}
        for i in range(self._algorithm_params["grid_search_folds"]):
            async_result[i] = {}

        outer_kernel = self._kernel[train_index, :][:, train_index]
        y_train = self._y[train_index]

        skf = StratifiedKFold(
            n_splits=self._algorithm_params["grid_search_folds"], shuffle=True
        )
        inner_cv = list(skf.split(np.zeros(len(y_train)), y_train))

        for i in range(len(inner_cv)):
            inner_train_index, inner_test_index = inner_cv[i]

            inner_kernel = outer_kernel[inner_train_index, :][:, inner_train_index]
            x_test_inner = outer_kernel[inner_test_index, :][:, inner_train_index]
            y_train_inner, y_test_inner = (
                y_train[inner_train_index],
                y_train[inner_test_index],
            )

            for c in self._algorithm_params["c_range"]:
                async_result[i][c] = inner_pool.apply_async(
                    self._grid_search,
                    (inner_kernel, x_test_inner, y_train_inner, y_test_inner, c),
                )
        inner_pool.close()
        inner_pool.join()

        best_parameter = self._select_best_parameter(async_result)
        x_test = self._kernel[test_index, :][:, train_index]
        y_train, y_test = self._y[train_index], self._y[test_index]

        _, y_hat, y_hat_train = self._launch_svc(
            outer_kernel, x_test, y_train, y_test, best_parameter["c"]
        )

        result = dict()
        result["best_parameter"] = best_parameter
        result["evaluation"] = utils.evaluate_prediction_multiclass(y_test, y_hat)
        result["evaluation_train"] = utils.evaluate_prediction_multiclass(
            y_train, y_hat_train
        )
        result["y_hat"] = y_hat
        result["y_hat_train"] = y_hat_train
        result["y"] = y_test
        result["y_train"] = y_train
        result["y_index"] = test_index
        result["x_index"] = train_index

        return result

    def apply_best_parameters(self, results_list):

        best_c_list = []
        bal_acc_list = []

        for result in results_list:
            best_c_list.append(result["best_parameter"]["c"])
            bal_acc_list.append(result["best_parameter"]["balanced_accuracy"])

        # 10^(mean of log10 of best Cs of each fold) is selected
        best_c = np.power(10, np.mean(np.log10(best_c_list)))
        # Mean balanced accuracy
        mean_bal_acc = np.mean(bal_acc_list)

        if self._algorithm_params["balanced"]:
            svc = OneVsOneClassifier(
                SVC(
                    C=best_c,
                    kernel="precomputed",
                    probability=True,
                    tol=1e-6,
                    class_weight="balanced",
                )
            )
        else:
            svc = OneVsOneClassifier(
                SVC(C=best_c, kernel="precomputed", probability=True, tol=1e-6)
            )

        svc.fit(self._kernel, self._y)

        return svc, {"c": best_c, "balanced_accuracy": mean_bal_acc}

    def save_classifier(self, classifier, output_dir):

        np.savetxt(
            path.join(output_dir, "support_vectors_indices.txt"), classifier.support_
        )
        np.savetxt(path.join(output_dir, "intersect.txt"), classifier.intercept_)

    def save_weights(self, classifier, x, output_dir):

        dual_coefficients = classifier.dual_coef_
        sv_indices = classifier.support_

        weighted_sv = dual_coefficients.transpose() * x[sv_indices]
        weights = np.sum(weighted_sv, 0)

        np.savetxt(path.join(output_dir, "weights.txt"), weights)

        return weights

    def save_parameters(self, parameters_dict, output_dir):
        with open(path.join(output_dir, "best_parameters.json"), "w") as f:
            json.dump(parameters_dict, f)

    @staticmethod
    def uses_kernel():
        return True

    @staticmethod
    def get_default_parameters():

        parameters_dict = {
            "balanced": True,
            "grid_search_folds": 10,
            "c_range": np.logspace(-6, 2, 17),
            "n_threads": 15,
        }

        return parameters_dict


class OneVsRestSVM(base.MLAlgorithm):
    def _launch_svc(self, kernel_train, x_test, y_train, y_test, c):

        if self._algorithm_params["balanced"]:
            svc = OneVsRestClassifier(
                SVC(
                    C=c,
                    kernel="precomputed",
                    probability=True,
                    tol=1e-6,
                    class_weight="balanced",
                )
            )
        else:
            svc = OneVsRestClassifier(
                SVC(C=c, kernel="precomputed", probability=True, tol=1e-6)
            )

        svc.fit(kernel_train, y_train)
        y_hat_train = svc.predict(kernel_train)
        y_hat = svc.predict(x_test)
        proba_test = svc.predict_proba(x_test)[:, 1]

        return svc, y_hat, y_hat_train

    def _grid_search(self, kernel_train, x_test, y_train, y_test, c):
        # y_hat is the value predicted

        _, y_hat, _ = self._launch_svc(kernel_train, x_test, y_train, y_test, c)
        res = utils.evaluate_prediction_multiclass(y_test, y_hat)

        return res["balanced_accuracy"]

    def _select_best_parameter(self, async_result):

        c_values = []
        accuracies = []
        for fold in async_result.keys():
            best_c = -1
            best_acc = -1

            for c, async_acc in async_result[fold].items():

                acc = async_acc.get()
                if acc > best_acc:
                    best_c = c
                    best_acc = acc
            c_values.append(best_c)
            accuracies.append(best_acc)

        best_acc = np.mean(accuracies)
        best_c = np.power(10, np.mean(np.log10(c_values)))

        return {"c": best_c, "balanced_accuracy": best_acc}

    def evaluate(self, train_index, test_index):

        inner_pool = ThreadPool(self._algorithm_params["n_threads"])
        async_result = {}
        for i in range(self._algorithm_params["grid_search_folds"]):
            async_result[i] = {}

        outer_kernel = self._kernel[train_index, :][:, train_index]
        y_train = self._y[train_index]

        skf = StratifiedKFold(
            n_splits=self._algorithm_params["grid_search_folds"], shuffle=True
        )
        inner_cv = list(skf.split(np.zeros(len(y_train)), y_train))

        for i in range(len(inner_cv)):
            inner_train_index, inner_test_index = inner_cv[i]

            inner_kernel = outer_kernel[inner_train_index, :][:, inner_train_index]
            x_test_inner = outer_kernel[inner_test_index, :][:, inner_train_index]
            y_train_inner, y_test_inner = (
                y_train[inner_train_index],
                y_train[inner_test_index],
            )

            for c in self._algorithm_params["c_range"]:
                async_result[i][c] = inner_pool.apply_async(
                    self._grid_search,
                    (inner_kernel, x_test_inner, y_train_inner, y_test_inner, c),
                )
        inner_pool.close()
        inner_pool.join()

        best_parameter = self._select_best_parameter(async_result)
        x_test = self._kernel[test_index, :][:, train_index]
        y_train, y_test = self._y[train_index], self._y[test_index]

        _, y_hat, y_hat_train = self._launch_svc(
            outer_kernel, x_test, y_train, y_test, best_parameter["c"]
        )

        result = dict()
        result["best_parameter"] = best_parameter
        result["evaluation"] = utils.evaluate_prediction_multiclass(y_test, y_hat)
        result["evaluation_train"] = utils.evaluate_prediction_multiclass(
            y_train, y_hat_train
        )
        result["y_hat"] = y_hat
        result["y_hat_train"] = y_hat_train
        result["y"] = y_test
        result["y_train"] = y_train
        result["y_index"] = test_index
        result["x_index"] = train_index

        return result

    def apply_best_parameters(self, results_list):

        best_c_list = []
        bal_acc_list = []

        for result in results_list:
            best_c_list.append(result["best_parameter"]["c"])
            bal_acc_list.append(result["best_parameter"]["balanced_accuracy"])

        # 10^(mean of log10 of best Cs of each fold) is selected
        best_c = np.power(10, np.mean(np.log10(best_c_list)))
        # Mean balanced accuracy
        mean_bal_acc = np.mean(bal_acc_list)

        if self._algorithm_params["balanced"]:
            svc = OneVsOneClassifier(
                SVC(
                    C=best_c,
                    kernel="precomputed",
                    probability=True,
                    tol=1e-6,
                    class_weight="balanced",
                )
            )
        else:
            svc = OneVsOneClassifier(
                SVC(C=best_c, kernel="precomputed", probability=True, tol=1e-6)
            )

        svc.fit(self._kernel, self._y)

        return svc, {"c": best_c, "balanced_accuracy": mean_bal_acc}

    def save_classifier(self, classifier, output_dir):

        np.savetxt(
            path.join(output_dir, "support_vectors_indices.txt"), classifier.support_
        )
        np.savetxt(path.join(output_dir, "intersect.txt"), classifier.intercept_)

    def save_weights(self, classifier, x, output_dir):

        dual_coefficients = classifier.dual_coef_
        sv_indices = classifier.support_

        weighted_sv = dual_coefficients.transpose() * x[sv_indices]
        weights = np.sum(weighted_sv, 0)

        np.savetxt(path.join(output_dir, "weights.txt"), weights)

        return weights

    def save_parameters(self, parameters_dict, output_dir):
        with open(path.join(output_dir, "best_parameters.json"), "w") as f:
            json.dump(parameters_dict, f)

    @staticmethod
    def uses_kernel():
        return True

    @staticmethod
    def get_default_parameters():

        parameters_dict = {
            "balanced": True,
            "grid_search_folds": 10,
            "c_range": np.logspace(-6, 2, 17),
            "n_threads": 15,
        }

        return parameters_dict
