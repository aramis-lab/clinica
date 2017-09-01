
import os
from os import path
import numpy as np

from clinica.pipeline.machine_learning import base, input, algorithm, validation


class VB_KFold_DualSVM(base.MLWorkflow):

    def __init__(self, caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id, image_type, output_dir, fwhm=0,
                 modulated="on", pvc=None, precomputed_kernel=None, mask_zeros=True, n_threads=15, n_folds=10,
                 grid_search_folds=10, balanced=True, c_range=np.logspace(-6, 2, 17)):
        self._output_dir = output_dir
        self._n_threads = n_threads
        self._n_folds = n_folds
        self._grid_search_folds = grid_search_folds
        self._balanced = balanced
        self._c_range = c_range

        self._input = input.CAPSVoxelBasedInput(caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id,
                                                image_type, fwhm, modulated, pvc, mask_zeros, precomputed_kernel)
        self._validation = None
        self._algorithm = None

    def run(self):

        x = self._input.get_x()
        y = self._input.get_y()
        kernel = self._input.get_kernel()

        self._algorithm = algorithm.DualSVMAlgorithm(kernel,
                                                     y,
                                                     balanced=self._balanced,
                                                     grid_search_folds=self._grid_search_folds,
                                                     c_range=self._c_range,
                                                     n_threads=self._n_threads)

        self._validation = validation.KFoldCV(self._algorithm)

        classifier, best_params, results = self._validation.validate(y, n_folds=self._n_folds, n_threads=self._n_threads)

        classifier_dir = path.join(self._output_dir, 'classifier')
        if not path.exists(classifier_dir):
            os.makedirs(classifier_dir)

        self._algorithm.save_classifier(classifier, classifier_dir)
        self._algorithm.save_weights(classifier, x, classifier_dir)
        self._algorithm.save_parameters(best_params, classifier_dir)

        self._validation.save_results(self._output_dir)

        # self._input.save_weights_as_nifti(weights)


class VB_RepKFold_DualSVM(base.MLWorkflow):

    def __init__(self, caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id, image_type, output_dir, fwhm=0,
                 modulated="on", precomputed_kernel=None, mask_zeros=True, n_threads=15, n_iterations=100, n_folds=10,
                 grid_search_folds=10, balanced=True, c_range=np.logspace(-6, 2, 17)):
        self._output_dir = output_dir
        self._n_threads = n_threads
        self._n_iterations = n_iterations
        self._n_folds = n_folds
        self._grid_search_folds = grid_search_folds
        self._balanced = balanced
        self._c_range = c_range

        self._input = input.CAPSVoxelBasedInput(caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id,
                                                image_type, fwhm, modulated, mask_zeros, precomputed_kernel)
        self._validation = None
        self._algorithm = None

    def run(self):

        x = self._input.get_x()
        y = self._input.get_y()
        kernel = self._input.get_kernel()

        self._algorithm = algorithm.DualSVMAlgorithm(kernel,
                                                     y,
                                                     balanced=self._balanced,
                                                     grid_search_folds=self._grid_search_folds,
                                                     c_range=self._c_range,
                                                     n_threads=self._n_threads)

        self._validation = validation.RepeatedKFoldCV(self._algorithm)

        classifier, best_params, results = self._validation.validate(y, n_iterations=self._n_iterations,
                                                                     n_folds=self._n_folds, n_threads=self._n_threads)

        classifier_dir = path.join(self._output_dir, 'classifier')
        if not path.exists(classifier_dir):
            os.makedirs(classifier_dir)

        self._algorithm.save_classifier(classifier, classifier_dir)
        weights = self._algorithm.save_weights(classifier, x, classifier_dir)
        self._algorithm.save_parameters(best_params, classifier_dir)

        self._validation.save_results(self._output_dir)

        self._input.save_weights_as_nifti(weights, classifier_dir)


class VB_RepHoldOut_DualSVM(base.MLWorkflow):

    def __init__(self, caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id, image_type, output_dir, fwhm=0,
                 modulated="on", pvc=None, precomputed_kernel=None, mask_zeros=True, n_threads=15, n_iterations=100, test_size=0.3,
                 grid_search_folds=10, balanced=True, c_range=np.logspace(-6, 2, 17)):
        self._output_dir = output_dir
        self._n_threads = n_threads
        self._n_iterations = n_iterations
        self._test_size = test_size
        self._grid_search_folds = grid_search_folds
        self._balanced = balanced
        self._c_range = c_range

        self._input = input.CAPSVoxelBasedInput(caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id,
                                                image_type, fwhm, modulated, pvc, mask_zeros, precomputed_kernel)

        self._validation = None
        self._algorithm = None

    def run(self):

        x = self._input.get_x()
        y = self._input.get_y()
        kernel = self._input.get_kernel()

        self._algorithm = algorithm.DualSVMAlgorithm(kernel,
                                                     y,
                                                     balanced=self._balanced,
                                                     grid_search_folds=self._grid_search_folds,
                                                     c_range=self._c_range,
                                                     n_threads=self._n_threads)

        self._validation = validation.RepeatedSplit(self._algorithm, n_iterations=self._n_iterations, test_size=self._test_size)

        classifier, best_params, results = self._validation.validate(y, n_threads=self._n_threads)
        classifier_dir = path.join(self._output_dir, 'classifier')
        if not path.exists(classifier_dir):
            os.makedirs(classifier_dir)

        self._algorithm.save_classifier(classifier, classifier_dir)
        weights = self._algorithm.save_weights(classifier, x, classifier_dir)
        self._algorithm.save_parameters(best_params, classifier_dir)

        self._validation.save_results(self._output_dir)

        self._input.save_weights_as_nifti(weights, classifier_dir)


class VB_RepHoldOut_LogisticRegression(base.MLWorkflow):
    
    def __init__(self, caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id, image_type,
                 output_dir, fwhm=0, modulated="on", mask_zeros=True, n_threads=15,
                 n_iterations=100, test_size=0.3,
                 grid_search_folds=10, balanced=True, voxelbased=False, c_range=np.logspace(-6, 2, 17)):
        self._output_dir = output_dir
        self._n_threads = n_threads
        self._n_iterations = n_iterations
        self._test_size = test_size
        self._grid_search_folds = grid_search_folds
        self._balanced = balanced
        self._c_range = c_range
        
        if voxelbased:
            self._input = input.CAPSVoxelBasedInput(caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id, image_type, fwhm, modulated, mask_zeros, None)
        else:
            self._input = input.CAPSRegionBasedInput(caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id, image_type, 'AAL2', modulated, mask_zeros, None)
        self._validation = None
        self._algorithm = None
    
    def run(self):
        
        x = self._input.get_x()
        y = self._input.get_y()
        
        self._algorithm = algorithm.LogisticReg(x, y, balanced=self._balanced,
                                                grid_search_folds=self._grid_search_folds,
                                                c_range=self._c_range,
                                                n_threads=self._n_threads)
            
        self._validation = validation.RepeatedSplit(self._algorithm, n_iterations=self._n_iterations, test_size=self._test_size)
        classifier, best_params, results = self._validation.validate(y, n_threads=self._n_threads)
                                                     
        classifier_dir = os.path.join(self._output_dir, 'classifier')
        if not path.exists(classifier_dir):
            os.makedirs(classifier_dir)

        self._algorithm.save_classifier(classifier, classifier_dir)
        self._algorithm.save_parameters(best_params, classifier_dir)
        self._validation.save_results(self._output_dir)


class RB_RepHoldOut_RandomForest(base.MLWorkflow):

    def __init__(self, caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id, image_type,
                 output_dir, fwhm=0, modulated="on", mask_zeros=True, n_threads=15,
                 n_splits=100, test_size=0.3,
                 grid_search_folds=10, balanced=True, voxelbased=False, n_estimators_range=range(5, 15, 1)):
        self._output_dir = output_dir
        self._n_threads = n_threads
        self._n_splits = n_splits
        self._test_size = test_size
        self._grid_search_folds = grid_search_folds
        self._balanced = balanced
        self._n_estimators_range = n_estimators_range
        self._input = input.CAPSRegionBasedInput(caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id, image_type, 'AAL2', modulated, mask_zeros, None)
        self._validation = None
        self._algorithm = None

    def run(self):

        x = self._input.get_x()
        y = self._input.get_y()

        self._algorithm = algorithm.RandomForest(x, y, balanced=self._balanced,
                                                 grid_search_folds=self._grid_search_folds,
                                                 n_estimators_range=self._n_estimators_range,
                                                 n_threads=self._n_threads)

        self._validation = validation.RepeatedSplit(self._algorithm, n_iterations=self._n_splits, test_size=self._test_size)
        classifier, best_params, results = self._validation.validate(y, n_threads=self._n_threads)

        classifier_dir = os.path.join(self._output_dir, 'classifier')
        if not path.exists(classifier_dir):
            os.makedirs(classifier_dir)

        self._algorithm.save_classifier(classifier, classifier_dir)
        self._algorithm.save_parameters(best_params, classifier_dir)
        self._validation.save_results(self._output_dir)





