
import os
from os import path
import numpy as np

from clinica.pipelines.machine_learning import base, input, algorithm, validation

__author__ = "Jorge Samper Gonzalez"
__copyright__ = "Copyright 2016, The Aramis Lab Team"
__credits__ = ["Jorge Samper Gonzalez"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"

#This code is an example of implementation of a machine learning pipelines

class VB_KFold_DualSVM(base.MLWorkflow):

    #First of all, input has to be chosen. According to it (CAPSVoxelBasedInput or CAPSRegionBasedInput), all the necessary inputs
    #can be found in input.py

    def __init__(self, caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id, image_type, output_dir, fwhm=0,
                 modulated="on", pvc=None, precomputed_kernel=None, mask_zeros=True, n_threads=15, n_folds=10,
                 grid_search_folds=10, balanced=True, c_range=np.logspace(-6, 2, 17)):

        #here some parameters selected for this task

        self._output_dir = output_dir
        self._n_threads = n_threads
        self._n_folds = n_folds
        self._grid_search_folds = grid_search_folds
        self._balanced = balanced
        self._c_range = c_range

        # In this case we are running a voxel based input approach
        #

        self._input = input.CAPSVoxelBasedInput(caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id,
                                                image_type, fwhm, modulated, pvc, mask_zeros, precomputed_kernel)

        #validation and algorithm will be selected in the next part of code

        self._validation = None
        self._algorithm = None

    def run(self):

        # call on parameters already computed

        x = self._input.get_x()
        y = self._input.get_y()
        kernel = self._input.get_kernel()

        #now algorithm has been selected, in this case Dual SVM algorithm. Look at algorithm.py to understand the input necessary for each method
        #input parameters were chosen previously

        self._algorithm = algorithm.DualSVMAlgorithm(kernel,
                                                     y,
                                                     balanced=self._balanced,
                                                     grid_search_folds=self._grid_search_folds,
                                                     c_range=self._c_range,
                                                     n_threads=self._n_threads)
        #here validation type is selected, it's the K fold cross-validation

        self._validation = validation.KFoldCV(self._algorithm)

        classifier, best_params, results = self._validation.validate(y, n_folds=self._n_folds, n_threads=self._n_threads)

        #creation of the path where all the results will be saved

        classifier_dir = path.join(self._output_dir, 'classifier')
        if not path.exists(classifier_dir):
            os.makedirs(classifier_dir)

        #here we have selected whant we wanted save
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
                 modulated="on", pvc=None, precomputed_kernel=None, mask_zeros=True, n_threads=15, n_iterations=100,
                 test_size=0.3, grid_search_folds=10, balanced=True, c_range=np.logspace(-6, 2, 17)):
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

        self._validation = validation.RepeatedHoldOut(self._algorithm, n_iterations=self._n_iterations, test_size=self._test_size)

        classifier, best_params, results = self._validation.validate(y, n_threads=self._n_threads)
        classifier_dir = path.join(self._output_dir, 'classifier')
        if not path.exists(classifier_dir):
            os.makedirs(classifier_dir)

        self._algorithm.save_classifier(classifier, classifier_dir)
        self._algorithm.save_parameters(best_params, classifier_dir)
        weights = self._algorithm.save_weights(classifier, x, classifier_dir)

        self._input.save_weights_as_nifti(weights, classifier_dir)

        self._validation.save_results(self._output_dir)


class RB_RepHoldOut_DualSVM(base.MLWorkflow):

    def __init__(self, caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id, image_type,  atlas,
                 output_dir, pvc=None, n_threads=15, n_iterations=100, test_size=0.3,
                 grid_search_folds=10, balanced=True, c_range=np.logspace(-6, 2, 17)):
        self._output_dir = output_dir
        self._n_threads = n_threads
        self._n_iterations = n_iterations
        self._test_size = test_size
        self._grid_search_folds = grid_search_folds
        self._balanced = balanced
        self._c_range = c_range

        self._input = input.CAPSRegionBasedInput(caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id,
                                                 image_type, atlas, pvc)
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

        self._validation = validation.RepeatedHoldOut(self._algorithm, n_iterations=self._n_iterations, test_size=self._test_size)

        classifier, best_params, results = self._validation.validate(y, n_threads=self._n_threads)
        classifier_dir = path.join(self._output_dir, 'classifier')
        if not path.exists(classifier_dir):
            os.makedirs(classifier_dir)

        self._algorithm.save_classifier(classifier, classifier_dir)
        self._algorithm.save_parameters(best_params, classifier_dir)
        weights = self._algorithm.save_weights(classifier, x, classifier_dir)

        self._input.save_weights_as_nifti(weights, classifier_dir)

        self._validation.save_results(self._output_dir)


class RB_RepHoldOut_LogisticRegression(base.MLWorkflow):
    
    def __init__(self, caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id, image_type, atlas,
                 output_dir, pvc=None, n_threads=15, n_iterations=100, test_size=0.3,
                 grid_search_folds=10, balanced=True, c_range=np.logspace(-6, 2, 17)):
        self._output_dir = output_dir
        self._n_threads = n_threads
        self._n_iterations = n_iterations
        self._test_size = test_size
        self._grid_search_folds = grid_search_folds
        self._balanced = balanced
        self._c_range = c_range

        self._input = input.CAPSRegionBasedInput(caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id,
                                                 image_type, atlas, pvc)
        self._validation = None
        self._algorithm = None
    
    def run(self):
        
        x = self._input.get_x()
        y = self._input.get_y()
        
        self._algorithm = algorithm.LogisticReg(x, y, balanced=self._balanced,
                                                grid_search_folds=self._grid_search_folds,
                                                c_range=self._c_range,
                                                n_threads=self._n_threads)
            
        self._validation = validation.RepeatedHoldOut(self._algorithm, n_iterations=self._n_iterations, test_size=self._test_size)
        classifier, best_params, results = self._validation.validate(y, n_threads=self._n_threads)
                                                     
        classifier_dir = os.path.join(self._output_dir, 'classifier')
        if not path.exists(classifier_dir):
            os.makedirs(classifier_dir)

        self._algorithm.save_classifier(classifier, classifier_dir)
        self._algorithm.save_parameters(best_params, classifier_dir)
        weights = self._algorithm.save_weights(classifier, classifier_dir)

        self._classifier = classifier

        self._input.save_weights_as_nifti(weights, classifier_dir)

        self._validation.save_results(self._output_dir)


class RB_RepHoldOut_RandomForest(base.MLWorkflow):

    def __init__(self, caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id, image_type, atlas,
                 output_dir, pvc=None, n_threads=15, n_iterations=100, test_size=0.3,
                 grid_search_folds=10, balanced=True, n_estimators_range=range(5, 15, 1)):
        self._output_dir = output_dir
        self._n_threads = n_threads
        self._n_iterations = n_iterations
        self._test_size = test_size
        self._grid_search_folds = grid_search_folds
        self._balanced = balanced
        self._n_estimators_range = n_estimators_range
        self._input = input.CAPSRegionBasedInput(caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id,
                                                 image_type, atlas, pvc)
        self._validation = None
        self._algorithm = None

    def run(self):

        x = self._input.get_x()
        y = self._input.get_y()

        self._algorithm = algorithm.RandomForest(x, y, balanced=self._balanced,
                                                 grid_search_folds=self._grid_search_folds,
                                                 n_estimators_range=self._n_estimators_range,
                                                 n_threads=self._n_threads)

        self._validation = validation.RepeatedHoldOut(self._algorithm, n_iterations=self._n_iterations, test_size=self._test_size)
        classifier, best_params, results = self._validation.validate(y, n_threads=self._n_threads)

        classifier_dir = os.path.join(self._output_dir, 'classifier')
        if not path.exists(classifier_dir):
            os.makedirs(classifier_dir)

        self._algorithm.save_classifier(classifier, classifier_dir)
        self._algorithm.save_parameters(best_params, classifier_dir)
        weights = self._algorithm.save_weights(classifier, classifier_dir)

        self._input.save_weights_as_nifti(weights, classifier_dir)

        self._validation.save_results(self._output_dir)


class RB_LearningCurveRepHoldOut_DualSVM(base.MLWorkflow):

    def __init__(self, caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id, image_type,  atlas,
                 output_dir, pvc=None, precomputed_kernel=None, n_threads=15, n_iterations=100, test_size=0.3,
                 n_learning_points=10, grid_search_folds=10, balanced=True, c_range=np.logspace(-6, 2, 17)):

        self._output_dir = output_dir
        self._n_threads = n_threads
        self._n_iterations = n_iterations
        self._test_size = test_size
        self._n_learning_points = n_learning_points
        self._grid_search_folds = grid_search_folds
        self._balanced = balanced
        self._c_range = c_range

        self._input = input.CAPSRegionBasedInput(caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id,
                                                 image_type, atlas, pvc, precomputed_kernel)
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

        self._validation = validation.LearningCurveRepeatedHoldOut(self._algorithm,
                                                                   n_iterations=self._n_iterations,
                                                                   test_size=self._test_size,
                                                                   n_learning_points=self._n_learning_points)

        classifier, best_params, results = self._validation.validate(y, n_threads=self._n_threads)

        for learning_point in range(self._n_learning_points):

            learning_point_dir = path.join(self._output_dir, 'learning_split-' + str(learning_point))

            classifier_dir = path.join(learning_point_dir, 'classifier')
            if not path.exists(classifier_dir):
                os.makedirs(classifier_dir)

            self._algorithm.save_classifier(classifier[learning_point], classifier_dir)
            self._algorithm.save_parameters(best_params[learning_point], classifier_dir)
            weights = self._algorithm.save_weights(classifier[learning_point], x, classifier_dir)

            self._input.save_weights_as_nifti(weights, classifier_dir)

        self._validation.save_results(self._output_dir)


class VB_LearningCurveRepHoldOut_DualSVM(base.MLWorkflow):

    def __init__(self, caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id, image_type, output_dir, fwhm=0,
                 modulated="on", pvc=None, precomputed_kernel=None, mask_zeros=True, n_threads=15, n_iterations=100,
                 test_size=0.3, n_learning_points=10, grid_search_folds=10, balanced=True, c_range=np.logspace(-6, 2, 17)):
        self._output_dir = output_dir
        self._n_threads = n_threads
        self._n_iterations = n_iterations
        self._test_size = test_size
        self._n_learning_points = n_learning_points
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

        self._validation = validation.LearningCurveRepeatedHoldOut(self._algorithm,
                                                                   n_iterations=self._n_iterations,
                                                                   test_size=self._test_size,
                                                                   n_learning_points=self._n_learning_points)

        classifier, best_params, results = self._validation.validate(y, n_threads=self._n_threads)

        for learning_point in range(self._n_learning_points):

            learning_point_dir = path.join(self._output_dir, 'learning_split-' + str(learning_point))

            classifier_dir = path.join(learning_point_dir, 'classifier')
            if not path.exists(classifier_dir):
                os.makedirs(classifier_dir)

            self._algorithm.save_classifier(classifier[learning_point], classifier_dir)
            self._algorithm.save_parameters(best_params[learning_point], classifier_dir)
            weights = self._algorithm.save_weights(classifier[learning_point], x, classifier_dir)

            self._input.save_weights_as_nifti(weights, classifier_dir)

        self._validation.save_results(self._output_dir)
