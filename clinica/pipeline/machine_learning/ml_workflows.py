
import numpy as np
from os import path
import json
from clinica.pipeline.machine_learning import base, input, algorithm, validation
import clinica.pipeline.machine_learning.svm_utils as utils


class VoxelBasedCVDualSVMWorkflow(base.MLWorkflow):

    def __init__(self, caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id, output_dir, fwhm=0, modulated="on",
                 precomputed_kernel=None, mask_zeros=True, n_threads=15, n_folds=10, grid_search_folds=10, balanced=True,
                 c_range=np.logspace(-6, 2, 17)):
        self._output_dir = output_dir
        self._n_threads = n_threads
        self._n_folds = n_folds
        self._grid_search_folds = grid_search_folds
        self._balanced = balanced
        self._c_range = c_range

        self._input = input.CAPST1Input(caps_directory, subjects_visits_tsv, diagnoses_tsv, group_id, fwhm, modulated, mask_zeros, precomputed_kernel)
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

        self._algorithm.save_classifier(classifier, self._output_dir)
        self._algorithm.save_weights(classifier, x, self._output_dir)
        self._algorithm.save_parameters(best_params, self._output_dir)

        self._validation.save_results(self._output_dir)
