
import abc

__author__ = "Jorge Samper Gonzalez"
__copyright__ = "Copyright 2016, The Aramis Lab Team"
__credits__ = ["Jorge Samper Gonzalez"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


class MLWorkflow:
    __metaclass__ = abc.ABCMeta

    # def __init__(self, ml_input, ml_validation, ml_algorithm, output_dir):
    #     self._ml_input = ml_input
    #     self._ml_validation = ml_validation
    #     self._ml_algorithm = ml_algorithm
    #     self._output_dir = output_dir

    @abc.abstractmethod
    def run(self):
        pass

    def save_image(self):
        import os
        import pandas as pd
        from neuropredict import visualize

        df = pd.io.parsers.read_csv(os.path.join(self._output_dir, 'results.tsv'), sep='\t')
        visualize.metric_distribution(df.as_matrix(['balanced_accuracy', 'auc', 'accuracy']),
                                      ['balanced_accuracy', 'auc', 'accuracy'],
                                      os.path.join(self._output_dir, 'metrics'))


class MLInput:
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def get_x(self):
        pass

    @abc.abstractmethod
    def get_y(self):
        pass


class MLValidation:
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def validate(self, y):
        pass


class MLAlgorithm:
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def evaluate(self, train_index, test_index):
        pass

    @abc.abstractmethod
    def save_classifier(self, classifier, output_dir):
        pass

    @abc.abstractmethod
    def save_parameters(self, parameters, output_dir):
        pass
