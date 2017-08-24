
import abc


class MLWorkflow:
    __metaclass__ = abc.ABCMeta

    # def __init__(self, ml_input, ml_validation, ml_algorithm):
    #     self._ml_input = ml_input
    #     self._ml_validation = ml_validation
    #     self._ml_algorithm = ml_algorithm

    @abc.abstractmethod
    def run(self):
        pass


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
