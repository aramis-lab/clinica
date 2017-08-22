
import abc


class MLWorkflow:
    __metaclass__ = abc.ABCMeta

    def __init__(self, ml_input, ml_validation, ml_algorithm):
        self._ml_input = ml_input
        self._ml_validation = ml_validation
        self._ml_algorithm = ml_algorithm


class MLInput:
    __metaclass__ = abc.ABCMeta
    pass


class MLValidation:
    __metaclass__ = abc.ABCMeta
    pass


class MLAlgorithm:
    __metaclass__ = abc.ABCMeta
    pass
