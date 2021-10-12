from abc import ABC, abstractmethod


class MLWorkflow(ABC):
    def __init__(
        self, input_class, validation_class, algorithm_class, all_params, output_dir
    ):

        self._input_class = input_class
        self._validation_class = validation_class
        self._algorithm_class = algorithm_class

        self._input_params = self.create_parameters_dict(all_params, input_class)
        self._validation_params = self.create_parameters_dict(
            all_params, validation_class
        )
        self._algorithm_params = self.create_parameters_dict(
            all_params, algorithm_class
        )

        self._output_dir = output_dir

        self._input = None
        self._validation = None
        self._algorithm = None

    def run(self):

        from os import makedirs, path

        # Instantiating input class
        self._input = self._input_class(self._input_params)

        # Computing input values
        x = self._input.get_x()
        y = self._input.get_y()

        # Instantiating classification algorithm
        if self._algorithm_class.uses_kernel():
            kernel = self._input.get_kernel()
            self._algorithm = self._algorithm_class(kernel, y, self._algorithm_params)
        else:
            self._algorithm = self._algorithm_class(x, y, self._algorithm_params)

        # Instantiating cross-validation method and classification algorithm
        self._validation = self._validation_class(
            self._algorithm, self._validation_params
        )

        # Launching classification with selected cross-validation
        classifier, best_params, results = self._validation.validate(y)

        # Creation of the directory to save results
        classifier_dir = path.join(self._output_dir, "classifier")
        makedirs(classifier_dir, exist_ok=True)

        # Saving algorithm trained classifier
        self._algorithm.save_classifier(classifier, classifier_dir)
        self._algorithm.save_weights(classifier, x, classifier_dir)
        self._algorithm.save_parameters(best_params, classifier_dir)

        # Saving validation trained classifier
        self._validation.save_results(self._output_dir)

    @staticmethod
    def create_parameters_dict(locals_dictionary, component_class):

        default_parameters = component_class.get_default_parameters()
        for key in locals_dictionary:
            if key in default_parameters:
                default_parameters[key] = locals_dictionary[key]
        return default_parameters


class MLInput(ABC):
    def __init__(self, input_params):

        self._input_params = self.get_default_parameters()
        self._input_params.update(input_params)

        self._x = None
        self._y = None
        self._kernel = None

    @abstractmethod
    def get_x(self):
        pass

    @abstractmethod
    def get_y(self):
        pass

    @staticmethod
    @abstractmethod
    def get_default_parameters():
        pass


class MLValidation(ABC):
    def __init__(self, ml_algorithm, validation_params):

        self._ml_algorithm = ml_algorithm

        self._validation_params = self.get_default_parameters()
        self._validation_params.update(validation_params)

        self._validation_results = []
        self._classifier = None
        self._best_params = None

    @abstractmethod
    def validate(self, y):
        pass

    @staticmethod
    @abstractmethod
    def get_default_parameters():
        pass


class MLAlgorithm(ABC):
    def __init__(self, input_data, y, algorithm_params):

        self._algorithm_params = self.get_default_parameters()
        self._algorithm_params.update(algorithm_params)

        if self.uses_kernel():
            self._kernel = input_data
        else:
            self._x = input_data

        self._y = y

    @staticmethod
    @abstractmethod
    def uses_kernel():
        pass

    @abstractmethod
    def evaluate(self, train_index, test_index):
        pass

    @abstractmethod
    def save_classifier(self, classifier, output_dir):
        pass

    @abstractmethod
    def save_parameters(self, parameters, output_dir):
        pass

    @staticmethod
    @abstractmethod
    def get_default_parameters():
        pass
