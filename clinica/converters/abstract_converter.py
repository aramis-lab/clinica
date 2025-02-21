import abc


class Converter:
    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def convert_images(self, src, dst):
        pass

    @abc.abstractmethod
    def convert_clinical_data(self, src, dst):
        pass
