# coding: utf8


import abc

__author__ = "Jorge Samper-Gonzalez"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Jorge Samper-Gonzalez"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper-Gonzalez"
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

        pd.io.parsers.read_csv(os.path.join(self._output_dir, 'results.tsv'), sep='\t')

    @staticmethod
    def metric_distribution(metric, labels, output_path, num_classes=2, metric_label='balanced accuracy'):
        """

        Distribution plots of various metrics such as balanced accuracy!

        metric is expected to be ndarray of size [num_repetitions, num_datasets]

        """
        # from __future__ import print_function, division

        import numpy as np
        import matplotlib.pyplot as plt
        from matplotlib import cm
        from matplotlib.backends.backend_pdf import PdfPages

        num_repetitions = metric.shape[0]
        num_datasets = metric.shape[1]
        assert len(labels) == num_datasets, "Differing number of features and labels!"
        method_ticks = 1.0 + np.arange(num_datasets)

        fig, ax = plt.subplots(figsize=[9, 9])
        line_coll = ax.violinplot(metric, widths=0.8, bw_method=0.2,
                                  showmedians=True, showextrema=False,
                                  positions=method_ticks)

        cmap = cm.get_cmap('Paired', num_datasets)
        for cc, ln in enumerate(line_coll['bodies']):
            ln.set_facecolor(cmap(cc))
            ln.set_label(labels[cc])

        plt.legend(loc=2, ncol=num_datasets)

        ax.tick_params(axis='both', which='major', labelsize=15)
        ax.grid(axis='y', which='major')

        lower_lim = np.round(np.min([np.float64(0.9 / num_classes), metric.min()]), 3)
        upper_lim = np.round(np.max([1.01, metric.max()]), 3)
        step_tick = 0.1
        ax.set_ylim(lower_lim, upper_lim)

        ax.set_xticks(method_ticks)
        ax.set_xlim(np.min(method_ticks) - 1, np.max(method_ticks) + 1)
        ax.set_xticklabels(labels, rotation=45)  # 'vertical'

        ax.set_yticks(np.arange(lower_lim, upper_lim, step_tick))
        ax.set_yticklabels(np.arange(lower_lim, upper_lim, step_tick))
        # plt.xlabel(xlabel, fontsize=16)
        plt.ylabel(metric_label, fontsize=16)

        fig.tight_layout()

        pp1 = PdfPages(output_path + '.pdf')
        pp1.savefig()
        pp1.close()

        return


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
