"""Statistics Surfstat3 - Clinica Command Line Interface.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details: https://gitlab.icm-institute.org/aramis/clinica/wikis/docs/InteractingWithClinica.
"""


import clinica.engine as ce

__author__ = "Junhao Wen"
__copyright__ = "Copyright 2016, The Aramis Lab Team"
__credits__ = ["Michael Bacci", "Junhao Wen"]
__license__ = "See LICENSE.txt file"
__version__ = "1.0.0"
__maintainer__ = "Junhao Wen"
__email__ = "junhao.Wen@inria.fr"
__status__ = "Development"

class StatisticsSurfstatCLI(ce.CmdParser):


    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """

        self._name = 'statistics-surfstat'


    def define_options(self):
        """Define the sub-command arguments
        """
        self._args.add_argument("caps_directory",
                                help='Path to the CAPS directory.')
        self._args.add_argument("tsv_file",
                                help='Path to the tsv containing the information for GLM.')
        self._args.add_argument("design_matrix",
                                help='A str to define the design matrix that fits into GLM, eg, 1 + group + sex + age')
        self._args.add_argument("contrast",
                                help='A str to define the contrast matrix for GLM, eg, group_label')
        self._args.add_argument("str_format",
                                help='A str to define the format string for the tsv column , eg, %%s %%s %%s %%f')
        self._args.add_argument("group_label",
                                help='A str for current group name')
        self._args.add_argument("glm_type",
                                help='A str based on glm type for the hypothesis, choose one between group_comparison and correlation')
        self._args.add_argument("-fwhm", "--full_width_at_half_maximum", type=int, default=20,
                                help='FWHM for the surface smoothing (default=20)')
        self._args.add_argument("-tup", "--threshold_uncorrected_pvalue", type=float, default=0.001,
                                help='Threshold to display the uncorrected Pvalue (default=0.001)')
        self._args.add_argument("-tcp", "--threshold_corrected_pvalue", type=float, default=0.05,
                                help='Threshold to display the corrected cluster (default=0.05)')
        self._args.add_argument("-ct", "--cluster_threshold", type=float, default=0.001,
                                help='Threshold to define a cluster in the process of cluster-wise correction (default=0.001)')
        self._args.add_argument("-np", "--n_procs", type=int, default=4,
                                help='Number of parallel processes to run (default=4)')
        self._args.add_argument("-wd", "--working_directory", type=str, default=None,
                                help='Temporary directory to run the workflow')


    def run_pipeline(self, args):
        """
        Run the pipeline with defined args
        """

        from statistics_surfstat_pipeline import StatisticsSurfstat

        pipeline = StatisticsSurfstat(
            # pass these args by the class attribute itself
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.tsv_file))
        pipeline.parameters = {
            # pass these args by using self.parameters in a dictionary
            'design_matrix': args.design_matrix,
            'contrast': args.contrast,
            'str_format': args.str_format,
            'group_label': args.group_label,
            'glm_type': args.glm_type,
            'full_width_at_half_maximum': args.full_width_at_half_maximum or 20,
            'threshold_uncorrected_pvalue': args.threshold_uncorrected_pvalue or 0.001,
            'threshold_corrected_pvalue': args.threshold_corrected_pvalue or 0.05,
            'cluster_threshold': args.cluster_threshold or 0.001
        }
        pipeline.base_dir = self.absolute_path(args.working_directory)

        # run the pipeline in n_procs cores based on your computation power.
        if args.n_procs:
            pipeline.write_graph()
            pipeline.run(plugin='MultiProc', plugin_args={'n_procs': args.n_procs})
        else:
            pipeline.write_graph()
            pipeline.run()