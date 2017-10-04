"""Statistics Surfstat3 - Clinica Command Line Interface.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details: https://gitlab.icm-institute.org/aramis/clinica/wikis/docs/InteractingWithClinica.
"""


import clinica.engine as ce

__author__ = "Junhao Wen"
__copyright__ = "Copyright 2016, The Aramis Lab Team"
__credits__ = ["Michael Bacci", "Junhao Wen"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
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
                                help='A str to define the contrast matrix for GLM, eg, group. Note, when you want to negative correlation, there is a bug for Clinica commandline currently')
        self._args.add_argument("str_format",
                                help='A str to define the format string for the tsv column , eg, %%s %%s %%s %%f')
        self._args.add_argument("group_label",
                                help='A str for current group name')
        self._args.add_argument("glm_type",
                                help='A str based on glm type for the hypothesis, choose one between group_comparison and correlation')
        self._args.add_argument("-ft", "--feature_type", type=str, default=None,
                                help='Feature type. Can be : cortical_thickness, pet_fdg_projection, pet_noddi_projection_ndi and pet_noddi_projection_odi. Default = cortical_thickness')
        self._args.add_argument("-cf", "--custom_file", type=str, default=None,
                                help='Pattern of file inside caps directory using @subject, @session, @fwhm, @hemi. No --feature_type must be specified in order to use this flag.')
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
        from statistics_surfstat_utils import check_inputs
        from clinica.utils.stream import cprint
        import os

        if args.feature_type is not None:
            if args.custom_file is not None:
                raise Exception('--feature_type and --custom_file are mutually exclusive : you must choose between one or the other. See documentation for more informations.')
            if args.feature_type == 'cortical_thickness':
                args.custom_file = '@subject/@session/t1/freesurfer_cross_sectional/@subject_@session/surf/@hemi.thickness.fwhm@fwhm.fsaverage.mgh'
            elif args.feature_type == 'pet_fdg_projection':
                args.custom_file = '@subject/@session/pet/surface/@subject_@session_task-rest_acq-FDG_pet_space-fsaverage_suvr-pons_pvc-iy_hemi-@hemi_fwhm-@fwhm_projection.mgh'
            elif args.feature_type == 'pet_noddi_projection_ndi':
                args.custom_file = '@subject/@session/noddi/postprocessing/noddi-register-vertex-fsaverage/cortex-projection/@subject_@session_OnFsaverage_fwhm-@fwhm_measure-ficvf_hemi-@hemi.mgh'
            elif args.feature_type == 'pet_noddi_projection_odi':
                args.custom_file = '@subject/@session/noddi/postprocessing/noddi-register-vertex-fsaverage/cortex-projection/@subject_@session_OnFsaverage_fwhm-@fwhm_measure-odi_hemi-@hemi.mgh'
            else:
                raise Exception('Feature type ' + args.feature_type + ' not recognized. Use --custom_file to specify your own files (without --feature_type).')
        elif args.feature_type is None:
            if args.custom_file is None:
                cprint('No feature type selected : using cortical thickness as default value')
                args.custom_file = '@subject/@session/t1/freesurfer_cross_sectional/@subject_@session/surf/@hemi.thickness.fwhm@fwhm.fsaverage.mgh'
            else:
                cprint('Using custom features.')

        # Check that group label is alphanumerical
        if not args.group_label.isalnum():
            raise ValueError('Not valid group_id value. It must be composed only by letters and/or numbers')

        # Check that group label does not already exists in CAPS folder
        #if os.path.exists(os.path.join(os.path.abspath(self.absolute_path(args.caps_directory)), 'groups', 'group-' + args.group_label)):
        #    error_message = 'group_id : ' + args.group_label + ' already exists, please choose an other one. Groups that exists in your CAPS directory are : \n'
        #    list_groups = os.listdir(os.path.join(os.path.abspath(self.absolute_path(args.caps_directory)), 'groups'))
        #    for e in list_groups:
        #       if e.startswith('group-'):
        #            error_message += e + ' \n'
        #   raise ValueError(error_message)

        pipeline = StatisticsSurfstat(
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.tsv_file))
        pipeline.parameters = {
            # pass these args by using self.parameters in a dictionary
            'design_matrix': args.design_matrix,
            'contrast': args.contrast,
            'str_format': args.str_format,
            'group_label': args.group_label,
            'glm_type': args.glm_type,
            'custom_file': args.custom_file,
            'full_width_at_half_maximum': args.full_width_at_half_maximum or 20,
            'threshold_uncorrected_pvalue': args.threshold_uncorrected_pvalue or 0.001,
            'threshold_corrected_pvalue': args.threshold_corrected_pvalue or 0.05,
            'cluster_threshold': args.cluster_threshold or 0.001
        }
        pipeline.base_dir = self.absolute_path(args.working_directory)

        check_inputs(pipeline.caps_directory,
                     pipeline.parameters['custom_file'],
                     pipeline.parameters['full_width_at_half_maximum'],
                     pipeline.tsv_file)

        # run the pipeline in n_procs cores based on your computation power.
        if args.n_procs:
            pipeline.write_graph()
            pipeline.run(plugin='MultiProc', plugin_args={'n_procs': args.n_procs})
        else:
            pipeline.write_graph()
            print(pipeline.parameters)
            pipeline.run()
