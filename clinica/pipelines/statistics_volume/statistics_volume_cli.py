# coding: utf-8

import clinica.engine as ce


class StatisticsVolumeCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline."""
        self._name = 'statistics-volume'

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = ('Volume-based mass-univariate analysis with SPM:\n'
                             'http://clinica.run/doc/Pipelines/Statistics_Volume/')

    def define_options(self):
        """Define the sub-command arguments."""
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES

        # Clinica compulsory arguments
        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')
        clinica_comp.add_argument("subject_visits_with_covariates_tsv",
                                  help='TSV file containing a list of subjects with their sessions and all '
                                       'the covariates and factors needed for the GLM.')
        clinica_comp.add_argument("contrast",
                                  help='Defines the contrast. Must be one of the column names form the TSV file.')
        clinica_comp.add_argument("feature_type",
                                  help='Type of volume-based feature: graymatter (from t1-volume pipeline) or '
                                       'fdg (from pet-volume pipeline). Use your own if you want to use '
                                       'the --custom_file flag. ')
        clinica_comp.add_argument("group_id",
                                  help='User-defined identifier for the provided group of subjects.')

        # Optional arguments (e.g. FWHM)
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES['OPTIONAL'])
        optional.add_argument("--custom_files", "-cf", type=str, default=None,
                              help='Custom file string. Specify filename using * when the subject or session name '
                                   'appears e.g. \'*_task-rest_acq-fdg_pet_space-Ixi549Space_pet.nii.gz\' will grab '
                                   'the corresponding file in all the subjects/sessions.')
        optional.add_argument("-gic", "--group_id_caps", type=str, default=None,
                              help='Name of the group that Clinica needs to use to grab input file.')
        optional.add_argument("-fwhm", "--full_width_at_half_maximum", type=int, default=8,
                              help='Full Width at Half Maximum (FWHM) of the smoothing used in your input file '
                                   '(default: --full_width_at_half_maximum %(default)s).')

        # Clinica standard arguments (e.g. --n_procs)
        self.add_clinica_standard_arguments(add_tsv_flag=False)

        # Advanced arguments (i.e. tricky parameters)
        advanced = self._args.add_argument_group(PIPELINE_CATEGORIES['ADVANCED'])
        advanced.add_argument("-ct", "--cluster_threshold",
                              type=float, default=0.001,
                              help='Threshold to define a cluster in the process of cluster-wise correction '
                                   '(default: --cluster_threshold %(default)s).')

    def run_command(self, args):
        from networkx import Graph
        from .statistics_volume_pipeline import StatisticsVolume
        from clinica.utils.ux import print_end_pipeline, print_crash_files_and_exit

        parameters = {
            'contrast': args.contrast,
            'feature_type': args.feature_type,
            'group_id': args.group_id,
            'custom_files': args.custom_files,
            'cluster_threshold': args.cluster_threshold,
            'group_id_caps': args.group_id_caps,
            'full_width_at_half_maximum': args.full_width_at_half_maximum
        }

        pipeline = StatisticsVolume(
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subject_visits_with_covariates_tsv),
            base_dir=self.absolute_path(args.working_directory),
            parameters=parameters,
            name=self.name
        )

        if args.n_procs:
            exec_pipeline = pipeline.run(plugin='MultiProc',
                                         plugin_args={'n_procs': args.n_procs})
        else:
            exec_pipeline = pipeline.run()

        if isinstance(exec_pipeline, Graph):
            print_end_pipeline(self.name, pipeline.base_dir, pipeline.base_dir_was_specified)
        else:
            print_crash_files_and_exit(args.logname, pipeline.base_dir)
