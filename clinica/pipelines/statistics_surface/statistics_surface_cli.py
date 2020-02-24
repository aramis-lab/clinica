# coding: utf8

import clinica.engine as ce


class StatisticsSurfaceCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline."""
        self._name = 'statistics-surface'

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = ('Surface-based mass-univariate analysis with SurfStat:\n'
                             'http://clinica.run/doc/Pipelines/Stats_Surface/')

    def define_options(self):
        """Define the sub-command arguments."""
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES
        from colorama import Fore
        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_id)
        clinica_comp = self._args.add_argument_group('%sMandatory arguments%s' % (Fore.BLUE, Fore.RESET))
        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')
        clinica_comp.add_argument("subject_visits_with_covariates_tsv",
                                  help='TSV file containing a list of subjects with their sessions and all '
                                       'the covariates and factors needed for the GLM.')
        clinica_comp.add_argument("design_matrix",
                                  help='String to define the design matrix that fits into the GLM, '
                                       'e.g. 1 + group + sex + age.')
        clinica_comp.add_argument("contrast",
                                  help='String to define the contrast matrix for the GLM, e.g. group. Please note '
                                       'that, when you want to perform negative correlation, the sign is ignored '
                                       'by the command line.')
        clinica_comp.add_argument("string_format",
                                  help='String to define the format of the columns in the TSV file, e.g., '
                                       '%%s %%s %%s %%f if the columns contain a string, a string, a string '
                                       'and a number, respectively.')
        clinica_comp.add_argument("group_id",
                                  help='User-defined identifier for the provided group of subjects.')
        clinica_comp.add_argument("glm_type",
                                  help='String based on the GLM type for the hypothesis. You can choose '
                                       'between group_comparison and correlation.')
        # Optional arguments (e.g. FWHM)
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES['OPTIONAL'])
        optional.add_argument("-fwhm", "--full_width_at_half_maximum",
                              type=int, default=20,
                              help='FWHM for the surface smoothing '
                                   '(default: --full_width_at_half_maximum %(default)s).')
        optional.add_argument("-ft", "--feature_type",
                              type=str, default='cortical_thickness',
                              help='Type of surface-based feature: cortical_thickness (from t1-freesurfer pipeline)'
                                   'or pet_fdg_projection (from pet-surface pipeline)'
                                   '(default: --feature_type %(default)s).')
        # Clinica standard arguments (e.g. --n_procs)
        self.add_clinica_standard_arguments(add_tsv_flag=False)
        # Advanced arguments (i.e. tricky parameters)
        advanced = self._args.add_argument_group(PIPELINE_CATEGORIES['ADVANCED'])
        advanced.add_argument("-cf", "--custom_file",
                              type=str, default=None,
                              help='Pattern of file inside CAPS directory using @subject, @session, @fwhm, @hemi. '
                                   'No --feature_type must be specified in order to use this flag. '
                                   'If you use this flag, you must specify a label with the --feature_label flag). '
                                   'See Wiki for an example.')
        advanced.add_argument("-fl", "--feature_label",
                              type=str, default=None,
                              help='Name of the feature type, it will be saved on the CAPS _measure-FEATURE_LABEL '
                                   'key-value association.')
        advanced.add_argument("-tup", "--threshold_uncorrected_pvalue",
                              type=float, default=0.001,
                              help='Threshold to display the uncorrected p-value '
                                   '(--threshold_uncorrected_pvalue %(default)s).')
        advanced.add_argument("-tcp", "--threshold_corrected_pvalue",
                              type=float, default=0.05,
                              help='Threshold to display the corrected p-value '
                                   '(default: --threshold_corrected_pvalue %(default)s)')
        advanced.add_argument("-ct", "--cluster_threshold",
                              type=float, default=0.001,
                              help='Threshold to define a cluster in the process of cluster-wise correction '
                                   '(default: --cluster_threshold %(default)s).')

    def run_command(self, args):
        """Run the pipeline with defined args."""
        import os
        from networkx import Graph
        from .statistics_surface_pipeline import StatisticsSurface
        from .statistics_surface_utils import (check_inputs,
                                               get_t1_freesurfer_custom_file,
                                               get_fdg_pet_surface_custom_file)
        from clinica.utils.stream import cprint
        from clinica.utils.ux import print_end_pipeline, print_crash_files_and_exit
        from clinica.utils.exceptions import ClinicaException

        if args.feature_type is not None:
            if args.custom_file is not None:
                raise ClinicaException('--feature_type and --custom_file are mutually exclusive: you must choose '
                                       'between one or the other. See documentation for more information.')
            if args.feature_label is not None:
                raise ClinicaException('--feature_label should not be used with --feature_type.')
            # FreeSurfer cortical thickness
            if args.feature_type == 'cortical_thickness':
                args.custom_file = get_t1_freesurfer_custom_file()
                args.feature_label = 'ct'
            # PET cortical projection
            elif args.feature_type == 'pet_fdg_projection':
                args.custom_file = get_fdg_pet_surface_custom_file()
                args.feature_label = 'fdg'
            else:
                raise ClinicaException('Feature type ' + args.feature_type + ' not recognized. Use --custom_file '
                                       'to specify your own files (without --feature_type).')
        elif args.feature_type is None:
            if args.custom_file is None:
                cprint('No feature type selected: using cortical thickness as default value')
                args.custom_file = get_t1_freesurfer_custom_file()
                args.feature_label = 'ct'
            else:
                cprint('Using custom features.')
                if args.feature_label is None:
                    raise ClinicaException('You must specify a --feature_label when using the --custom_files flag.')

        # Check if the group label has been existed, if yes, give an error to the users
        # Note(AR): if the user wants to compare Cortical Thickness measure with PET measure
        # using the group_id, Clinica won't allow it.
        # TODO: Modify this behaviour
        if os.path.exists(os.path.join(os.path.abspath(self.absolute_path(args.caps_directory)), 'groups', 'group-' + args.group_id)):
            error_message = 'group_id: ' + args.group_id + ' already exists, please choose another one or delete ' \
                            'the existing folder and also the working directory and rerun the pipeline'
            raise ClinicaException(error_message)
        parameters = {
            'group_label': args.group_id,
            'design_matrix': args.design_matrix,
            'contrast': args.contrast,
            'str_format': args.string_format,
            'glm_type': args.glm_type,
            'custom_file': args.custom_file,
            'feature_label': args.feature_label,
            'full_width_at_half_maximum': args.full_width_at_half_maximum,
            'threshold_uncorrected_pvalue': args.threshold_uncorrected_pvalue,
            'threshold_corrected_pvalue': args.threshold_corrected_pvalue,
            'cluster_threshold': args.cluster_threshold,
        }
        pipeline = StatisticsSurface(
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subject_visits_with_covariates_tsv),
            base_dir=self.absolute_path(args.working_directory),
            parameters=parameters,
            name=self.name
        )

        check_inputs(pipeline.caps_directory,
                     pipeline.parameters['custom_file'],
                     pipeline.parameters['full_width_at_half_maximum'],
                     pipeline.tsv_file)

        cprint("Parameters used for this pipeline:")
        cprint(pipeline.parameters)

        if args.n_procs:
            exec_pipeline = pipeline.run(plugin='MultiProc',
                                         plugin_args={'n_procs': args.n_procs})
        else:
            exec_pipeline = pipeline.run()

        if isinstance(exec_pipeline, Graph):
            print_end_pipeline(self.name, pipeline.base_dir, pipeline.base_dir_was_specified)
        else:
            print_crash_files_and_exit(args.logname, pipeline.base_dir)
