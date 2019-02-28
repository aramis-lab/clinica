# coding: utf8

import clinica.engine as ce


class StatisticsSurfaceCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipelines.
        """
        self._name = 'statistics-surface'

    def define_description(self):
        """Define a description of this pipeline.
        """
        self._description = 'Surface-based mass-univariate analysis with SurfStat:\nhttp://clinica.run/doc/Pipelines/Stats_Surface/'  # noqa

    def define_options(self):
        """Define the sub-command arguments
        """
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES
        from colorama import Fore
        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_id)
        clinica_comp = self._args.add_argument_group('%sMandatory arguments%s' % (Fore.BLUE , Fore.RESET))
        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')
        clinica_comp.add_argument("subject_visits_with_covariates_tsv",
                                  help='TSV file containing a list of subjects with their sessions and all the covariates and factors needed for the GLM.')  # noqa
        clinica_comp.add_argument("design_matrix",
                                  help='String to define the design matrix that fits into the GLM, e.g. 1 + group + sex + age.')
        clinica_comp.add_argument("contrast",
                                  help='String to define the contrast matrix for the GLM, e.g. group. Please note that, when you want to perform negative correlation, the sign is ignored by the command line.')  # noqa
        clinica_comp.add_argument("string_format",
                                  help='String to define the format of the columns in the TSV file, e.g., %%s %%s %%s %%f if the columns contain a string, a string, a string and a number, respectively.')  # noqa
        clinica_comp.add_argument("group_id",
                                  help='User-defined identifier for the provided group of subjects.')
        clinica_comp.add_argument("glm_type",
                                  help='String based on the GLM type for the hypothesis. You can choose between group_comparison and correlation.')  # noqa
        # Optional arguments (e.g. FWHM)
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES['OPTIONAL'])
        optional.add_argument("-fwhm", "--full_width_at_half_maximum",
                              type=int, default=20,
                              help='FWHM for the surface smoothing (default: --full_width_at_half_maximum 20).')
        optional.add_argument("-ft", "--feature_type",
                              type=str, default=None,
                              help='Type of surface-based feature: cortical_thickness or pet_fdg_projection (default: --feature_type cortical_thickness). Note: noddi_projection_ndi, noddi_projection_odi, noddi_projection_fiso, dti_projection_fa, dti_projection_md, dti_projection_rd or dti_projection_ad exists but are not yet officially released.')  # noqa
        # Clinica standard arguments (e.g. --n_procs)
        clinica_opt = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_OPTIONAL'])
        clinica_opt.add_argument("-wd", "--working_directory",
                                 help='Temporary directory to store pipelines intermediate results')
        clinica_opt.add_argument("-np", "--n_procs",
                                 metavar=('N'), type=int,
                                 help='Number of cores used to run in parallel')
        # Advanced arguments (i.e. tricky parameters)
        advanced = self._args.add_argument_group(PIPELINE_CATEGORIES['ADVANCED'])
        advanced.add_argument("-cf", "--custom_file",
                              type=str, default=None,
                              help='Pattern of file inside CAPS directory using @subject, @session, @fwhm, @hemi. No --feature_type must be specified in order to use this flag. If you use this flag, you must specify a label with the --feature_label flag). See Wiki for an example.')  # noqa
        advanced.add_argument("-fl", "--feature_label",
                              type=str, default=None,
                              help='Name of the feature type, it will be saved on the CAPS _measure-FEATURE_LABEL key-value association.')  # noqa
        advanced.add_argument("-tup", "--threshold_uncorrected_pvalue",
                              type=float, default=0.001,
                              help='Threshold to display the uncorrected p-value (--threshold_uncorrected_pvalue 0.001).')  # noqa
        advanced.add_argument("-tcp", "--threshold_corrected_pvalue",
                              type=float, default=0.05,
                              help='Threshold to display the corrected p-value (default: --threshold_corrected_pvalue 0.05)')  # noqa
        advanced.add_argument("-ct", "--cluster_threshold",
                              type=float, default=0.001,
                              help='Threshold to define a cluster in the process of cluster-wise correction (default: --cluster_threshold 0.001).')  # noqa

    def run_command(self, args):
        """
        Run the pipelines with defined args
        """
        from clinica.pipelines.statistics_surface.statistics_surface_pipeline import StatisticsSurface
        from clinica.pipelines.statistics_surface.statistics_surface_utils import check_inputs
        from clinica.utils.stream import cprint
        import os

        if args.feature_type is not None:
            if args.custom_file is not None:
                raise Exception('--feature_type and --custom_file are mutually exclusive : you must choose between one or the other. See documentation for more informations.')
            if args.feature_label is not None:
                raise Exception('--feature_label should not be used with --feature_type')
            # FreeSurfer cortical thickness
            if args.feature_type == 'cortical_thickness':
                args.custom_file = '@subject/@session/t1/freesurfer_cross_sectional/@subject_@session/surf/@hemi.thickness.fwhm@fwhm.fsaverage.mgh'
                args.feature_label = 'ct'
            # PET cortical projection
            elif args.feature_type == 'pet_fdg_projection':
                args.custom_file = '@subject/@session/pet/surface/@subject_@session_task-rest_acq-fdg_pet_space-fsaverage_suvr-pons_pvc-iy_hemi-@hemi_fwhm-@fwhm_projection.mgh'
                args.feature_label = 'fdg'
            # NODDI, NDI, ODI and FISO
            elif args.feature_type == 'noddi_projection_ndi':
                args.custom_file = '@subject/@session/noddi/postprocessing/noddi-register-vertex-fsaverage/cortex-projection/@subject_@session_OnFsaverage_fwhm-@fwhm_measure-ficvf_hemi-@hemi.mgh'
                args.feature_label = 'NDI'
            elif args.feature_type == 'noddi_projection_fiso':
                args.custom_file = '@subject/@session/noddi/postprocessing/noddi-register-vertex-fsaverage/cortex-projection/@subject_@session_OnFsaverage_fwhm-@fwhm_measure-fiso_hemi-@hemi.mgh'
                args.feature_label = 'FISO'
            elif args.feature_type == 'noddi_projection_odi':
                args.custom_file = '@subject/@session/noddi/postprocessing/noddi-register-vertex-fsaverage/cortex-projection/@subject_@session_OnFsaverage_fwhm-@fwhm_measure-odi_hemi-@hemi.mgh'
                args.feature_label = 'ODI'
            # DTI fa, md, rd and ad
            elif args.feature_type == 'dti_projection_fa':
                args.custom_file = '@subject/@session/noddi/postprocessing/dti-register-vertex-fsaverage/cortex-projection/@subject_@session_OnFsaverage_fwhm-@fwhm_measure-fa_hemi-@hemi.mgh'
                args.feature_label = 'FA'
            elif args.feature_type == 'dti_projection_md':
                args.custom_file = '@subject/@session/noddi/postprocessing/dti-register-vertex-fsaverage/cortex-projection/@subject_@session_OnFsaverage_fwhm-@fwhm_measure-md_hemi-@hemi.mgh'
                args.feature_label = 'MD'
            elif args.feature_type == 'dti_projection_rd':
                args.custom_file = '@subject/@session/noddi/postprocessing/dti-register-vertex-fsaverage/cortex-projection/@subject_@session_OnFsaverage_fwhm-@fwhm_measure-rd_hemi-@hemi.mgh'
                args.feature_label = 'RD'
            elif args.feature_type == 'dti_projection_ad':
                args.custom_file = '@subject/@session/noddi/postprocessing/dti-register-vertex-fsaverage/cortex-projection/@subject_@session_OnFsaverage_fwhm-@fwhm_measure-ad_hemi-@hemi.mgh'
                args.feature_label = 'AD'
            else:
                raise Exception('Feature type ' + args.feature_type + ' not recognized. Use --custom_file to specify your own files (without --feature_type).')
        elif args.feature_type is None:
            if args.custom_file is None:
                cprint('No feature type selected : using cortical thickness as default value')
                args.custom_file = '@subject/@session/t1/freesurfer_cross_sectional/@subject_@session/surf/@hemi.thickness.fwhm@fwhm.fsaverage.mgh'
                args.feature_label = 'ct'
            else:
                cprint('Using custom features.')
                if args.feature_label is None:
                    raise Exception('You must specify a --feature_label when using the --custom_files flag')

        # Check if the group label has been existed, if yes, give an error to the users
        if os.path.exists(os.path.join(os.path.abspath(self.absolute_path(args.caps_directory)), 'groups', 'group-' + args.group_id)):
            error_message = 'group_id : ' + args.group_id + ' already exists, please choose another one or delete the existing folder and also the working directory and rerun the pipelines'
            raise Exception(error_message)

        pipeline = StatisticsSurface(
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subject_visits_with_covariates_tsv))
        pipeline.parameters = {
            # pass these args by using self.parameters in a dictionary
            'design_matrix': args.design_matrix,
            'contrast': args.contrast,
            'str_format': args.string_format,
            'group_label': args.group_id,
            'glm_type': args.glm_type,
            'custom_file': args.custom_file,
            'feature_label': args.feature_label,
            'full_width_at_half_maximum': args.full_width_at_half_maximum,
            'threshold_uncorrected_pvalue': args.threshold_uncorrected_pvalue,
            'threshold_corrected_pvalue': args.threshold_corrected_pvalue,
            'cluster_threshold': args.cluster_threshold
        }
        pipeline.base_dir = self.absolute_path(args.working_directory)

        check_inputs(pipeline.caps_directory,
                     pipeline.parameters['custom_file'],
                     pipeline.parameters['full_width_at_half_maximum'],
                     pipeline.tsv_file)

        if args.n_procs:
            pipeline.run(plugin='MultiProc',
                         plugin_args={'n_procs': args.n_procs})
        else:
            print(pipeline.parameters)
            pipeline.run()

        cprint("The " + self._name + " pipeline has completed. You can now delete the working directory (" + args.working_directory + ").")
