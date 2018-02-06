# coding: utf8

import clinica.engine as ce

__author__ = "Junhao Wen"
__copyright__ = "Copyright 2016-2018, The Aramis Lab Team"
__credits__ = ["Junhao Wen"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Junhao Wen"
__email__ = "Junhao.Wen@inria.fr"
__status__ = "Development"


class StatisticsSurfaceCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipelines.
        """
        self._name = 'statistics-surface'

    def define_description(self):
        """Define a description of this pipeline.
        """
        self._description = 'Surface-based mass-univariate analysis with SurfStat:\nhttp://clinica.run/doc/Pipelines/Stats_Surface/'

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
                                help='Feature type. Can be : cortical_thickness, pet_fdg_projection, noddi_projection_ndi, noddi_projection_odi, noddi_projection_fiso, dti_projection_fa, dti_projection_md, dti_projection_rd and dti_projection_ad. Default = cortical_thickness')
        self._args.add_argument("-cf", "--custom_file", type=str, default=None,
                                help='Pattern of file inside caps directory using @subject, @session, @fwhm, @hemi. No --feature_type must be specified in order to use this flag.')
        self._args.add_argument("-fwhm", "--full_width_at_half_maximum", type=int, default=20,
                                help='FWHM for the surface smoothing (default=20)')
        self._args.add_argument("-tup", "--threshold_uncorrected_pvalue", type=float, default=0.001,
                                help='Threshold to display the uncorrected pvalue (default=0.001)')
        self._args.add_argument("-tcp", "--threshold_corrected_pvalue", type=float, default=0.05,
                                help='Threshold to display the corrected pvalue (default=0.05)')
        self._args.add_argument("-ct", "--cluster_threshold", type=float, default=0.001,
                                help='Threshold to define a cluster in the process of cluster-wise correction (default=0.001)')
        self._args.add_argument("-np", "--n_procs", type=int, default=4,
                                help='Number of parallel processes to run (default=4)')
        self._args.add_argument("-wd", "--working_directory", type=str, default=None,
                                help='Temporary directory to run the workflow')

    def run_command(self, args):
        """
        Run the pipelines with defined args
        """

        from statistics_surface_pipeline import StatisticsSurface
        from statistics_surface_utils import check_inputs
        from clinica.utils.stream import cprint
        import os

        if args.feature_type is not None:
            if args.custom_file is not None:
                raise Exception('--feature_type and --custom_file are mutually exclusive : you must choose between one or the other. See documentation for more informations.')
            # freesurfer cortical thickness
            if args.feature_type == 'cortical_thickness':
                args.custom_file = '@subject/@session/t1/freesurfer_cross_sectional/@subject_@session/surf/@hemi.thickness.fwhm@fwhm.fsaverage.mgh'
            # pet cortical projection
            elif args.feature_type == 'pet_fdg_projection':
                args.custom_file = '@subject/@session/pet/surface/@subject_@session_task-rest_acq-FDG_pet_space-fsaverage_suvr-pons_pvc-iy_hemi-@hemi_fwhm-@fwhm_projection.mgh'
            # NODDI NDI, ODI and FISO
            elif args.feature_type == 'noddi_projection_ndi':
                args.custom_file = '@subject/@session/noddi/postprocessing/noddi-register-vertex-fsaverage/cortex-projection/@subject_@session_OnFsaverage_fwhm-@fwhm_measure-ficvf_hemi-@hemi.mgh'
            elif args.feature_type == 'noddi_projection_fiso':
                args.custom_file = '@subject/@session/noddi/postprocessing/noddi-register-vertex-fsaverage/cortex-projection/@subject_@session_OnFsaverage_fwhm-@fwhm_measure-fiso_hemi-@hemi.mgh'
            elif args.feature_type == 'noddi_projection_odi':
                args.custom_file = '@subject/@session/noddi/postprocessing/noddi-register-vertex-fsaverage/cortex-projection/@subject_@session_OnFsaverage_fwhm-@fwhm_measure-odi_hemi-@hemi.mgh'
            # DTI fa, md, rd and ad
            elif args.feature_type == 'dti_projection_fa':
                args.custom_file = '@subject/@session/noddi/postprocessing/dti-register-vertex-fsaverage/cortex-projection/@subject_@session_OnFsaverage_fwhm-@fwhm_measure-fa_hemi-@hemi.mgh'
            elif args.feature_type == 'dti_projection_md':
                args.custom_file = '@subject/@session/noddi/postprocessing/dti-register-vertex-fsaverage/cortex-projection/@subject_@session_OnFsaverage_fwhm-@fwhm_measure-md_hemi-@hemi.mgh'
            elif args.feature_type == 'dti_projection_rd':
                args.custom_file = '@subject/@session/noddi/postprocessing/dti-register-vertex-fsaverage/cortex-projection/@subject_@session_OnFsaverage_fwhm-@fwhm_measure-rd_hemi-@hemi.mgh'
            elif args.feature_type == 'dti_projection_ad':
                args.custom_file = '@subject/@session/noddi/postprocessing/dti-register-vertex-fsaverage/cortex-projection/@subject_@session_OnFsaverage_fwhm-@fwhm_measure-ad_hemi-@hemi.mgh'

            else:
                raise Exception('Feature type ' + args.feature_type + ' not recognized. Use --custom_file to specify your own files (without --feature_type).')
        elif args.feature_type is None:
            if args.custom_file is None:
                cprint('No feature type selected : using cortical thickness as default value')
                args.custom_file = '@subject/@session/t1/freesurfer_cross_sectional/@subject_@session/surf/@hemi.thickness.fwhm@fwhm.fsaverage.mgh'
            else:
                cprint('Using custom features.')

        #Check if the group label has been existed, if yes, give the warning to the users
        if os.path.exists(os.path.join(os.path.abspath(self.absolute_path(args.caps_directory)), 'groups', 'group-' + args.group_label)):
            error_message = 'group_id : ' + args.group_label + ' already exists, please choose another one or delete the existing folder and also the working directory and rerun the pipelines'
            raise Exception(error_message)

        pipeline = StatisticsSurface(
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

        # run the pipelines in n_procs cores based on your computation power.
        if args.n_procs:
            pipeline.write_graph()
            pipeline.run(plugin='MultiProc', plugin_args={'n_procs': args.n_procs})
        else:
            pipeline.write_graph()
            print(pipeline.parameters)
            pipeline.run()
