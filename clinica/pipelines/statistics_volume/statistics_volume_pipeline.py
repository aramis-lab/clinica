# coding: utf-8

import clinica.pipelines.engine as cpe


class StatisticsVolume(cpe.Pipeline):
    """StatisticsVolume - Volume-based mass-univariate analysis with SPM.

    Returns:
        A clinica pipeline object containing the StatisticsVolume pipeline.
    """
    def check_pipeline_parameters(self):
        """Check pipeline parameters."""
        from clinica.utils.exceptions import ClinicaException

        if self.parameters['cluster_threshold'] < 0 or self.parameters['cluster_threshold'] > 1:
            raise ClinicaException("Cluster threshold should be between 0 and 1 "
                                   "(given value: %s)." % self.parameters['cluster_threshold'])

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file."""
        pass

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Returns:
            A list of (string) input fields name.
        """

        return ['input_files']

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """

        return ['spmT_0001',
                'spmT_0002',
                'spm_figures',
                'variance_of_error',
                'resels_per_voxels',
                'mask',
                'regression_coeff',
                'contrast']

    def build_input_node(self):
        """Build and connect an input node to the pipeline."""
        import os
        from colorama import Fore
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from clinica.iotools.utils.data_handling import check_volume_location_in_world_coordinate_system
        from clinica.utils.exceptions import ClinicaBIDSError, ClinicaException
        from clinica.utils.filemanip import extract_subjects_sessions_from_filename, save_participants_sessions
        from clinica.utils.inputs import clinica_file_reader
        from clinica.utils.input_files import T1W_NII
        from clinica.utils.stream import cprint
        from clinica.utils.ux import print_images_to_process, print_begin_image

        gic = '*'
        if self.parameters['group_id_caps'] is not None:
            gic = self.parameters['group_id_caps']

        all_errors = []
        if self.parameters['feature_type'] == 'fdg':
            try:
                input_files = clinica_file_reader(self.subjects,
                                                  self.sessions,
                                                  self.caps_directory,
                                                  {'pattern': '*_pet_space-Ixi549Space_suvr-pons_mask-brain_fwhm-' + str(self.parameters['full_width_at_half_maximum']) + 'mm_pet.nii*',
                                                   'description': 'pons normalized FDG PET image in MNI space (brain masked)',
                                                   'needed_pipeline': 'pet-volume'})
            except ClinicaException as e:
                all_errors.append(e)
        elif self.parameters['feature_type'] == 'graymatter':
            try:
                input_files = clinica_file_reader(self.subjects,
                                                  self.sessions,
                                                  self.caps_directory,
                                                  {'pattern': 't1/spm/dartel/group-' + gic + '/*_T1w_segm-graymatter_space-Ixi549Space_modulated-on_fwhm-' + str(self.parameters['full_width_at_half_maximum']) + 'mm_probability.nii.*',
                                                   'description': 'probability map of gray matter segmentation based on T1w image in MNI space',
                                                   'needed_pipeline': 't1-volume or t1-volume-existing-template'})
            except ClinicaException as e:
                all_errors.append(e)

        else:
            if not self.parameters['custom_files']:
                raise ClinicaException(Fore.RED + '[Error] You did not specify the --custom_files flag in the command line for the feature type '
                                       + Fore.Blue + self.parameters['feature_type'] + Fore.RED + '! Clinica can\'t '
                                       + 'know what file to use in your analysis ! Type: \n\t' + Fore.BLUE + 'clinica run statistics-volume\n'
                                       + Fore.RED + ' to have help on how to use the command line.' + Fore.RESET)
            try:
                # If custom file are grabbed, information of fwhm is irrelevant and should not appear on final filenames
                self.parameters['full_width_at_half_maximum'] = None
                input_files = clinica_file_reader(self.subjects,
                                                  self.sessions,
                                                  self.caps_directory,
                                                  {'pattern': self.parameters['custom_files'],
                                                   'description': 'custom file provided by user'})
            except ClinicaException as e:
                all_errors.append(e)

        if len(all_errors) > 0:
            error_message = 'Clinica faced errors while trying to read files in your BIDS or CAPS directories.\n'
            for msg in all_errors:
                error_message += str(msg)
            raise ClinicaException(error_message)

        read_parameters_node = npe.Node(name="LoadingCLIArguments",
                                        interface=nutil.IdentityInterface(
                                            fields=self.get_input_fields(),
                                            mandatory_inputs=True),
                                        synchronize=True)
        read_parameters_node.inputs.input_files = input_files

        if len(self.subjects):
            print_images_to_process(self.subjects, self.sessions)
            cprint('The pipeline will last a few minutes. Images generated by SPM will popup during the pipeline.')
            print_begin_image('group-' + self.parameters['group_id'])

        self.connect([
            (read_parameters_node,      self.input_node,    [('input_files',    'input_files')])
        ])

    def build_output_node(self):
        """Build and connect an output node to the pipeline."""
        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio
        from os.path import join, pardir

        relative_path = join('groups',
                             'group-' + self.parameters['group_id'],
                             'statistics_volume',
                             'group_comparison_measure-' + self.parameters['feature_type'])

        datasink = npe.Node(nio.DataSink(), name='sinker')
        datasink.inputs.base_directory = join(self.caps_directory, relative_path)

        datasink.inputs.parameterization = True
        if self.parameters['full_width_at_half_maximum']:
            datasink.inputs.regexp_substitutions = [
                # t-stat map
                (join(self.caps_directory, relative_path) + r'/spm_results_analysis_./(.*)',
                 join(self.caps_directory, relative_path) + r'/\1'),

                # contrasts
                (join(self.caps_directory, relative_path) + r'/contrasts/(.*)',
                 join(self.caps_directory, relative_path) + r'/\1'),

                # resels per voxels
                (join(self.caps_directory, relative_path) + '/resels_per_voxels/resels_per_voxel.nii',
                 join(self.caps_directory, relative_path) + '/group-' + self.parameters['group_id'] + '_RPV.nii'),

                # mask
                (join(self.caps_directory, relative_path) + '/mask/included_voxel_mask.nii',
                 join(self.caps_directory, relative_path) + '/group-' + self.parameters['group_id'] + '_mask.nii'),

                # variance of error
                (join(self.caps_directory, relative_path) + '/variance_of_error/(.*)',
                 join(self.caps_directory, relative_path) + r'/\1'),

                # tsv file
                (join(self.caps_directory, relative_path) + r'/tsv_file/.*',
                 join(self.caps_directory, relative_path) + '/' + pardir + '/group-' + self.parameters['group_id'] + '_participants.tsv'),

                # report (figures)
                (join(self.caps_directory, relative_path) + r'/figures/(.*)',
                 join(self.caps_directory, relative_path) + r'/\1'),

                # regression coefficient
                (join(self.caps_directory, relative_path) + r'/regression_coeff/(.*).nii',
                 join(self.caps_directory, relative_path) + '/group-' + self.parameters['group_id'] + r'_covariate-\1'
                 + '_measure-' + self.parameters['feature_type'] + '_fwhm-'
                 + str(self.parameters['full_width_at_half_maximum']) + '_regressionCoefficient.nii'),

            ]

        datasink.inputs.tsv_file = self.tsv_file

        self.connect([
            (self.output_node, datasink, [('spmT_0001', 'spm_results_analysis_1')]),
            (self.output_node, datasink, [('spmT_0002', 'spm_results_analysis_2')]),
            (self.output_node, datasink, [('spm_figures', 'figures')]),
            (self.output_node, datasink, [('variance_of_error', 'variance_of_error')]),
            (self.output_node, datasink, [('resels_per_voxels', 'resels_per_voxels')]),
            (self.output_node, datasink, [('mask', 'mask')]),
            (self.output_node, datasink, [('regression_coeff', 'regression_coeff')]),
            (self.output_node, datasink, [('contrasts', 'contrasts')])
        ])

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline."""
        import clinica.pipelines.statistics_volume.statistics_volume_utils as utils
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from clinica.utils.filemanip import unzip_nii
        from os.path import join, dirname

        # SPM cannot handle zipped files
        unzip_node = npe.Node(nutil.Function(input_names=['in_file'],
                                             output_names=['output_files'],
                                             function=unzip_nii),
                              name='unzip_node')

        # Get indexes of the 2 groups, based on the contrast column of the tsv file
        get_groups = npe.Node(nutil.Function(input_names=['tsv', 'contrast'],
                                             output_names=['idx_group1', 'idx_group2', 'class_names'],
                                             function=utils.get_group_1_and_2),
                              name='get_groups')

        get_groups.inputs.contrast = self.parameters['contrast']
        get_groups.inputs.tsv = self.tsv_file

        # Run SPM nodes are all a copy of a generic SPM script launcher
        run_spm_script_node = npe.Node(nutil.Function(input_names=['m_file'],
                                                      output_names=['spm_mat'],
                                                      function=utils.run_m_script),
                                       name='run_spm_script_node')

        run_spm_model_creation = run_spm_script_node.clone(name='run_spm_model_creation')
        run_spm_model_estimation = run_spm_script_node.clone(name='run_spm_model_estimation')
        run_spm_model_contrast = run_spm_script_node.clone(name='run_spm_model_contrast')
        run_spm_model_result_no_correction = run_spm_script_node.clone(name='run_spm_model_result_no_correction')

        # All the following node are creating the correct (.m) script for the different SPM steps
        # 1. Model creation
        # 2. Model estimation
        # 3. Creation of contrast with covariates
        # 4. Creation of results

        # 1. Model creation
        # We use overwrite option to be sure this node is always run so that it can delete the output dir if it
        # already exists (this may cause error in output files otherwise)
        model_creation = npe.Node(nutil.Function(input_names=['tsv',
                                                              'contrast',
                                                              'idx_group1',
                                                              'idx_group2',
                                                              'file_list',
                                                              'template_file'],
                                                 output_names=['script_file', 'covariates'],
                                                 function=utils.model_creation),
                                  name='model_creation',
                                  overwrite=True)
        model_creation.inputs.tsv = self.tsv_file
        model_creation.inputs.contrast = self.parameters['contrast']
        model_creation.inputs.template_file = join(dirname(__file__), 'template_model_creation.m')

        # 2. Model estimation
        model_estimation = npe.Node(nutil.Function(input_names=['mat_file', 'template_file'],
                                                   output_names=['script_file'],
                                                   function=utils.estimate),
                                    name='model_estimation')
        model_estimation.inputs.template_file = join(dirname(__file__), 'template_model_estimation.m')

        # 3. Contrast
        model_contrast = npe.Node(nutil.Function(input_names=['mat_file',
                                                              'template_file',
                                                              'covariates',
                                                              'class_names'],
                                                 output_names=['script_file'],
                                                 function=utils.contrast),
                                  name='model_contrast')
        model_contrast.inputs.template_file = join(dirname(__file__), 'template_model_contrast.m')

        # 4. Results
        model_result_no_correction = npe.Node(nutil.Function(input_names=['mat_file', 'template_file',  'method', 'threshold'],
                                                             output_names=['script_file'],
                                                             function=utils.results),
                                              name='model_result_no_correction')
        model_result_no_correction.inputs.template_file = join(dirname(__file__), 'template_model_results.m')

        model_result_no_correction.inputs.method = 'none'
        model_result_no_correction.inputs.threshold = self.parameters['cluster_threshold']

        # Print result to txt file if spm

        # Export results to output node
        read_output_node = npe.Node(nutil.Function(input_names=['spm_mat', 'class_names', 'covariates', 'group_id', 'fwhm', 'measure'],
                                                   output_names=['spmT_0001',
                                                                 'spmT_0002',
                                                                 'spm_figures',
                                                                 'variance_of_error',
                                                                 'resels_per_voxels',
                                                                 'mask',
                                                                 'regression_coeff',
                                                                 'contrasts'],
                                                   function=utils.read_output),
                                    name='read_output_node')
        read_output_node.inputs.group_id = self.parameters['group_id']
        read_output_node.inputs.fwhm = self.parameters['full_width_at_half_maximum']
        read_output_node.inputs.measure = self.parameters['feature_type']

        # Connection
        # ==========
        self.connect([
            (self.input_node, unzip_node, [('input_files', 'in_file')]),
            (unzip_node, model_creation, [('output_files', 'file_list')]),
            (get_groups, model_creation, [('idx_group1', 'idx_group1')]),
            (get_groups, model_creation, [('idx_group2', 'idx_group2')]),

            (model_creation, run_spm_model_creation, [('script_file', 'm_file')]),
            (run_spm_model_creation, model_estimation, [('spm_mat', 'mat_file')]),
            (model_estimation, run_spm_model_estimation, [('script_file', 'm_file')]),

            (get_groups, model_contrast, [('class_names', 'class_names')]),
            (run_spm_model_estimation, model_contrast, [('spm_mat', 'mat_file')]),
            (model_creation, model_contrast, [('covariates', 'covariates')]),

            (model_contrast, run_spm_model_contrast, [('script_file', 'm_file')]),

            (run_spm_model_contrast, model_result_no_correction, [('spm_mat', 'mat_file')]),

            (model_result_no_correction, run_spm_model_result_no_correction, [('script_file', 'm_file')]),

            (run_spm_model_result_no_correction, read_output_node, [('spm_mat', 'spm_mat')]),

            (get_groups, read_output_node, [('class_names', 'class_names')]),
            (model_creation, read_output_node, [('covariates', 'covariates')]),

            (read_output_node, self.output_node, [('spmT_0001', 'spmT_0001')]),
            (read_output_node, self.output_node, [('spmT_0002', 'spmT_0002')]),
            (read_output_node, self.output_node, [('spm_figures', 'spm_figures')]),
            (read_output_node, self.output_node, [('variance_of_error', 'variance_of_error')]),
            (read_output_node, self.output_node, [('resels_per_voxels', 'resels_per_voxels')]),
            (read_output_node, self.output_node, [('mask', 'mask')]),
            (read_output_node, self.output_node, [('regression_coeff', 'regression_coeff')]),
            (read_output_node, self.output_node, [('contrasts', 'contrasts')]),
        ])
