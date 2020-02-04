
__author__ = "Arnaud Marcoux"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Arnaud Marcoux"]
__license__ = "See LICENSE.txt file"
__version__ = "0.3.0"
__maintainer__ = "Arnaud Marcoux"
__email__ = "arnaud.marcoux@icm-institute.org"
__status__ = "Development"


import clinica.pipelines.engine as cpe


class StatisticsVolume(cpe.Pipeline):
    """Statistics_Volume SHORT DESCRIPTION.


    Args:
        input_dir: A BIDS directory.
        output_dir: An empty output directory where CAPS structured data will be written.
        subjects_sessions_list: The Subjects-Sessions list file (in .tsv format).

    Returns:
        A clinica pipeline object containing the Statistics_Volume pipeline.

    Raises:


    """

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file.
        """
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
                'new_figure_names',
                'variance_of_error',
                'resels_per_voxels',
                'mask',
                'regression_coeff',
                'contrast ']

    def build_input_node(self):
        """Build and connect an input node to the pipeline.
        """

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from clinica.utils.inputs import clinica_file_reader
        from clinica.utils.exceptions import ClinicaException
        from colorama import Fore

        gic = '*'
        if self.parameters['group_id_caps'] is not None:
            gic = self.parameters['group_id_caps']

        all_errors = []
        if self.parameters['feature_type'] == 'fdg-pet':
            try:
                input_files = clinica_file_reader(self.subjects,
                                                  self.sessions,
                                                  self.caps_directory,
                                                  {'pattern': '*_pet_space-Ixi549Space_suvr-pons_mask-brain_fwhm-' + str(self.parameters['smoothing']) + 'mm_pet.nii*',
                                                   'description': 'pons normalized FDG PET image in MNI space (brain masked)',
                                                   'needed_pipeline': 'pet-volume'})
            except ClinicaException as e:
                all_errors.append(e)
        elif self.parameters['feature_type'] == 't1-gm':
            try:
                input_files = clinica_file_reader(self.subjects,
                                                  self.sessions,
                                                  self.caps_directory,
                                                  {'pattern': 't1/spm/dartel/group-' + gic + '/*_T1w_segm-graymatter_space-Ixi549Space_modulated-on_fwhm-' + str(self.parameters['smoothing']) + 'mm_probability.nii.*',
                                                   'description': 'probability map of gray matter segmentation based on T1w image in MNI space',
                                                   'needed_pipeline': 't1-volume or t1-volume-existing-template'})
            except ClinicaException as e:
                all_errors.append(e)

        elif self.parameters['feature_type'] == 't1-wm':
            try:
                input_files = clinica_file_reader(self.subjects,
                                                  self.sessions,
                                                  self.caps_directory,
                                                  {'pattern': 't1/spm/dartel/group-' + gic + '/*_T1w_segm-whitematter_space-Ixi549Space_modulated-on_fwhm-' + str(self.parameters['smoothing']) + 'mm_probability.nii.*',
                                                   'description': 'probability map of white matter segmentation based on T1w image in MNI space',
                                                   'needed_pipeline': 't1-volume or t1-volume-existing-template'})
            except ClinicaException as e:
                all_errors.append(e)
        elif self.parameters['feature_type'] == 'custom':
            if not self.parameters['custom_files']:
                raise ClinicaException(Fore.RED + '[Error] You did not specify the --custom_files flag in the command line ! Clinica can\'t '
                                       + 'know what file to use in your analysis ! Type: \n\t' + Fore.BLUE + 'clinica run statistics-volume\n'
                                       + Fore.RED + ' to have help on how to use the command line.' + Fore.RESET)
            try:
                input_files = clinica_file_reader(self.subjects,
                                                  self.sessions,
                                                  self.caps_directory,
                                                  {'pattern': self.parameters['custom_files'],
                                                   'description': 'custom file provided by user'})
            except ClinicaException as e:
                all_errors.append(e)
        else:
            raise ClinicaException(Fore.RED + '[Error] ' + Fore.YELLOW + self.parameters['feature_type']
                                   + Fore.RED + ' modality is not currently supported for this analysis' + Fore.RESET)

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

        self.connect([
            (read_parameters_node,      self.input_node,    [('input_files',    'input_files')])
        ])

    def build_output_node(self):
        """Build and connect an output node to the pipeline.
        """
        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio
        from os.path import join

        datasink = npe.Node(nio.DataSink(),
                            name='sinker')
        datasink.inputs.base_directory = join(self.caps_directory,
                                              'groups',
                                              'group-' + self.parameters['group_id'],
                                              'statistics',
                                              'spm_2_sample_t_test')
        datasink.inputs.parameterization = True
        datasink.inputs.regexp_substitutions = [
            (r'(.*)/group-' + self.parameters['group_id'] + r'/statistics/spm_2_sample_t_test/spm_results_analysis_./(.*)',
             r'\1/group-' + self.parameters['group_id'] + r'/statistics/spm_2_sample_t_test/\2'),
            (r'(.*)/group-' + self.parameters['group_id'] + r'/statistics/spm_2_sample_t_test/tsv_file/.*',
             r'\1/group-' + self.parameters['group_id'] + r'/statistics/participant.tsv')
            # Uncomment if you need all the files in the same folder
            # ,
            # (r'(.*)/group-' + self.parameters['group_id'] + r'.*/(.*)',
            # r'\1/group-' + self.parameters['group_id'] + r'/\2')
        ]

        datasink.inputs.tsv_file = self.tsv_file

        self.connect([
            (self.output_node, datasink, [('spmT_0001', 'spm_results_analysis_1')]),
            (self.output_node, datasink, [('spmT_0002', 'spm_results_analysis_2')]),
            (self.output_node, datasink, [('new_figure_names', 'figures')]),
            (self.output_node, datasink, [('variance_of_error', 'variance_of_error')]),
            (self.output_node, datasink, [('resels_per_voxels', 'resels_per_voxels')]),
            (self.output_node, datasink, [('mask', 'mask')]),
            (self.output_node, datasink, [('regression_coeff', 'regression_coeff')]),
            (self.output_node, datasink, [('contrasts', 'contrasts')])
        ])

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline.
        """

        import clinica.pipelines.statistics_volume.statistics_volume_utils as utils
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from clinica.utils.filemanip import unzip_nii
        from os.path import join, dirname, abspath

        # SPM cannot handle zipped files
        unzip_node = npe.Node(nutil.Function(input_names=['in_file'],
                                             output_names=['output_files'],
                                             function=unzip_nii),
                              name='unzip_node')

        # Get indexes of the 2 groups, based on the contrast column of the tsv file
        get_groups = npe.Node(nutil.Function(input_names=['csv', 'contrast'],
                                             output_names=['idx_group1', 'idx_group2', 'class_names'],
                                             function=utils.get_group_1_and_2),
                              name='get_groups')

        get_groups.inputs.contrast = self.parameters['contrast']
        get_groups.inputs.csv = self.tsv_file

        # Run SPM nodes are all a copy of a generic SPM script launcher

        run_spm_script_node = npe.Node(nutil.Function(input_names=['m_file'],
                                                      output_names=['spm_mat'],
                                                      function=utils.run_m_script),
                                       name='run_spm_script_node')

        run_spm_model_creation = run_spm_script_node.clone(name='run_spm_model_creation')
        run_spm_model_estimation = run_spm_script_node.clone(name='run_spm_model_estimation')
        run_spm_model_contrast = run_spm_script_node.clone(name='run_spm_model_contrast')
        run_spm_model_result_no_correction = run_spm_script_node.clone(name='run_spm_model_result_no_correction')
        run_spm_model_result_FWE_correction = run_spm_script_node.clone(name='run_spm_model_result_FWE_correction')

        # All the following node are creating the correct (.m) script for the different SPM steps
        # 1. Model creation
        # 2. Model estimation
        # 3. Creation of contrast with covariables
        # 4. Creation of results

        # 1. Model creation
        # We use overwritte option to be sure this node is always run so that it can delete the output dir if it
        # already exists (this may cause error in output files otherwise)
        model_creation = npe.Node(nutil.Function(input_names=['csv',
                                                              'contrast',
                                                              'idx_group1',
                                                              'idx_group2',
                                                              'file_list',
                                                              'template_file'],
                                                 output_names=['script_file', 'covariables'],
                                                 function=utils.model_creation),
                                  name='model_creation',
                                  overwrite=True)
        model_creation.inputs.csv = self.tsv_file
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
                                                              'covariables',
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

        model_result_FWE_correction = model_result_no_correction.clone(name='model_result_FWE_correction')

        model_result_no_correction.inputs.method = 'none'
        model_result_no_correction.inputs.threshold = self.parameters['threshold_uncorrected_pvalue']

        model_result_FWE_correction.inputs.method = 'FWE'
        model_result_FWE_correction.inputs.threshold = self.parameters['threshold_corrected_pvalue']


        # Print result to txt file if spm

        # Export results to output node
        read_output_node = npe.Node(nutil.Function(input_names=['spm_mat', 'spm_mat_2', 'class_names', 'covariables'],
                                                   output_names=['spmT_0001',
                                                                 'spmT_0002',
                                                                 'new_figure_names',
                                                                 'variance_of_error',
                                                                 'resels_per_voxels',
                                                                 'mask',
                                                                 'regression_coeff',
                                                                 'contrasts'],
                                                   function=utils.read_output),
                                    name='read_output_node')

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
            (model_creation, model_contrast, [('covariables', 'covariables')]),

            (model_contrast, run_spm_model_contrast, [('script_file', 'm_file')]),

            (run_spm_model_contrast, model_result_no_correction, [('spm_mat', 'mat_file')]),
            (run_spm_model_contrast, model_result_FWE_correction, [('spm_mat', 'mat_file')]),

            (model_result_no_correction, run_spm_model_result_no_correction, [('script_file', 'm_file')]),
            (model_result_FWE_correction, run_spm_model_result_FWE_correction, [('script_file', 'm_file')]),

            (run_spm_model_result_no_correction, read_output_node, [('spm_mat', 'spm_mat')]),
            (run_spm_model_result_FWE_correction, read_output_node, [('spm_mat', 'spm_mat_2')]), # Only used for synchronisation

            (get_groups, read_output_node, [('class_names', 'class_names')]),
            (model_creation, read_output_node, [('covariables', 'covariables')]),

            (read_output_node, self.output_node, [('spmT_0001', 'spmT_0001')]),
            (read_output_node, self.output_node, [('spmT_0002', 'spmT_0002')]),
            (read_output_node, self.output_node, [('new_figure_names', 'new_figure_names')]),
            (read_output_node, self.output_node, [('variance_of_error', 'variance_of_error')]),
            (read_output_node, self.output_node, [('resels_per_voxels', 'resels_per_voxels')]),
            (read_output_node, self.output_node, [('mask', 'mask')]),
            (read_output_node, self.output_node, [('regression_coeff', 'regression_coeff')]),
            (read_output_node, self.output_node, [('contrasts', 'contrasts')]),
        ])
