
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

        return []

    def build_input_node(self):
        """Build and connect an input node to the pipeline.
        """

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from clinica.utils.inputs import clinica_file_reader
        from clinica.utils.exceptions import ClinicaException
        from colorama import Fore

        all_errors = []
        if self.parameters['file_id'] == 'fdg-pet':
            try:
                input_files = clinica_file_reader(self.subjects,
                                                  self.sessions,
                                                  self.caps_directory,
                                                  {'pattern': '*_pet_space-Ixi549Space_suvr-pons_mask-brain_fwhm-8mm_pet.nii*',
                                                   'description': 'pons normalized FDG PET image in MNI space (brain masked)',
                                                   'needed_pipeline': 'pet-volume'})
            except ClinicaException as e:
                all_errors.append(e)
        elif self.parameters['file_id'] == 't1':
            try:
                input_files = clinica_file_reader(self.subjects,
                                                  self.sessions,
                                                  self.caps_directory,
                                                  {'pattern': 't1/spm/segmentation/normalized_space/*_T1w_segm-graymatter_space-Ixi549Space_modulated-off_probability.nii*',
                                                   'description': 'probability map of gray matter segmentation based on T1w image in MNI space',
                                                   'needed_pipeline': 't1-volume or t1-volume-existing-template'})
            except ClinicaException as e:
                all_errors.append(e)
        else:
            raise ClinicaException(Fore.RED + '[Error] ' + Fore.YELLOW + self.parameters['file_id']
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

        # In the same idea as the input node, this output node is supposedly
        # used to write the output fields in a CAPS. It should be executed only
        # if this pipeline output is not already connected to a next Clinica
        # pipeline.

        pass

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline.
        """

        import clinica.pipelines.statistics_volume.statistics_volume_utils as utils
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from clinica.utils.filemanip import unzip_nii
        from os.path import join, dirname

        unzip_node = npe.Node(nutil.Function(input_names=['in_file'],
                                             output_names=['output_files'],
                                             function=unzip_nii),
                              name='unzip_node')

        get_groups = npe.Node(nutil.Function(input_names=['csv', 'contrast'],
                                             output_names=['idx_group1', 'idx_group2'],
                                             function=utils.get_group_1_and_2),
                              name='get_groups')

        get_groups.inputs.contrast = self.parameters['contrast']
        get_groups.inputs.csv = self.tsv_file

        # Script for model creation of 2 sample t tests
        # - scans1 (list of file for first group)
        # - scans2 (list of file for second group)
        # - covariables with their name, and values (first group and second group concatenated)
        # - output directory

        model_creation = npe.Node(nutil.Function(input_names=['csv',
                                                              'contrast',
                                                              'idx_group1',
                                                              'idx_group2',
                                                              'file_list',
                                                              'template_file'],
                                                 output_names=['script_file'],
                                                 function=utils.model_creation),
                                  name='model_creation')
        model_creation.inputs.csv = self.tsv_file
        model_creation.inputs.contrast = self.parameters['contrast']
        model_creation.inputs.template_file = join(dirname(__file__), 'template_model_creation.m')

        run_spm_script_node = npe.Node(nutil.Function(input_names=['m_file'],
                                                      output_names=['spm_mat'],
                                                      function=utils.run_m_script),
                                       name='run_spm_script_node')

        run_spm_model_creation = run_spm_script_node.clone(name='run_spm_model_creation')

        model_estimation = npe.Node(nutil.Function(input_names=['mat_file', 'template_file'],
                                                   output_names=['script_file'],
                                                   function=utils.estimate),
                                    name='model_estimation')
        model_estimation.inputs.template_file = join(dirname(__file__), 'template_model_estimation.m')
        run_spm_model_estimation = run_spm_script_node.clone(name='run_spm_model_estimation')

        model_result = npe.Node(nutil.Function(input_names=['mat_file', 'template_file'],
                                               output_names=['script_file'],
                                               function=utils.results),
                                name='model_result')
        model_result.inputs.template_file = join(dirname(__file__), 'template_model_results.m')
        run_spm_model_result = run_spm_script_node.clone(name='run_spm_model_result')

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
            (run_spm_model_estimation, model_result, [('spm_mat', 'mat_file')]),
            (model_result, run_spm_model_result, [('script_file', 'm_file')])
        ])
