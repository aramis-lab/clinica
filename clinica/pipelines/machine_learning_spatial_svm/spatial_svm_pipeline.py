# coding: utf8

# WARNING: Don't put any import statement here except if it's absolutely
# necessary. Put it *inside* the different methods.
# Otherwise it will slow down the dynamic loading of the pipelines list by the
# command line tool.
import clinica.pipelines.engine as cpe

__author__ = "Simona Bottani"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Simona Bottani"
__email__ = "simona.bottani@icm-institute.org"
__status__ = "Development"


class SpatialSVM(cpe.Pipeline):
    """SpatialSVM - Prepare input data for SVM with spatial and anatomical regularization.

    Todos:
        - [ ] Final version of CAPS.
        - [ ] Remove --voxel_size flag and detect automatically this parameter.

    Args:
        input_dir: A BIDS directory.
        output_dir: An empty output directory where CAPS structured data will be written.
        subjects_sessions_list: The Subjects-Sessions list file (in .tsv format).

    Returns:
        A clinica pipeline object containing the SpatialSVM pipeline.

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

        return ['dartel_input', 'input_image']

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """

        return ['regularized_image']

    def build_input_node(self):
        """Build and connect an input node to the pipeline.
        """

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from clinica.utils.stream import cprint
        from os.path import exists, join, abspath
        from os import listdir
        from clinica.utils.exceptions import ClinicaCAPSError, ClinicaException
        from clinica.utils.inputs import clinica_file_reader, clinica_group_reader
        import clinica.utils.input_files as input_files

        # Check that group-id already exists
        if not exists(join(abspath(self.caps_directory), 'groups', 'group-' + self.parameters['group_id'])):
            error_message = 'group_id : ' + self.parameters['group_id'] + ' does not exists, ' \
                            + 'please choose an other one. Groups that exist' \
                            + 's in your CAPS directory are : \n'
            list_groups = listdir(join(abspath(self.caps_directory), 'groups'))
            has_one_group = False
            for e in list_groups:
                if e.startswith('group-'):
                    error_message += e + ' \n'
                    has_one_group = True
            if not has_one_group:
                error_message = error_message + 'No group found ! ' \
                                + 'Use t1-volume pipeline if you do not ' \
                                + 'have a template yet ! '
            raise ValueError(error_message)

        read_parameters_node = npe.Node(name="LoadingCLIArguments",
                                        interface=nutil.IdentityInterface(fields=self.get_input_fields(),
                                                                          mandatory_inputs=True))
        image_type = self.parameters['image_type']
        pet_type = self.parameters['pet_type']
        no_pvc = self.parameters['no_pvc']

        all_errors = []
        if image_type == 't1':
            try:
                input_image = clinica_file_reader(self.subjects,
                                                  self.sessions,
                                                  self.caps_directory,
                                                  {'pattern': 't1/spm/dartel/group-' + self.parameters['group_id']
                                                              + '/*_T1w_segm-graymatter_space-Ixi549Space_modulated-on_probability.nii.gz',
                                                   'description': 'graymatter tissue segmented in T1w MRI in Ixi549 space',
                                                   'needed_pipeline': 't1-volume-tissue-segmentation'})
            except ClinicaException as e:
                all_errors.append(e)

        elif image_type is 'pet':
            if no_pvc.lower() == 'true':
                try:
                    input_image = clinica_file_reader(self.subjects,
                                                      self.sessions,
                                                      self.caps_directory,
                                                      {'pattern': 'pet/preprocessing/group-' + self.parameters['group_id']
                                                                  + '/*_pet_space-Ixi549Space_suvr-pons_pet.nii.gz',
                                                       'description': pet_type + ' PET in Ixi549 space',
                                                       'needed_pipeline': 'pet-volume'})
                except ClinicaException as e:
                    all_errors.append(e)

            elif no_pvc.lower() == 'false':
                try:
                    input_image = clinica_file_reader(self.subjects,
                                                      self.sessions,
                                                      self.caps_directory,
                                                      {'pattern': 'pet/preprocessing/group-' + self.parameters['group_id']
                                                                  + '_pet_space-Ixi549Space_pvc-rbv_suvr-pons_pet.nii.gz',
                                                       'description': pet_type + ' PET partial volume corrected (RBV) in Ixi549 space',
                                                       'needed_pipeline': 'pet-volume with PVC'})
                except ClinicaException as e:
                    all_errors.append(e)

            else:
                raise ValueError(no_pvc + ' is not a valid keyword for -no_pvc'
                                 + ':only True or False are accepted (string)')
        else:
            raise ValueError('Image type ' + image_type + ' unknown')
        try:
            dartel_input = clinica_group_reader(self.caps_directory,
                                                {'pattern': 'group-' + self.parameters['group_id'] + '_template.nii*',
                                                 'description': 'template file of group ' + self.parameters['group_id'],
                                                 'needed_pipeline': 't1-volume or t1-volume-create-dartel'})
        except ClinicaException as e:
            all_errors.append(e)

        # Raise all errors if some happened
        if len(all_errors) > 0:
            error_message = 'Clinica faced errors while trying to read files in your CAPS directories.\n'
            for msg in all_errors:
                error_message += str(msg)
            raise ClinicaCAPSError(error_message)

        read_parameters_node.inputs.dartel_input = dartel_input
        read_parameters_node.inputs.input_image = input_image

        self.connect([
            (read_parameters_node,      self.input_node,    [('dartel_input',    'dartel_input')]),
            (read_parameters_node,      self.input_node,    [('input_image',    'input_image')])

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

        import clinica.pipelines.machine_learning_spatial_svm.spatial_svm_utils as utils
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio

        fisher_tensor_generation = npe.Node(name="obtain_g_fisher_tensor",
                                            interface=nutil.Function(input_names=['dartel_input', 'FWHM'],
                                                                     output_names=['fisher_tensor', 'fisher_tensor_path'],
                                                                     function=utils.obtain_g_fisher_tensor))
        fisher_tensor_generation.inputs.FWHM = self.parameters['fwhm']

        time_step_generation = npe.Node(name='estimation_time_step',
                                        interface=nutil.Function(input_names=['dartel_input', 'FWHM', 'g'],
                                                                 output_names=['t_step', 'json_file'],
                                                                 function=utils.obtain_time_step_estimation))
        time_step_generation.inputs.FWHM = self.parameters['fwhm']

        heat_solver_equation = npe.MapNode(name='heat_solver_equation',
                                           interface=nutil.Function(input_names=['input_image', 'g',
                                                                                 'FWHM', 't_step', 'dartel_input'],
                                                                    output_names=['regularized_image'],
                                                                    function=utils.heat_solver_equation),
                                           iterfield=['input_image'])
        heat_solver_equation.inputs.FWHM = self.parameters['fwhm']

        datasink = npe.Node(nio.DataSink(),
                            name='sinker')
        datasink.inputs.base_directory = self.caps_directory
        datasink.inputs.parameterization = True
        if self.parameters['image_type'] == 't1':
            datasink.inputs.regexp_substitutions = [
                (r'(.*)/regularized_image/.*/(.*(sub-(.*)_ses-(.*))_T1w(.*)_probability(.*))$',
                 r'\1/subjects/sub-\4/ses-\5/machine_learning/input_spatial_svm/group-' + self.parameters[
                     'group_id'] + r'/\3_T1w\6_spatialregularization\7'),

                (r'(.*)json_file/(output_data.json)$',
                 r'\1/groups/group-' + self.parameters['group_id'] + r'/machine_learning/input_spatial_svm/group-' + self.parameters[
                     'group_id'] + r'_space-Ixi549Space_parameters.json'),

                (r'(.*)fisher_tensor_path/(output_fisher_tensor.npy)$',
                 r'\1/groups/group-' + self.parameters['group_id'] + r'/machine_learning/input_spatial_svm/group-' + self.parameters[
                     'group_id'] + r'_space-Ixi549Space_gram.npy')
            ]

        elif self.parameters['image_type'] == 'pet':
            datasink.inputs.regexp_substitutions = [
                (r'(.*)/regularized_image/.*/(.*(sub-(.*)_ses-(.*))_(task.*)_pet(.*))$',
                 r'\1/subjects/sub-\4/ses-\5/machine_learning/input_spatial_svm/group-' + self.parameters[
                     'group_id'] + r'/\3_\6_spatialregularization\7'),
                (r'(.*)json_file/(output_data.json)$',
                 r'\1/groups/group-' + self.parameters['group_id'] + r'/machine_learning/input_spatial_svm/group-' +
                 self.parameters['group_id'] + r'_space-Ixi549Space_parameters.json'),
                (r'(.*)fisher_tensor_path/(output_fisher_tensor.npy)$',
                 r'\1/groups/group-' + self.parameters['group_id'] + r'/machine_learning/input_spatial_svm/group-' +
                 self.parameters[
                     'group_id'] + r'_space-Ixi549Space_gram.npy')
            ]
        # Connection
        # ==========
        self.connect([
            (self.input_node,      fisher_tensor_generation,    [('dartel_input',    'dartel_input')]),
            (fisher_tensor_generation,      time_step_generation,    [('fisher_tensor',    'g')]),

            (self.input_node, time_step_generation, [('dartel_input', 'dartel_input')]),
            (self.input_node, heat_solver_equation, [('input_image', 'input_image')]),
            (fisher_tensor_generation, heat_solver_equation, [('fisher_tensor', 'g')]),
            (time_step_generation, heat_solver_equation, [('t_step', 't_step')]),
            (self.input_node, heat_solver_equation, [('dartel_input', 'dartel_input')]),

            (fisher_tensor_generation, datasink, [('fisher_tensor_path', 'fisher_tensor_path')]),
            (time_step_generation, datasink, [('json_file', 'json_file')]),
            (heat_solver_equation, datasink, [('regularized_image', 'regularized_image')])
        ])
