# coding: utf-8

import clinica.pipelines.engine as cpe


class PetSurface(cpe.Pipeline):
    """PetSurface - Surface-based processing of PET images.

    Returns:
        A clinica pipeline object containing the PetSurface pipeline.
    """

    def check_pipeline_parameters(self):
        """Check pipeline parameters."""
        if 'pet_tracer' not in self.parameters.keys():
            self.parameters['pet_tracer'] = 'fdg'

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file."""
        pass

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Returns:
            A list of (string) input fields name.
        """
        return ['orig_nu',
                'pet',
                'psf',
                'white_surface_left',
                'white_surface_right',
                'destrieux_left',
                'destrieux_right',
                'desikan_left',
                'desikan_right']

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline."""
        return []

    def build_input_node(self):
        """Build and connect an input node to the pipeline."""
        import os
        from clinica.utils.ux import print_images_to_process
        from clinica.utils.filemanip import save_participants_sessions
        from clinica.utils.stream import cprint

        if self.parameters['longitudinal']:
            self.build_input_node_longitudinal()
        else:
            self.build_input_node_cross_sectional()

        # Save subjects to process in <WD>/<Pipeline.name>/participants.tsv
        folder_participants_tsv = os.path.join(self.base_dir, self.name)
        save_participants_sessions(self.subjects, self.sessions, folder_participants_tsv)

        if len(self.subjects):
            print_images_to_process(self.subjects, self.sessions)
            cprint('List available in %s' % os.path.join(folder_participants_tsv, 'participants.tsv'))
            cprint('The pipeline will last approximately a few hours per image.')

    def build_input_node_longitudinal(self):
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from clinica.utils.inputs import clinica_file_reader
        from clinica.utils.exceptions import ClinicaException
        from clinica.iotools.utils.data_handling import check_relative_volume_location_in_world_coordinate_system
        import clinica.utils.input_files as input_files

        read_parameters_node = npe.Node(name="LoadingCLIArguments",
                                        interface=nutil.IdentityInterface(
                                            fields=self.get_input_fields(),
                                            mandatory_inputs=True),
                                        synchronize=True)

        if self.parameters['pet_tracer'].lower() == 'fdg':
            pet_file_to_grab = input_files.PET_FDG_NII
            pet_json_file_to_grab = input_files.PET_FDG_JSON
        elif self.parameters['pet_tracer'].lower() == 'av45':
            pet_file_to_grab = input_files.PET_AV45_NII
            pet_json_file_to_grab = input_files.PET_AV45_JSON
        else:
            raise NotImplementedError('Only "fdg" or "av45" tracers are currently accepted (given tracer "%s").' %
                                      self.parameters['pet_tracer'])
        all_errors = []
        try:

            read_parameters_node.inputs.pet = clinica_file_reader(self.subjects,
                                                                  self.sessions,
                                                                  self.bids_directory,
                                                                  pet_file_to_grab)
        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.orig_nu = clinica_file_reader(self.subjects,
                                                                      self.sessions,
                                                                      self.caps_directory,
                                                                      input_files.T1_FS_LONG_ORIG_NU)

        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.psf = clinica_file_reader(self.subjects,
                                                                  self.sessions,
                                                                  self.bids_directory,
                                                                  pet_json_file_to_grab)
        except ClinicaException as e:
            all_errors.append(e)

        try:

            read_parameters_node.inputs.white_surface_right = clinica_file_reader(self.subjects,
                                                                                  self.sessions,
                                                                                  self.caps_directory,
                                                                                  input_files.T1_FS_LONG_SURF_R)
        except ClinicaException as e:
            all_errors.append(e)

        try:

            read_parameters_node.inputs.white_surface_left = clinica_file_reader(self.subjects,
                                                                                 self.sessions,
                                                                                 self.caps_directory,
                                                                                 input_files.T1_FS_LONG_SURF_L)

        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.destrieux_left = clinica_file_reader(self.subjects,
                                                                             self.sessions,
                                                                             self.caps_directory,
                                                                             input_files.T1_FS_LONG_DESTRIEUX_PARC_L)

        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.destrieux_right = clinica_file_reader(self.subjects,
                                                                              self.sessions,
                                                                              self.caps_directory,
                                                                              input_files.T1_FS_LONG_DESTRIEUX_PARC_R)

        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.desikan_left = clinica_file_reader(self.subjects,
                                                                           self.sessions,
                                                                           self.caps_directory,
                                                                           input_files.T1_FS_LONG_DESIKAN_PARC_L)

        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.desikan_right = clinica_file_reader(self.subjects,
                                                                            self.sessions,
                                                                            self.caps_directory,
                                                                            input_files.T1_FS_LONG_DESIKAN_PARC_R)

        except ClinicaException as e:
            all_errors.append(e)

        if len(all_errors) > 0:
            error_message = 'Clinica faced errors while trying to read files in your BIDS or CAPS directories.\n'
            for msg in all_errors:
                error_message += str(msg)
            raise ClinicaException(error_message)

        check_relative_volume_location_in_world_coordinate_system('T1w-MRI (orig_nu.mgz)',
                                                                  read_parameters_node.inputs.orig_nu,
                                                                  self.parameters['pet_tracer'].upper() + ' PET',
                                                                  read_parameters_node.inputs.pet,
                                                                  self.bids_directory,
                                                                  self.parameters['pet_tracer'].lower())

        self.connect([
            (read_parameters_node, self.input_node, [('pet', 'pet')]),
            (read_parameters_node, self.input_node, [('orig_nu', 'orig_nu')]),
            (read_parameters_node, self.input_node, [('psf', 'psf')]),
            (read_parameters_node, self.input_node, [('white_surface_left', 'white_surface_left')]),
            (read_parameters_node, self.input_node, [('white_surface_right', 'white_surface_right')]),
            (read_parameters_node, self.input_node, [('destrieux_left', 'destrieux_left')]),
            (read_parameters_node, self.input_node, [('destrieux_right', 'destrieux_right')]),
            (read_parameters_node, self.input_node, [('desikan_left', 'desikan_left')]),
            (read_parameters_node, self.input_node, [('desikan_right', 'desikan_right')])
        ])

    def build_input_node_cross_sectional(self):
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from clinica.utils.inputs import clinica_file_reader
        from clinica.utils.exceptions import ClinicaException
        from clinica.iotools.utils.data_handling import check_relative_volume_location_in_world_coordinate_system
        import clinica.utils.input_files as input_files

        read_parameters_node = npe.Node(name="LoadingCLIArguments",
                                        interface=nutil.IdentityInterface(
                                            fields=self.get_input_fields(),
                                            mandatory_inputs=True),
                                        synchronize=True)

        if self.parameters['pet_tracer'].lower() == 'fdg':
            pet_file_to_grab = input_files.PET_FDG_NII
            pet_json_file_to_grab = input_files.PET_FDG_JSON
        elif self.parameters['pet_tracer'].lower() == 'av45':
            pet_file_to_grab = input_files.PET_AV45_NII
            pet_json_file_to_grab = input_files.PET_AV45_JSON
        else:
            raise NotImplementedError('Only "fdg" or "av45" tracers are currently accepted (given tracer "%s").' %
                                      self.parameters['pet_tracer'])
        all_errors = []
        try:

            read_parameters_node.inputs.pet = clinica_file_reader(self.subjects,
                                                                  self.sessions,
                                                                  self.bids_directory,
                                                                  pet_file_to_grab)
        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.orig_nu = clinica_file_reader(self.subjects,
                                                                      self.sessions,
                                                                      self.caps_directory,
                                                                      input_files.T1_FS_ORIG_NU)
        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.psf = clinica_file_reader(self.subjects,
                                                                  self.sessions,
                                                                  self.bids_directory,
                                                                  pet_json_file_to_grab)
        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.white_surface_right = clinica_file_reader(self.subjects,
                                                                                  self.sessions,
                                                                                  self.caps_directory,
                                                                                  input_files.T1_FS_WM_SURF_R)
        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.white_surface_left = clinica_file_reader(self.subjects,
                                                                                 self.sessions,
                                                                                 self.caps_directory,
                                                                                 input_files.T1_FS_WM_SURF_L)
        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.destrieux_left = clinica_file_reader(self.subjects,
                                                                             self.sessions,
                                                                             self.caps_directory,
                                                                             input_files.T1_FS_DESTRIEUX_PARC_L)
        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.destrieux_right = clinica_file_reader(self.subjects,
                                                                              self.sessions,
                                                                              self.caps_directory,
                                                                              input_files.T1_FS_DESTRIEUX_PARC_R)
        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.desikan_left = clinica_file_reader(self.subjects,
                                                                           self.sessions,
                                                                           self.caps_directory,
                                                                           input_files.T1_FS_DESIKAN_PARC_L)
        except ClinicaException as e:
            all_errors.append(e)

        try:
            read_parameters_node.inputs.desikan_right = clinica_file_reader(self.subjects,
                                                                            self.sessions,
                                                                            self.caps_directory,
                                                                            input_files.T1_FS_DESIKAN_PARC_R)
        except ClinicaException as e:
            all_errors.append(e)

        if len(all_errors) > 0:
            error_message = 'Clinica faced errors while trying to read files in your BIDS or CAPS directories.\n'
            for msg in all_errors:
                error_message += str(msg)
            raise ClinicaException(error_message)

        check_relative_volume_location_in_world_coordinate_system('T1w-MRI (orig_nu.mgz)', read_parameters_node.inputs.orig_nu,
                                                                  self.parameters['pet_tracer'].upper() + ' PET', read_parameters_node.inputs.pet,
                                                                  self.bids_directory,
                                                                  self.parameters['pet_tracer'].lower())

        self.connect([
            (read_parameters_node,      self.input_node,    [('pet',                    'pet')]),
            (read_parameters_node,      self.input_node,    [('orig_nu',                'orig_nu')]),
            (read_parameters_node,      self.input_node,    [('psf',                    'psf')]),
            (read_parameters_node,      self.input_node,    [('white_surface_left',     'white_surface_left')]),
            (read_parameters_node,      self.input_node,    [('white_surface_right',    'white_surface_right')]),
            (read_parameters_node,      self.input_node,    [('destrieux_left',         'destrieux_left')]),
            (read_parameters_node,      self.input_node,    [('destrieux_right',        'destrieux_right')]),
            (read_parameters_node,      self.input_node,    [('desikan_left',           'desikan_left')]),
            (read_parameters_node,      self.input_node,    [('desikan_right',          'desikan_right')])
        ])

    def build_output_node(self):
        """Build and connect an output node to the pipeline."""
        pass

    def build_core_nodes(self):
        """The function get_wf constructs a pipeline for one subject (in pet_surface_utils.py) and runs it.
        We use iterables to give to the node all the files and information needed.
        """
        # TODO(@arnaud.marcoux): Convert it to a Node with iterables + MapNodes.
        #   I'm experimenting something to avoid the "MapNode of MapNode" case
        #   with iterables. I'll try to apply it on the tractography pipeline.
        #   Check it out to get inspiration from it when it's ready.

        import os
        import nipype.pipeline.engine as npe
        import nipype.interfaces.utility as niu
        import clinica.pipelines.pet_surface.pet_surface_utils as utils
        from colorama import Fore
        from nipype.interfaces import spm
        import platform

        full_pipe = npe.MapNode(niu.Function(input_names=['subject_id',
                                                          'session_id',
                                                          'caps_dir',
                                                          'psf',
                                                          'pet',
                                                          'orig_nu',
                                                          'white_surface_left',
                                                          'white_surface_right',
                                                          'working_directory_subjects',
                                                          'pet_tracer',
                                                          'csv_segmentation',
                                                          'subcortical_eroded_mask',
                                                          'matscript_folder_inverse_deformation',
                                                          'desikan_left',
                                                          'desikan_right',
                                                          'destrieux_left',
                                                          'destrieux_right',
                                                          'use_spm_standalone',
                                                          'is_longitudinal'],
                                             output_names=[],
                                             function=utils.get_wf),
                                name='full_pipeline_mapnode',
                                iterfield=['subject_id',
                                           'session_id',
                                           'psf',
                                           'pet',
                                           'orig_nu',
                                           'white_surface_left',
                                           'white_surface_right',
                                           'desikan_left',
                                           'desikan_right',
                                           'destrieux_left',
                                           'destrieux_right'])

        full_pipe.inputs.subject_id = self.subjects
        full_pipe.inputs.caps_dir = self.caps_directory
        full_pipe.inputs.session_id = self.sessions
        full_pipe.inputs.working_directory_subjects = self.base_dir
        full_pipe.inputs.pet_tracer = self.parameters['pet_tracer']
        full_pipe.inputs.csv_segmentation = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                                                         '..',
                                                                         '..',
                                                                         'resources',
                                                                         'label_conversion_gtmsegmentation.csv'))
        if self.parameters['pet_tracer'].lower() == 'fdg':
            full_pipe.inputs.subcortical_eroded_mask = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                                                                    '..',
                                                                                    '..',
                                                                                    'resources',
                                                                                    'masks',
                                                                                    'region-pons_eroded-6mm_mask.nii.gz'))
        elif self.parameters['pet_tracer'].lower() == 'av45':
            full_pipe.inputs.subcortical_eroded_mask = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                                                                    '..',
                                                                                    '..',
                                                                                    'resources',
                                                                                    'masks',
                                                                                    'region-cerebellumPons_eroded-6mm_mask.nii.gz'))

        full_pipe.inputs.matscript_folder_inverse_deformation = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))
        full_pipe.inputs.is_longitudinal = self.parameters['longitudinal']

        # This section of code determines whether to use SPM standalone or not
        full_pipe.inputs.use_spm_standalone = False
        if all(elem in os.environ.keys() for elem in ['SPMSTANDALONE_HOME', 'MCR_HOME']):
            if os.path.exists(os.path.expandvars('$SPMSTANDALONE_HOME')) and os.path.exists(os.path.expandvars('$MCR_HOME')):
                print(Fore.GREEN + 'SPM standalone has been found and will be used in this pipeline' + Fore.RESET)
                if platform.system() == 'Darwin':
                    matlab_cmd = ('cd ' + os.path.expandvars('$SPMSTANDALONE_HOME') + ' && ' + './run_spm12.sh'
                                  + ' ' + os.environ['MCR_HOME']
                                  + ' script')
                elif platform.system() == 'Linux':
                    matlab_cmd = (os.path.join(os.path.expandvars('$SPMSTANDALONE_HOME'), 'run_spm12.sh')
                                  + ' ' + os.environ['MCR_HOME']
                                  + ' script')
                spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)
                full_pipe.inputs.use_spm_standalone = True
            else:
                raise FileNotFoundError('$SPMSTANDALONE_HOME and $MCR_HOME are defined, but linked to non existent folder ')

        # Connection
        # ==========
        self.connect([
            (self.input_node, full_pipe, [('pet', 'pet')]),
            (self.input_node, full_pipe, [('psf', 'psf')]),
            (self.input_node, full_pipe, [('white_surface_left', 'white_surface_left')]),
            (self.input_node, full_pipe, [('white_surface_right', 'white_surface_right')]),
            (self.input_node, full_pipe, [('orig_nu', 'orig_nu')]),
            (self.input_node, full_pipe, [('destrieux_left', 'destrieux_left')]),
            (self.input_node, full_pipe, [('destrieux_right', 'destrieux_right')]),
            (self.input_node, full_pipe, [('desikan_left', 'desikan_left')]),
            (self.input_node, full_pipe, [('desikan_right', 'desikan_right')])
        ])
