# coding: utf-8

__author__ = "Arnaud Marcoux"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Arnaud Marcoux", "Michael Bacci"]
__license__ = "See LICENSE.txt file"
__version__ = "1.0.0"
__maintainer__ = "Arnaud Marcoux"
__email__ = "arnaud.marcoux@inria.fr"
__status__ = "Development"

import clinica.pipelines.engine as cpe


class PetSurface(cpe.Pipeline):
    """Project PET signal onto the surface of the cortex.

    Args:
        input_dir: A BIDS directory.
        output_dir: An empty output directory where CAPS structured data will be
            written.
        subjects_sessions_list: The Subjects-Sessions list file (in .tsv
            format).

    Returns:
        A clinica pipeline object containing the PetSurface pipeline.

    """

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file.
        """
        pass

    def get_input_fields(self):
        """Possible inputs of this pipeline.
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
        """Specify the list of possible outputs of this pipeline.
        """

        return []

    def build_input_node(self):
        """We iterate over subjects to get all the files needed to run the pipeline
        """
        from clinica.utils.stream import cprint
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from clinica.utils.inputs import clinica_file_reader

        read_parameters_node = npe.Node(name="LoadingCLIArguments",
                                        interface=nutil.IdentityInterface(
                                            fields=self.get_input_fields(),
                                            mandatory_inputs=True),
                                        synchronize=True)
        all_errors = []
        warning_freesurfer = False

        read_parameters_node.inputs.pet, error_str = clinica_file_reader(self.subjects,
                                                                         self.sessions,
                                                                         self.bids_directory,
                                                                         '*' + self.parameters['pet_type'] + '_pet.nii*')
        if error_str:
            all_errors.append(error_str)

        read_parameters_node.inputs.orig_nu, error_str = clinica_file_reader(self.subjects,
                                                                             self.sessions,
                                                                             self.caps_directory,
                                                                             'orig_nu.mgz')
        if error_str:
            warning_freesurfer = True
            all_errors.append(error_str)
        read_parameters_node.inputs.psf, error_str = clinica_file_reader(self.subjects,
                                                                         self.sessions,
                                                                         self.bids_directory,
                                                                         '*' + self.parameters['pet_type'] + '_pet.json')
        if error_str:
            warning_freesurfer = True
            all_errors.append(error_str)
        read_parameters_node.inputs.white_surface_right, error_str = clinica_file_reader(self.subjects,
                                                                                         self.sessions,
                                                                                         self.caps_directory,
                                                                                         'sub-*_ses-*/surf/rh.white')
        if error_str:
            warning_freesurfer = True
            all_errors.append(error_str)
        read_parameters_node.inputs.white_surface_left, error_str = clinica_file_reader(self.subjects,
                                                                                        self.sessions,
                                                                                        self.caps_directory,
                                                                                        'sub-*_ses-*/surf/lh.white')
        if error_str:
            warning_freesurfer = True
            all_errors.append(error_str)
        read_parameters_node.inputs.destrieux_left, error_str = clinica_file_reader(self.subjects,
                                                                                    self.sessions,
                                                                                    self.caps_directory,
                                                                                    'sub-*_ses-*/*/lh.aparc.a2009s.annot')
        if error_str:
            warning_freesurfer = True
            all_errors.append(error_str)
        read_parameters_node.inputs.destrieux_right, error_str = clinica_file_reader(self.subjects,
                                                                                     self.sessions,
                                                                                     self.caps_directory,
                                                                                     'sub-*_ses-*/*/rh.aparc.a2009s.annot')
        if error_str:
            warning_freesurfer = True
            all_errors.append(error_str)
        read_parameters_node.inputs.desikan_left, error_str = clinica_file_reader(self.subjects,
                                                                                  self.sessions,
                                                                                  self.caps_directory,
                                                                                  'sub-*_ses-*/*/lh.aparc.annot')
        if error_str:
            warning_freesurfer = True
            all_errors.append(error_str)
        read_parameters_node.inputs.desikan_right, error_str = clinica_file_reader(self.subjects,
                                                                                   self.sessions,
                                                                                   self.caps_directory,
                                                                                   'sub-*_ses-*/*/rh.aparc.annot')
        if error_str:
            all_errors.append(error_str)

        if len(all_errors) > 0:
            error_message = 'Clinica faced errors while trying to read files in your BIDS or CAPS directories.\n'
            if warning_freesurfer:
                error_message += 'Please note that the t1-freesurfer pipeline needs to have run on all the subjects.\n'
            for msg in all_errors:
                error_message += msg
            raise RuntimeError(error_message)

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
        """Build and connect an output node to the pipeline.
        """

        # In the same idea as the input node, this output node is supposedly
        # used to write the output fields in a CAPS. It should be executed only
        # if this pipeline output is not already connected to a next Clinica
        # pipeline.

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
                                                          'pet_type',
                                                          'csv_segmentation',
                                                          'subcortical_eroded_mask',
                                                          'matscript_folder_inverse_deformation',
                                                          'desikan_left',
                                                          'desikan_right',
                                                          'destrieux_left',
                                                          'destrieux_right',
                                                          'use_spm_standalone'],
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
        full_pipe.inputs.working_directory_subjects = self.parameters['wd']
        full_pipe.inputs.pet_type = self.parameters['pet_type']
        full_pipe.inputs.csv_segmentation = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                                                         '..',
                                                                         '..',
                                                                         'resources',
                                                                         'label_conversion_gtmsegmentation.csv'))
        if self.parameters['pet_type'].lower() == 'fdg':
            full_pipe.inputs.subcortical_eroded_mask = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                                                                    '..',
                                                                                    '..',
                                                                                    'resources',
                                                                                    'masks',
                                                                                    'region-pons_eroded-6mm_mask.nii.gz'))
        elif self.parameters['pet_type'].lower() == 'av45':
            full_pipe.inputs.subcortical_eroded_mask = os.path.abspath(os.path.join(os.path.dirname(os.path.realpath(__file__)),
                                                                                    '..',
                                                                                    '..',
                                                                                    'resources',
                                                                                    'masks',
                                                                                    'region-cerebellumPons_eroded-6mm_mask.nii.gz'))

        full_pipe.inputs.matscript_folder_inverse_deformation = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

        # This section of code determines wether to use SPM standalone or not
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
