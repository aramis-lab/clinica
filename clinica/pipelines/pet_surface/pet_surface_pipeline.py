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
        from clinica.lib.pycaps.caps_layout import CAPSLayout

        def handle_files(subject, session, files, description, list_missing_files):
            from clinica.utils.stream import cprint
            missing_f = True
            for i in range(len(files)):
                if subject + '_' + session in files[i] and missing_f is True:
                    result = files[i]
                    missing_f = False
            if missing_f:
                list_missing_files.append('(' + subject + ' - ' + session + ')' + 'file for the following category was not found : ' + description)
                result = ' * missing *'
            cprint(description + ' = ' + result)
            return result

        read_parameters_node = npe.Node(name="LoadingCLIArguments",
                                        interface=nutil.IdentityInterface(
                                            fields=self.get_input_fields(),
                                            mandatory_inputs=True),
                                        synchronize=True)

        read_parameters_node.inputs.pet = []
        read_parameters_node.inputs.orig_nu = []
        read_parameters_node.inputs.psf = []
        read_parameters_node.inputs.white_surface_right = []
        read_parameters_node.inputs.white_surface_left = []
        read_parameters_node.inputs.destrieux_left = []
        read_parameters_node.inputs.destrieux_right = []
        read_parameters_node.inputs.desikan_left = []
        read_parameters_node.inputs.desikan_right = []
        caps_layout = CAPSLayout(self.caps_directory)

        cprint('------- INPUT FILES FOR EACH SUBJECTS -------')
        subjects_regex = '|'.join(sub[4:] for sub in self.subjects)
        unique_session = set(list(self.sessions))
        sessions_regex = '|'.join(sub[4:] for sub in unique_session)

        cprint('  * grabbing all pet files from BIDS folder')
        pet_file = self.bids_layout.get(return_type='file',
                                        type='pet',
                                        extensions=[self.parameters['pet_type'] + '_pet.nii',
                                                    self.parameters['pet_type'].upper() + '_pet.nii',
                                                    self.parameters['pet_type'] + '_pet.nii.gz',
                                                    self.parameters['pet_type'].upper() + '_pet.nii.gz'],
                                        subject=subjects_regex,
                                        session=sessions_regex)
        cprint('  * grabbing all orig_nu files from CAPS folder')
        orig_nu_file = caps_layout.get(return_type='file',
                                       freesurfer_file='orig_nu',
                                       extensions='mgz|MGZ',
                                       regex_search=True,
                                       subject=subjects_regex,
                                       session=sessions_regex)

        cprint('  * grabbing all point spread function files from BIDS folder')
        psf_file = self.bids_layout.get(return_type='file',
                                        type='pet',
                                        extensions=[self.parameters['pet_type'] + '_pet.json',
                                                    self.parameters['pet_type'].upper() + '_pet.json'],
                                        subject=subjects_regex,
                                        session=sessions_regex)

        cprint('  * grabbing all right hemisphere white surface files from CAPS folder')
        white_surface_left_file = caps_layout.get(return_type='file',
                                                  freesurfer_file='lh',
                                                  extensions='\\.white',
                                                  regex_search=True,
                                                  subject=subjects_regex,
                                                  session=sessions_regex)

        cprint('  * grabbing all left hemisphere white surface files from CAPS folder')
        white_surface_right_file = caps_layout.get(return_type='file',
                                                   freesurfer_file='rh',
                                                   extensions='\\.white',
                                                   regex_search=True,
                                                   subject=subjects_regex,
                                                   session=sessions_regex)

        cprint('  * grabbing all left hemisphere Destrieux atlases in subject\'s space from CAPS folder')
        destrieux_left_files = caps_layout.get(return_type='file',
                                               freesurfer_file='lh.aparc.a2009s.annot',
                                               regex_search=True,
                                               subject=subjects_regex,
                                               session=sessions_regex)

        cprint('  * grabbing all right hemisphere Destrieux atlases in subject\'s space from CAPS folder')
        destrieux_right_files = caps_layout.get(return_type='file',
                                                freesurfer_file='rh.aparc.a2009s.annot',
                                                regex_search=True,
                                                subject=subjects_regex,
                                                session=sessions_regex)

        cprint('  * grabbing all left hemisphere Desikan-Killiany atlases in subject\'s space from CAPS folder')
        desikan_left_files = caps_layout.get(return_type='file',
                                             freesurfer_file='lh.aparc.annot',
                                             regex_search=True,
                                             subject=subjects_regex,
                                             session=sessions_regex)

        cprint('  * grabbing all right hemisphere Desikan-Killiany atlases in subject\'s space from CAPS folder')
        desikan_right_files = caps_layout.get(return_type='file',
                                              freesurfer_file='rh.aparc.annot',
                                              regex_search=True,
                                              subject=subjects_regex,
                                              session=sessions_regex)

        missing_files = []
        for i in range(len(self.subjects)):
            read_parameters_node.inputs.pet.append(handle_files(self.subjects[i],
                                                                self.sessions[i],
                                                                pet_file,
                                                                'PET file',
                                                                missing_files))
            read_parameters_node.inputs.orig_nu.append(handle_files(self.subjects[i],
                                                                    self.sessions[i],
                                                                    orig_nu_file,
                                                                    'orig_nu',
                                                                    missing_files))
            read_parameters_node.inputs.psf.append(handle_files(self.subjects[i],
                                                                self.sessions[i],
                                                                psf_file,
                                                                'pet json',
                                                                missing_files))
            read_parameters_node.inputs.white_surface_left.append(handle_files(self.subjects[i],
                                                                               self.sessions[i],
                                                                               white_surface_left_file,
                                                                               'lh.white',
                                                                               missing_files))
            read_parameters_node.inputs.white_surface_right.append(handle_files(self.subjects[i],
                                                                                self.sessions[i],
                                                                                white_surface_right_file,
                                                                                'rh.white',
                                                                                missing_files))
            read_parameters_node.inputs.destrieux_left.append(handle_files(self.subjects[i],
                                                                           self.sessions[i],
                                                                           destrieux_left_files,
                                                                           'destrieux lh',
                                                                           missing_files))
            read_parameters_node.inputs.destrieux_right.append(handle_files(self.subjects[i],
                                                                            self.sessions[i],
                                                                            destrieux_right_files,
                                                                            'destrieux rh',
                                                                            missing_files))
            read_parameters_node.inputs.desikan_left.append(handle_files(self.subjects[i],
                                                                         self.sessions[i],
                                                                         desikan_left_files,
                                                                         'desikan lh',
                                                                         missing_files))
            read_parameters_node.inputs.desikan_right.append(handle_files(self.subjects[i],
                                                                          self.sessions[i],
                                                                          desikan_right_files,
                                                                          'desikan rh',
                                                                          missing_files))

        if len(missing_files) > 0:
            error_string = 'The following elements are missing : \n'
            for missing_f in missing_files:
                error_string += missing_f + ' \n'
            raise Exception(error_string)

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
                                                          'destrieux_right'],
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
