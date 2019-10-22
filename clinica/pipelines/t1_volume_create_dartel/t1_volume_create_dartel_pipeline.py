# coding: utf8

import clinica.pipelines.engine as cpe

__author__ = "Jorge Samper-Gonzalez"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Jorge Samper-Gonzalez"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper-Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


class T1VolumeCreateDartel(cpe.Pipeline):
    """T1VolumeCreateDartel - Create new Dartel template.

    Args:
        input_dir: A BIDS directory.
        output_dir: An empty output directory where CAPS structured data will be written.
        subjects_sessions_list: The Subjects-Sessions list file (in .tsv format).

    Returns:
        A clinica pipeline object containing the T1VolumeCreateDartel pipeline.
    """
    def __init__(self, bids_directory=None, caps_directory=None, tsv_file=None, name=None, group_id='default'):
        import os

        super(T1VolumeCreateDartel, self).__init__(bids_directory, caps_directory, tsv_file, name)

        if not group_id.isalnum():
            raise ValueError('Not valid group_id value. It must be composed only by letters and/or numbers')
        self._group_id = group_id

        # Check that group does not already exists
        if os.path.exists(os.path.join(os.path.abspath(caps_directory), 'groups', 'group-' + group_id)):
            error_message = 'group_id : ' + group_id + ' already exists, please choose another one.' \
                            + ' Groups that exist in your CAPS directory are: \n'
            list_groups = os.listdir(os.path.join(os.path.abspath(caps_directory), 'groups'))
            for e in list_groups:
                if e.startswith('group-'):
                    error_message += e + ' \n'
            raise ValueError(error_message)

        # Check that there is at least 2 subjects
        if len(self.subjects) <= 1:
            raise ValueError('This pipeline needs at least 2 subjects to perform DARTEL, and found '
                             + str(len(self.subjects)) + ' only in ' + self.tsv_file + '.')

        # Default parameters
        self._parameters = {'dartel_tissues': [1, 2, 3],
                            'iteration_parameters': None,
                            'optimization_parameters': None,
                            'regularization_form': None,
                            'template_prefix': None
                            }

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file.
        """
        pass

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipelines.

        Returns:
            A list of (string) input fields name.
        """

        return ['dartel_inputs']

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipelines.

        Returns:
            A list of (string) output fields name.
        """

        return ['final_template_file', 'template_files', 'dartel_flow_fields']

    def build_input_node(self):
        """Build and connect an input node to the pipelines.
        """

        import nipype.pipeline.engine as npe
        from clinica.utils.inputs import clinica_file_reader
        from clinica.utils.exceptions import ClinicaBIDSError, ClinicaException
        import nipype.interfaces.utility as nutil

        tissue_names = {1: 'graymatter',
                        2: 'whitematter',
                        3: 'csf',
                        4: 'bone',
                        5: 'softtissue',
                        6: 'background'
                        }

        read_parameters_node = npe.Node(name="LoadingCLIArguments",
                                        interface=nutil.IdentityInterface(fields=self.get_input_fields(),
                                                                          mandatory_inputs=True))
        all_errors = []
        d_input = []
        for tissue_number in self.parameters['dartel_tissues']:
            try:
                current_file = clinica_file_reader(self.subjects,
                                                   self.sessions,
                                                   self.caps_directory,
                                                   {'pattern': 't1/spm/segmentation/dartel_input/*_*_T1w_segm-'
                                                               + tissue_names[tissue_number] + '_dartelinput.nii*',
                                                    'description': 'Dartel input for tissue ' + tissue_names[tissue_number]
                                                                   + ' from T1w MRI',
                                                    'needed_pipeline': 't1-volume-tissue-segmentation'})
                d_input.append(current_file)
            except ClinicaException as e:
                all_errors.append(e)

        # Raise all errors if some happened
        if len(all_errors) > 0:
            error_message = 'Clinica faced errors while trying to read files in your BIDS or CAPS directories.\n'
            for msg in all_errors:
                error_message += str(msg)
            raise RuntimeError(error_message)

        # d_input is a list of size len(self.parameters['dartel_tissues'])
        #     Each element of this list is a list of size len(self.subjects)
        read_parameters_node.inputs.dartel_inputs = d_input

        self.connect([
            (read_parameters_node, self.input_node, [('dartel_inputs', 'dartel_input_images')])
        ])

    def build_output_node(self):
        """Build and connect an output node to the pipelines.
        """
        import os.path as op
        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio
        import re
        from clinica.utils.io import zip_nii

        # Writing flowfields into CAPS
        # ============================
        write_flowfields_node = npe.MapNode(name='write_flowfields_node',
                                            iterfield=['container', 'flow_fields'],
                                            interface=nio.DataSink(infields=['flow_fields']))
        write_flowfields_node.inputs.base_directory = self.caps_directory
        write_flowfields_node.inputs.parameterization = False
        write_flowfields_node.inputs.container = ['subjects/' + self.subjects[i] + '/' + self.sessions[i] +
                                                  '/t1/spm/dartel/group-' + self._group_id
                                                  for i in range(len(self.subjects))]
        write_flowfields_node.inputs.regexp_substitutions = [
            (r'(.*)_Template(\.nii(\.gz)?)$', r'\1\2'),
            (r'(.*)c1(sub-.*)(\.nii(\.gz)?)$', r'\1\2_segm-graymatter\3'),
            (r'(.*)c2(sub-.*)(\.nii(\.gz)?)$', r'\1\2_segm-whitematter\3'),
            (r'(.*)c3(sub-.*)(\.nii(\.gz)?)$', r'\1\2_segm-csf\3'),
            (r'(.*)c4(sub-.*)(\.nii(\.gz)?)$', r'\1\2_segm-bone\3'),
            (r'(.*)c5(sub-.*)(\.nii(\.gz)?)$', r'\1\2_segm-softtissue\3'),
            (r'(.*)c6(sub-.*)(\.nii(\.gz)?)$', r'\1\2_segm-background\3'),
            (r'(.*)r(sub-.*)(\.nii(\.gz)?)$', r'\1\2\3'),
            (r'(.*)_dartelinput(\.nii(\.gz)?)$', r'\1\2'),
            (r'(.*)flow_fields/u_(sub-.*)_segm-.*(\.nii(\.gz)?)$',
             r'\1\2_target-' + re.escape(self._group_id) + r'_transformation-forward_deformation\3'),
            (r'trait_added', r'')
        ]

        # Writing templates into CAPS
        # ===========================
        write_template_node = npe.Node(nio.DataSink(), name='write_template_node')
        write_template_node.inputs.parameterization = False
        write_template_node.inputs.base_directory = self.caps_directory
        write_template_node.inputs.container = op.join('groups/group-' + self._group_id, 't1')
        write_template_node.inputs.regexp_substitutions = [
            (r'(.*)final_template_file/.*(\.nii(\.gz)?)$',
             r'\1group-' + re.escape(self._group_id) + r'_template\2'),
            (r'(.*)template_files/.*([0-9])(\.nii(\.gz)?)$',
             r'\1group-' + re.escape(self._group_id) + r'_iteration-\2_template\3')
        ]

        self.connect([
            (self.output_node, write_flowfields_node, [(('dartel_flow_fields', zip_nii, True), 'flow_fields')]),
            (self.output_node, write_template_node, [(('final_template_file', zip_nii, True), 'final_template_file'),
                                                     (('template_files', zip_nii, True), 'template_files')])
        ])

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipelines.
        """
        import nipype.interfaces.spm as spm
        import nipype.pipeline.engine as npe
        import nipype.interfaces.utility as nutil
        from clinica.utils.io import unzip_nii
        from clinica.utils.spm import get_tpm

        # Get Tissue Probability Map from SPM
        tissue_map = get_tpm()

        # Unzipping
        # =========
        unzip_node = npe.MapNode(nutil.Function(input_names=['in_file'],
                                                output_names=['out_file'],
                                                function=unzip_nii),
                                 name='unzip_node', iterfield=['in_file'])

        # DARTEL template
        # ===============
        dartel_template = npe.Node(spm.DARTEL(),
                                   name='dartel_template')

        if self.parameters['iteration_parameters'] is not None:
            dartel_template.inputs.iteration_parameters = self.parameters['iteration_parameters']
        if self.parameters['optimization_parameters'] is not None:
            dartel_template.inputs.optimization_parameters = self.parameters['optimization_parameters']
        if self.parameters['regularization_form'] is not None:
            dartel_template.inputs.regularization_form = self.parameters['regularization_form']
        if self.parameters['template_prefix'] is not None:
            dartel_template.inputs.template_prefix = self.parameters['template_prefix']

        # Connection
        # ==========
        self.connect([
            (self.input_node, unzip_node,    [('dartel_input_images', 'in_file')]),
            (unzip_node, dartel_template,    [('out_file', 'image_files')]),
            (dartel_template, self.output_node, [('dartel_flow_fields', 'dartel_flow_fields'),
                                                 ('final_template_file', 'final_template_file'),
                                                 ('template_files', 'template_files')])
        ])
