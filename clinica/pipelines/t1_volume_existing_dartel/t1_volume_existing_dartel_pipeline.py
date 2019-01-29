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


class T1VolumeExistingDartel(cpe.Pipeline):
    """T1VolumeExistingDartel - Reuse existing Dartel template.

    Args:
        input_dir: A BIDS directory.
        output_dir: An empty output directory where CAPS structured data will
        be written.  subjects_sessions_list: The Subjects-Sessions list file
        (in .tsv format).

    Returns:
        A clinica pipeline object containing the T1VolumeExistingDartel pipeline.

    Raises:


    """

    def __init__(self, bids_directory=None, caps_directory=None, tsv_file=None, name=None, group_id='default'):
        super(T1VolumeExistingDartel, self).__init__(bids_directory, caps_directory, tsv_file, name)

        if not group_id.isalnum():
            raise ValueError('Not valid group_id value. It must be composed only by letters and/or numbers')

        self._group_id = group_id

        # Default parameters
        self._parameters = {'tissues': [1, 2, 3],
                            'iteration_parameters': None,
                            'optimization_parameters': None,
                            'regularization_form': None
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

        return ['dartel_input_images', 'dartel_iteration_templates']

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipelines.

        Returns:
            A list of (string) output fields name.
        """

        return ['dartel_flow_fields']

    def build_input_node(self):
        """Build and connect an input node to the pipelines.
        """

        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio

        tissue_names = {1: 'graymatter',
                        2: 'whitematter',
                        3: 'csf',
                        4: 'bone',
                        5: 'softtissue',
                        6: 'background'
                        }

        # Dartel Input Tissues DataGrabber
        # =================================
        dartel_input_reader = npe.MapNode(nio.DataGrabber(infields=['subject_id', 'session',
                                                                    'subject_repeat', 'session_repeat',
                                                                    'tissue'],
                                                          outfields=['out_files']),
                                          name="dartel_input_reader",
                                          iterfield=['tissue'])

        dartel_input_reader.inputs.base_directory = self.caps_directory
        dartel_input_reader.inputs.template = 'subjects/%s/%s/t1/spm/segmentation/dartel_input/%s_%s_T1w_segm-%s_dartelinput.nii*'
        dartel_input_reader.inputs.subject_id = self.subjects
        dartel_input_reader.inputs.session = self.sessions
        dartel_input_reader.inputs.tissue = [tissue_names[t] for t in self.parameters['tissues']]
        dartel_input_reader.inputs.subject_repeat = self.subjects
        dartel_input_reader.inputs.session_repeat = self.sessions
        dartel_input_reader.inputs.sort_filelist = False

        # Dartel Templates DataGrabber
        # ============================
        templates_reader = npe.MapNode(nio.DataGrabber(infields=['iteration'],
                                                       outfields=['out_files']),
                                       name="templates_reader",
                                       iterfield=['iteration'])
        templates_reader.inputs.base_directory = self.caps_directory
        templates_reader.inputs.template = 'groups/group-' + self._group_id + '/t1/group-' + \
                                           self._group_id + '_iteration-%d_template.nii*'
        templates_reader.inputs.iteration = range(1, 7)
        templates_reader.inputs.sort_filelist = False

        self.connect([
            (dartel_input_reader, self.input_node, [('out_files', 'dartel_input_images')]),
            (templates_reader, self.input_node, [('out_files', 'dartel_iteration_templates')])
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

        self.connect([
            (self.output_node, write_flowfields_node, [(('dartel_flow_fields', zip_nii, True), 'flow_fields')])
        ])

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipelines.
        """

        import nipype.pipeline.engine as npe
        import nipype.interfaces.utility as nutil
        from clinica.utils.io import unzip_nii
        import clinica.pipelines.t1_volume_existing_dartel.t1_volume_existing_dartel_utils as utils

        # Unzipping
        # =========
        unzip_dartel_input_node = npe.MapNode(nutil.Function(input_names=['in_file'],
                                                             output_names=['out_file'],
                                                             function=unzip_nii),
                                              name='unzip_dartel_input_node',
                                              iterfield=['in_file'])
        unzip_templates_node = npe.Node(nutil.Function(input_names=['in_file'],
                                                       output_names=['out_file'],
                                                       function=unzip_nii),
                                        name='unzip_templates_node')
        # DARTEL with existing template
        # =============================
        dartel_existing_template = npe.MapNode(utils.DARTELExistingTemplate(),
                                               name='dartel_existing_template',
                                               iterfield=['image_files'])
        if self.parameters['optimization_parameters'] is not None:
            dartel_existing_template.inputs.optimization_parameters = self.parameters['optimization_parameters']
        if self.parameters['regularization_form'] is not None:
            dartel_existing_template.inputs.regularization_form = self.parameters['regularization_form']

        # Connection
        # ==========
        self.connect([
            (self.input_node, unzip_dartel_input_node, [('dartel_input_images', 'in_file')]),
            (self.input_node, unzip_templates_node, [('dartel_iteration_templates', 'in_file')]),
            (unzip_dartel_input_node, dartel_existing_template, [(('out_file', utils.prepare_dartel_input_images),
                                                                  'image_files')]),
            (unzip_templates_node, dartel_existing_template, [(('out_file', utils.create_iteration_parameters,
                                                                self.parameters['iteration_parameters']),
                                                               'iteration_parameters')]),
            (dartel_existing_template, self.output_node, [('dartel_flow_fields', 'dartel_flow_fields')])
        ])
