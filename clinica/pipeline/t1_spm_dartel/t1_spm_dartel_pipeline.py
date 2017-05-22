"""T1 SPM Dartel - Clinica Pipeline.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details: https://gitlab.icm-institute.org/aramis/clinica/wikis/docs/InteractingWithClinica.
"""

import clinica.pipeline.engine as cpe


class T1SPMDartel(cpe.Pipeline):
    """T1 SPM Dartel SHORT DESCRIPTION.

    Warnings:
        - A WARNING.

    Todos:
        - [x] A FILLED TODO ITEM.
        - [ ] AN ON-GOING TODO ITEM.

    Args:
        input_dir: A BIDS directory.
        output_dir: An empty output directory where CAPS structured data will be written.
        subjects_sessions_list: The Subjects-Sessions list file (in .tsv format).

    Returns:
        A clinica pipeline object containing the T1 SPM Dartel pipeline.

    Raises:


    Example:
        >>> from t1_spm_dartel import T1SPMDartel
        >>> pipeline = T1SPMDartel('~/MYDATASET_BIDS', '~/MYDATASET_CAPS')
        >>> pipeline.parameters.update({
        >>>     # ...
        >>> })
        >>> pipeline.base_dir = '/tmp/'
        >>> pipeline.run()
    """

    def __init__(self, group_id, bids_directory=None, caps_directory=None, tsv_file=None, name=None):
        super(T1SPMDartel, self).__init__(bids_directory, caps_directory, tsv_file, name)

        self._group_id = group_id
        # Default parameters
        self._parameters = {'dartel_tissues': [1, 2, 3],
                            'iteration_parameters': None,
                            'optimization_parameters': None,
                            'regularization_form': None,
                            'template_prefix': None
                            }

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Returns:
            A list of (string) input fields name.
        """

        return ['dartel_input_images']

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """

        return ['final_template_file', 'template_files', 'dartel_flow_fields']

    def build_input_node(self):
        """Build and connect an input node to the pipeline.
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

        # DataGrabber
        caps_reader = npe.MapNode(nio.DataGrabber(infields=['subject_id', 'session',
                                                            'subject_repeat', 'session_repeat',
                                                            'tissue'],
                                                  outfields=['out_files']),
                                  name="caps_reader", iterfield=['tissue'])

        caps_reader.inputs.base_directory = self.caps_directory
        caps_reader.inputs.template = 'subjects/%s/%s/t1/spm/segmentation/dartel_input/%s_%s_T1w_segm-%s_dartelinput.nii*'
        caps_reader.inputs.subject_id = self.subjects
        caps_reader.inputs.session = self.sessions
        caps_reader.inputs.tissue = [tissue_names[t] for t in self.parameters['dartel_tissues']]
        caps_reader.inputs.subject_repeat = self.subjects
        caps_reader.inputs.session_repeat = self.sessions
        caps_reader.inputs.sort_filelist = False

        self.connect([
            (caps_reader, self.input_node, [('out_files', 'dartel_input_images')])
        ])

    def build_output_node(self):
        """Build and connect an output node to the pipeline.
        """
        import os.path as op
        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio
        import re
        from clinica.utils.io import zip_nii

        # Writing flowfields CAPS
        # =======================
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
            (r'(.*)u_(sub-.*)(\.nii(\.gz)?)$', r'\1\2_target-' + re.escape(self._group_id) + r'_deformation\3'),
            (r'trait_added', r'')
        ]

        # Writing templates CAPS
        # ======================
        write_template_node = npe.Node(nio.DataSink(), name='write_template_node')
        write_template_node.inputs.parameterization = False
        write_template_node.inputs.base_directory = op.join(self.caps_directory, 'groups/group-' + self._group_id, 't1')
        write_template_node.inputs.regexp_substitutions = [
            (r'(.*)final_template_file/.*(\.nii(\.gz)?)$', r'\1group-' + re.escape(self._group_id) + r'_template\2'),
            (r'(.*)template_files/.*([0-9])(\.nii(\.gz)?)$', r'\1group-' + re.escape(self._group_id) + r'_iteration-\2_template\3')
        ]

        self.connect([
            (self.output_node, write_flowfields_node, [(('dartel_flow_fields', zip_nii, True), 'flow_fields')]),
            (self.output_node, write_template_node, [(('final_template_file', zip_nii, True), 'final_template_file'),
                                                     (('template_files', zip_nii, True), 'template_files')])
        ])

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline.
        """

        import os
        import nipype.interfaces.spm as spm
        import nipype.interfaces.matlab as mlab
        import nipype.pipeline.engine as npe
        import nipype.interfaces.utility as nutil
        from clinica.utils.io import unzip_nii

        spm_home = os.getenv("SPM_HOME")
        mlab_home = os.getenv("MATLABCMD")
        mlab.MatlabCommand.set_default_matlab_cmd(mlab_home)
        mlab.MatlabCommand.set_default_paths(spm_home)

        version = spm.Info.version()

        if version:
            if version['name'] == 'SPM8':
                print 'You are using SPM version 8. The recommended version to use with Clinica is SPM 12. ' \
                      'Please upgrade your SPM toolbox.'
            elif version['name'] != 'SPM12':
                raise RuntimeError('SPM version 8 or 12 could not be found. Please upgrade your SPM toolbox.')
        else:
            raise RuntimeError('SPM could not be found. Please verify your SPM_HOME environment variable.')

        # Unzipping
        # ===============================
        unzip_node = npe.MapNode(nutil.Function(input_names=['in_file'],
                                                output_names=['out_file'],
                                                function=unzip_nii),
                                 name='unzip_node', iterfield=['in_file'])

        # DARTEL template
        # ===============================
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
            (self.unzip_node, dartel_template,    [('out_file', 'image_files')]),
            (dartel_template, self.output_node, [('dartel_flow_fields', 'dartel_flow_fields'),
                                                 ('final_template_file', 'final_template_file'),
                                                 ('template_files', 'template_files')])
        ])
