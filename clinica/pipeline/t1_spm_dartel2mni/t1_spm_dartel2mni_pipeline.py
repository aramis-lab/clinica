"""T1 SPM Dartel2MNI - Clinica Pipeline.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details: https://gitlab.icm-institute.org/aramis/clinica/wikis/docs/InteractingWithClinica.
"""

import clinica.pipeline.engine as cpe


class T1SPMDartel2MNI(cpe.Pipeline):
    """T1 SPM Dartel2MNI SHORT DESCRIPTION.

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
        A clinica pipeline object containing the T1 SPM Dartel2MNI pipeline.

    Raises:


    Example:
        >>> from t1_spm_dartel2mni import T1SPMDartel2MNI
        >>> pipeline = T1SPMDartel2MNI('myGroup', '~/MYDATASET_BIDS', '~/MYDATASET_CAPS')
        >>> pipeline.parameters.update({
        >>>     # ...
        >>> })
        >>> pipeline.base_dir = '/tmp/'
        >>> pipeline.run()
    """

    def __init__(self, group_id, bids_directory=None, caps_directory=None, tsv_file=None, name=None):
        super(T1SPMDartel2MNI, self).__init__(bids_directory, caps_directory, tsv_file, name)

        self._group_id = group_id

        # Default parameters
        self._parameters = {'tissues': [1, 2, 3],
                            'bounding_box': None,
                            'voxel_size': None,
                            'modulation': True,
                            'fwhm': 0
                            }

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file.
        """
        pass

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Returns:
            A list of (string) input fields name.
        """

        return ['apply_to_files', 'flowfield_files', 'template_file']

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """

        return ['normalized_files']

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

        # DataGrabbers
        tissues_caps_reader = npe.MapNode(nio.DataGrabber(infields=['subject_id', 'session',
                                                                    'subject_repeat', 'session_repeat',
                                                                    'tissue'],
                                                          outfields=['out_files']),
                                          name="tissues_caps_reader",
                                          iterfield=['subject_id', 'session',
                                                     'subject_repeat', 'session_repeat'])
        tissues_caps_reader.inputs.base_directory = self.caps_directory
        tissues_caps_reader.inputs.template = 'subjects/%s/%s/t1/spm/segmentation/dartel_input/%s_%s_T1w_segm-%s_dartelinput.nii*'
        tissues_caps_reader.inputs.subject_id = self.subjects
        tissues_caps_reader.inputs.session = self.sessions
        tissues_caps_reader.inputs.tissue = [tissue_names[t] for t in self.parameters['tissues']]
        tissues_caps_reader.inputs.subject_repeat = self.subjects
        tissues_caps_reader.inputs.session_repeat = self.sessions
        tissues_caps_reader.inputs.sort_filelist = False

        flowfields_caps_reader = npe.Node(nio.DataGrabber(infields=['subject_id', 'session',
                                                                    'subject_repeat', 'session_repeat'],
                                                          outfields=['out_files']),
                                          name="flowfields_caps_reader")
        flowfields_caps_reader.inputs.base_directory = self.caps_directory
        flowfields_caps_reader.inputs.template = 'subjects/%s/%s/t1/spm/dartel/group-' + self._group_id + '/%s_%s_T1w_target-' + self._group_id + '_deformation.nii*'
        flowfields_caps_reader.inputs.subject_id = self.subjects
        flowfields_caps_reader.inputs.session = self.sessions
        flowfields_caps_reader.inputs.subject_repeat = self.subjects
        flowfields_caps_reader.inputs.session_repeat = self.sessions
        flowfields_caps_reader.inputs.sort_filelist = False

        template_caps_reader = npe.Node(nio.DataGrabber(infields=['group_id', 'group_id_repeat'],
                                                        outfields=['out_files']),
                                        name="template_caps_reader")
        template_caps_reader.inputs.base_directory = self.caps_directory
        template_caps_reader.inputs.template = 'groups/group-%s/t1/group-%s_template.nii*'
        template_caps_reader.inputs.group_id = self._group_id
        template_caps_reader.inputs.group_id_repeat = self._group_id

        self.connect([
            (tissues_caps_reader, self.input_node, [('out_files', 'apply_to_files')]),
            (flowfields_caps_reader, self.input_node, [('out_files', 'flowfield_files')]),
            (template_caps_reader, self.input_node, [('out_files', 'template_file')])
        ])

    def build_output_node(self):
        """Build and connect an output node to the pipeline.
        """

        import os.path as op
        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio
        import re
        from clinica.utils.io import zip_nii

        # Writing flowfields into CAPS
        # =======================
        write_normalized_node = npe.MapNode(name='write_normalized_node',
                                            iterfield=['container', 'normalized_files'],
                                            interface=nio.DataSink(infields=['normalized_files']))
        write_normalized_node.inputs.base_directory = self.caps_directory
        write_normalized_node.inputs.parameterization = False
        write_normalized_node.inputs.container = ['subjects/' + self.subjects[i] + '/' + self.sessions[i] +
                                                  '/t1/spm/dartel/group-' + self._group_id
                                                  for i in range(len(self.subjects))]
        write_normalized_node.inputs.regexp_substitutions = [
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
            (r'(.*)_deformation(\.nii(\.gz)?)$', r'\1\2'),
            # TODO Check which MNI space
            (r'(.*)mw(sub-.*)(\.nii(\.gz)?)$', r'\1\2_space-MNIXXX_modulated-on_probability\3'),
            (r'(.*)w(sub-.*)(\.nii(\.gz)?)$', r'\1\2_space-MNIXXX_modulated-off_probability\3'),
            (r'trait_added', r'')
        ]

        self.connect([
            (self.output_node, write_normalized_node, [(('normalized_files', zip_nii, True), 'normalized_files')])
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
        unzip_tissues_node = npe.MapNode(nutil.Function(input_names=['in_file'],
                                                        output_names=['out_file'],
                                                        function=unzip_nii),
                                         name='unzip_tissues_node',
                                         iterfield=['in_file'])
        unzip_flowfields_node = npe.MapNode(nutil.Function(input_names=['in_file'],
                                                           output_names=['out_file'],
                                                           function=unzip_nii),
                                            name='unzip_flowfields_node',
                                            iterfield=['in_file'])
        unzip_template_node = npe.Node(nutil.Function(input_names=['in_file'],
                                                      output_names=['out_file'],
                                                      function=unzip_nii),
                                       name='unzip_template_node')

        # DARTEL2MNI Registration
        # =======================
        dartel2mni_node = npe.MapNode(spm.DARTELNorm2MNI(),
                                      name='dartel2MNI',
                                      iterfield=['apply_to_files', 'flowfield_files'])

        if self.parameters['bounding_box'] is not None:
            dartel2mni_node.inputs.bounding_box = self.parameters['bounding_box']
        if self.parameters['voxel_size'] is not None:
            dartel2mni_node.inputs.voxel_size = self.parameters['voxel_size']
        #Modulation
        dartel2mni_node.inputs.modulate = self.parameters['modulate']
        #Smoothing
        dartel2mni_node.inputs.fwhm = self.parameters['fwhm']

        # Connection
        # ==========
        self.connect([
            # STEP 1
            (self.input_node, unzip_tissues_node, [('apply_to_files', 'in_file')]),
            (self.input_node, unzip_flowfields_node, [('flowfield_files', 'in_file')]),
            (self.input_node, unzip_template_node, [('template_file', 'in_file')]),
            (unzip_tissues_node, dartel2mni_node, [('out_file', 'apply_to_files')]),
            (unzip_flowfields_node, dartel2mni_node, [('out_file', 'flowfield_files')]),
            (unzip_template_node, dartel2mni_node, [('out_file', 'template_file')]),
            (dartel2mni_node, self.output_node, [('normalized_files', 'normalized_files')])
        ])
