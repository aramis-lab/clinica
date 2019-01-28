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


class T1VolumeDartel2MNI(cpe.Pipeline):
    """T1VolumeDartel2MNI - Dartel template to MNI.

    Args:
        input_dir: A BIDS directory.
        output_dir: An empty output directory where CAPS structured data will
        be written.  subjects_sessions_list: The Subjects-Sessions list file
        (in .tsv format).

    Returns:
        A clinica pipeline object containing the T1VolumeDartel2MNI pipeline.

    Raises:


    """

    def __init__(self, bids_directory=None, caps_directory=None, tsv_file=None, name=None, group_id='default'):
        import os
        super(T1VolumeDartel2MNI, self).__init__(bids_directory, caps_directory, tsv_file, name)

        if not group_id.isalnum():
            raise ValueError('Not valid group_id value. It must be composed only by letters and/or numbers')

        self._group_id = group_id
        # Check that group already exists
        if not os.path.exists(os.path.join(os.path.abspath(caps_directory), 'groups', 'group-' + group_id)):
            error_message = group_id \
                            + ' does not exists, please choose an other one (or maybe you need to run t1-spm-dartel).' \
                            + '\nGroups that already exists in your CAPS directory are : \n'
            list_groups = os.listdir(os.path.join(os.path.abspath(caps_directory), 'groups'))
            is_empty = True
            for e in list_groups:
                if e.startswith('group-'):
                    error_message += e + ' \n'
                    is_empty = False
            if is_empty is True:
                error_message += 'NO GROUP FOUND'
            raise ValueError(error_message)

        # Default parameters
        self._parameters = {'tissues': [1, 2, 3],
                            'bounding_box': None,
                            'voxel_size': None,
                            'modulation': True,
                            'fwhm': [8],
                            'atlas_list': ['AAL2', 'LPBA40', 'Neuromorphometrics', 'AICHA', 'Hammers']
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

        return ['apply_to_files', 'flowfield_files', 'template_file']

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipelines.

        Returns:
            A list of (string) output fields name.
        """

        return ['normalized_files', 'smoothed_normalized_files', 'atlas_statistics']

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

        # Segmented Tissues DataGrabber
        # ==============================
        tissues_caps_reader = npe.MapNode(nio.DataGrabber(infields=['subject_id', 'session',
                                                                    'subject_repeat', 'session_repeat',
                                                                    'tissue'],
                                                          outfields=['out_files']),
                                          name="tissues_caps_reader",
                                          iterfield=['subject_id', 'session',
                                                     'subject_repeat', 'session_repeat'])
        tissues_caps_reader.inputs.base_directory = self.caps_directory
        tissues_caps_reader.inputs.template = 'subjects/%s/%s/t1/spm/segmentation/native_space/%s_%s_T1w_segm-%s_probability.nii*'
        tissues_caps_reader.inputs.subject_id = self.subjects
        tissues_caps_reader.inputs.session = self.sessions
        tissues_caps_reader.inputs.tissue = [tissue_names[t] for t in self.parameters['tissues']]
        tissues_caps_reader.inputs.subject_repeat = self.subjects
        tissues_caps_reader.inputs.session_repeat = self.sessions
        tissues_caps_reader.inputs.sort_filelist = False

        # Flow Fields DataGrabber
        # ==============================
        flowfields_caps_reader = npe.Node(nio.DataGrabber(infields=['subject_id', 'session',
                                                                    'subject_repeat', 'session_repeat'],
                                                          outfields=['out_files']),
                                          name="flowfields_caps_reader")
        flowfields_caps_reader.inputs.base_directory = self.caps_directory
        flowfields_caps_reader.inputs.template = 'subjects/%s/%s/t1/spm/dartel/group-' + self._group_id \
                                                 + '/%s_%s_T1w_target-' + self._group_id \
                                                 + '_transformation-forward_deformation.nii*'
        flowfields_caps_reader.inputs.subject_id = self.subjects
        flowfields_caps_reader.inputs.session = self.sessions
        flowfields_caps_reader.inputs.subject_repeat = self.subjects
        flowfields_caps_reader.inputs.session_repeat = self.sessions
        flowfields_caps_reader.inputs.sort_filelist = False

        # Dartel Template DataGrabber
        # ===========================
        template_caps_reader = npe.Node(nio.DataGrabber(infields=['group_id', 'group_id_repeat'],
                                                        outfields=['out_files']),
                                        name="template_caps_reader")
        template_caps_reader.inputs.base_directory = self.caps_directory
        template_caps_reader.inputs.template = 'groups/group-%s/t1/group-%s_template.nii*'
        template_caps_reader.inputs.group_id = self._group_id
        template_caps_reader.inputs.group_id_repeat = self._group_id
        template_caps_reader.inputs.sort_filelist = False

        self.connect([
            (tissues_caps_reader, self.input_node, [('out_files', 'apply_to_files')]),
            (flowfields_caps_reader, self.input_node, [('out_files', 'flowfield_files')]),
            (template_caps_reader, self.input_node, [('out_files', 'template_file')])
        ])

    def build_output_node(self):
        """Build and connect an output node to the pipelines.
        """
        import os.path as op
        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio
        import re
        from clinica.utils.io import zip_nii

        # Writing normalized images (and smoothed) into CAPS
        # ==================================================
        write_normalized_node = npe.MapNode(name='write_normalized_node',
                                            iterfield=['container', 'normalized_files', 'smoothed_normalized_files'],
                                            interface=nio.DataSink(infields=['normalized_files',
                                                                             'smoothed_normalized_files']))
        write_normalized_node.inputs.base_directory = self.caps_directory
        write_normalized_node.inputs.parameterization = False
        write_normalized_node.inputs.container = ['subjects/' + self.subjects[i] + '/' + self.sessions[i] +
                                                  '/t1/spm/dartel/group-' + self._group_id
                                                  for i in range(len(self.subjects))]
        write_normalized_node.inputs.regexp_substitutions = [
            (r'(.*)c1(sub-.*)(\.nii(\.gz)?)$', r'\1\2_segm-graymatter_probability\3'),
            (r'(.*)c2(sub-.*)(\.nii(\.gz)?)$', r'\1\2_segm-whitematter_probability\3'),
            (r'(.*)c3(sub-.*)(\.nii(\.gz)?)$', r'\1\2_segm-csf_probability\3'),
            (r'(.*)c4(sub-.*)(\.nii(\.gz)?)$', r'\1\2_segm-bone_probability\3'),
            (r'(.*)c5(sub-.*)(\.nii(\.gz)?)$', r'\1\2_segm-softtissue_probability\3'),
            (r'(.*)c6(sub-.*)(\.nii(\.gz)?)$', r'\1\2_segm-background_probability\3'),
            (r'(.*)mw(sub-.*)_probability(\.nii(\.gz)?)$', r'\1\2_space-Ixi549Space_modulated-on_probability\3'),
            (r'(.*)w(sub-.*)_probability(\.nii(\.gz)?)$', r'\1\2_space-Ixi549Space_modulated-off_probability\3'),
            (r'(.*)/normalized_files/(sub-.*)$', r'\1/\2'),
            (r'(.*)/smoothed_normalized_files/(fwhm-[0-9]+mm)_(sub-.*)_probability(\.nii(\.gz)?)$',
             r'\1/\3_\2_probability\4'),
            (r'trait_added', r'')
        ]

        # Writing atlas statistics into CAPS
        # ==================================
        write_atlas_node = npe.MapNode(name='write_atlas_node',
                                            iterfield=['container', 'atlas_statistics'],
                                            interface=nio.DataSink(infields=['atlas_statistics']))
        write_atlas_node.inputs.base_directory = self.caps_directory
        write_atlas_node.inputs.parameterization = False
        write_atlas_node.inputs.container = ['subjects/' + self.subjects[i] + '/' + self.sessions[i] +
                                             '/t1/spm/dartel/group-' + self._group_id + '/atlas_statistics'
                                             for i in range(len(self.subjects))]
        write_atlas_node.inputs.regexp_substitutions = [
            (r'(.*atlas_statistics)/atlas_statistics/mwc1(sub-.*)(_space-.*_map-graymatter_statistics\.tsv)$',
             r'\1/\2\3'),
            (r'(.*atlas_statistics)/atlas_statistics/(m?w)?(sub-.*_T1w).*(_space-.*_map-graymatter_statistics).*(\.tsv)$',
             r'\1/\3\4\5'),
            (r'trait_added', r'')
        ]

        self.connect([
            (self.output_node, write_normalized_node, [(('normalized_files', zip_nii, True), 'normalized_files'),
                                                       (('smoothed_normalized_files', zip_nii, True),
                                                        'smoothed_normalized_files')]),
            (self.output_node, write_atlas_node, [('atlas_statistics', 'atlas_statistics')])
        ])

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipelines.
        """

        import os
        import nipype.interfaces.spm as spm
        import nipype.interfaces.matlab as mlab
        import nipype.pipeline.engine as npe
        import nipype.interfaces.utility as nutil
        from clinica.utils.io import unzip_nii
        from clinica.pipelines.t1_volume_dartel2mni.t1_volume_dartel2mni_utils import prepare_flowfields, join_smoothed_files, atlas_statistics, select_gm_images

        spm_home = os.getenv("SPM_HOME")
        mlab_home = os.getenv("MATLABCMD")
        mlab.MatlabCommand.set_default_matlab_cmd(mlab_home)
        mlab.MatlabCommand.set_default_paths(spm_home)

        if 'SPMSTANDALONE_HOME' in os.environ:
            if 'MCR_HOME' in os.environ:
                matlab_cmd = os.path.join(os.environ['SPMSTANDALONE_HOME'],
                                          'run_spm12.sh') \
                             + ' ' + os.environ['MCR_HOME'] \
                             + ' script'
                spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)
                version = spm.SPMCommand().version
            else:
                raise EnvironmentError('MCR_HOME variable not in environnement. Althought, '
                                       + 'SPMSTANDALONE_HOME has been found')
        else:
            version = spm.Info.getinfo()

        if version:
            if isinstance(version, dict):
                spm_path = version['path']
                if version['name'] == 'SPM8':
                    print('You are using SPM version 8. The recommended version to use with Clinica is SPM 12. '
                          + 'Please upgrade your SPM toolbox.')
                    tissue_map = os.path.join(spm_path, 'toolbox/Seg/TPM.nii')
                elif version['name'] == 'SPM12':
                    tissue_map = os.path.join(spm_path, 'tpm/TPM.nii')
                else:
                    raise RuntimeError('SPM version 8 or 12 could not be found. Please upgrade your SPM toolbox.')
            if isinstance(version, str):
                if float(version) >= 12.7169:
                    tissue_map = os.path.join(str(spm_home), 'spm12_mcr/spm/spm12/tpm/TPM.nii')
                else:
                    raise RuntimeError('SPM standalone version not supported. Please upgrade SPM standalone.')
        else:
            raise RuntimeError('SPM could not be found. Please verify your SPM_HOME environment variable.')

        # Unzipping
        # =========
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
        dartel2mni_node.inputs.modulate = self.parameters['modulation']
        dartel2mni_node.inputs.fwhm = 0

        # Smoothing
        # =========
        if self.parameters['fwhm'] is not None and len(self.parameters['fwhm']) > 0:
            smoothing_node = npe.MapNode(spm.Smooth(),
                                         name='smoothing_node',
                                         iterfield=['in_files'])

            smoothing_node.iterables = [('fwhm', [[x, x, x] for x in self.parameters['fwhm']]),
                                        ('out_prefix', ['fwhm-' + str(x) + 'mm_' for x in self.parameters['fwhm']])]
            smoothing_node.synchronize = True

            join_smoothing_node = npe.JoinNode(interface=nutil.Function(input_names=['smoothed_normalized_files'],
                                                                        output_names=['smoothed_normalized_files'],
                                                                        function=join_smoothed_files),
                                               joinsource='smoothing_node',
                                               joinfield='smoothed_normalized_files',
                                               name='join_smoothing_node')
            self.connect([
                (dartel2mni_node, smoothing_node, [('normalized_files', 'in_files')]),
                (smoothing_node, join_smoothing_node, [('smoothed_files', 'smoothed_normalized_files')]),
                (join_smoothing_node, self.output_node, [('smoothed_normalized_files', 'smoothed_normalized_files')])
            ])
        else:
            self.output_node.inputs.smoothed_normalized_files = []

        # Atlas Statistics
        # ================
        atlas_stats_node = npe.MapNode(nutil.Function(input_names=['in_image',
                                                                   'in_atlas_list'],
                                                      output_names=['atlas_statistics'],
                                                      function=atlas_statistics),
                                       name='atlas_stats_node',
                                       iterfield=['in_image'])
        atlas_stats_node.inputs.in_atlas_list = self.parameters['atlas_list']

        # Connection
        # ==========
        self.connect([
            (self.input_node, unzip_tissues_node, [('apply_to_files', 'in_file')]),
            (self.input_node, unzip_flowfields_node, [('flowfield_files', 'in_file')]),
            (self.input_node, unzip_template_node, [('template_file', 'in_file')]),
            (unzip_tissues_node, dartel2mni_node, [('out_file', 'apply_to_files')]),
            (unzip_flowfields_node, dartel2mni_node, [(('out_file', prepare_flowfields, self.parameters['tissues']),
                                                       'flowfield_files')]),
            (unzip_template_node, dartel2mni_node, [('out_file', 'template_file')]),
            (dartel2mni_node, self.output_node, [('normalized_files', 'normalized_files')]),
            (dartel2mni_node, atlas_stats_node, [(('normalized_files', select_gm_images), 'in_image')]),
            (atlas_stats_node, self.output_node, [('atlas_statistics', 'atlas_statistics')])
        ])
