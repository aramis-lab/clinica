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


class T1VolumeTissueSegmentation(cpe.Pipeline):
    """T1VolumeTissueSegmentation - Tissue segmentation, bias correction and
    spatial normalization to MNI space.

    Args:
        bids_directory: A BIDS directory.
        caps_directory: An empty output directory where CAPS structured data will be written.
        subjects_sessions_list: The Subjects-Sessions list file (in .tsv format).

    Returns:
        A clinica pipeline object containing the T1VolumeTissueSegmentation pipeline.

    Raises:

    """

    def __init__(self, bids_directory=None, caps_directory=None, tsv_file=None, name=None):
        super(T1VolumeTissueSegmentation, self).__init__(bids_directory, caps_directory, tsv_file, name)
        # Default parameters
        self._parameters = {'tissue_classes': [1, 2, 3],
                            'dartel_tissues': [1, 2, 3],
                            'tpm': None,
                            'save_warped_unmodulated': True,
                            'save_warped_modulated': False,
                            'affine_regularization': None,
                            'channel_info': None,
                            'sampling_distance': None,
                            'warping_regularization': None,
                            'write_deformation_fields': None,
                            'save_t1_mni': True
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

        return ['input_images']

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipelines.

        Returns:
            A list of (string) output fields name.
        """

        return ['bias_corrected_images',
                'bias_field_images',
                'dartel_input_images',
                'forward_deformation_field',
                'inverse_deformation_field',
                'modulated_class_images',
                'native_class_images',
                'normalized_class_images',
                'transformation_mat',
                't1_mni']

    def build_input_node(self):
        """Build and connect an input node to the pipelines.
        """

        import nipype.pipeline.engine as npe
        import nipype.interfaces.utility as nutil
        import clinica.pipelines.t1_volume_tissue_segmentation.t1_volume_tissue_segmentation_utils as utils

        # Reading BIDS
        # ============
        read_node = npe.Node(name="read_node",
                             interface=nutil.IdentityInterface(fields=['bids_images'],
                                                               mandatory_inputs=True))
        read_node.inputs.bids_images = utils.select_bids_images(self.subjects, self.sessions, 'T1w', self.bids_layout)

        self.connect([
            (read_node, self.input_node, [('bids_images', 'input_images')])
        ])

    def build_output_node(self):
        """Build and connect an output node to the pipelines.
        """

        import nipype.pipeline.engine as npe
        import clinica.pipelines.t1_volume_tissue_segmentation.t1_volume_tissue_segmentation_utils as utils
        import nipype.interfaces.io as nio
        from clinica.utils.io import zip_nii

        # Writing CAPS
        # ============
        datasink_infields = ['native_space', 'dartel_input']

        datasink_connections = [(('native_class_images', utils.group_nested_images_by_subject, True), 'native_space'),
                                (('dartel_input_images', utils.group_nested_images_by_subject, True), 'dartel_input')]

        if self.parameters['save_warped_unmodulated']:
            datasink_connections.append(
                (('normalized_class_images', utils.group_nested_images_by_subject, True), 'normalized'))
            datasink_infields.append('normalized')

        if self.parameters['save_warped_modulated']:
            datasink_connections.append(
                (('modulated_class_images', utils.group_nested_images_by_subject, True), 'modulated_normalized'))
            datasink_infields.append('modulated_normalized')

        if self.parameters['write_deformation_fields'] is not None:
            if self.parameters['write_deformation_fields'][0]:
                datasink_connections.append((('inverse_deformation_field', zip_nii, True), 'inverse_deformation_field'))
                datasink_infields.append('inverse_deformation_field')
            if self.parameters['write_deformation_fields'][1]:
                datasink_connections.append((('forward_deformation_field', zip_nii, True), 'forward_deformation_field'))
                datasink_infields.append('forward_deformation_field')

        if self.parameters['save_t1_mni'] is not None and self.parameters['save_t1_mni']:
            datasink_connections.append((('t1_mni', zip_nii, True), 't1_mni'))
            datasink_infields.append('t1_mni')

        datasink_iterfields = ['container'] + datasink_infields
        write_node = npe.MapNode(name='WritingCAPS',
                                 iterfield=datasink_iterfields,
                                 interface=nio.DataSink(infields=datasink_infields))
        write_node.inputs.base_directory = self.caps_directory
        write_node.inputs.parameterization = False
        write_node.inputs.container = ['subjects/' + self.subjects[i] + '/' + self.sessions[i] + '/t1/spm/segmentation'
                                       for i in range(len(self.subjects))]

        write_node.inputs.regexp_substitutions = [
            (r'(.*)c1(sub-.*)(\.nii(\.gz)?)$', r'\1\2_segm-graymatter\3'),
            (r'(.*)c2(sub-.*)(\.nii(\.gz)?)$', r'\1\2_segm-whitematter\3'),
            (r'(.*)c3(sub-.*)(\.nii(\.gz)?)$', r'\1\2_segm-csf\3'),
            (r'(.*)c4(sub-.*)(\.nii(\.gz)?)$', r'\1\2_segm-bone\3'),
            (r'(.*)c5(sub-.*)(\.nii(\.gz)?)$', r'\1\2_segm-softtissue\3'),
            (r'(.*)c6(sub-.*)(\.nii(\.gz)?)$', r'\1\2_segm-background\3'),
            (r'(.*)(/native_space/sub-.*)(\.nii(\.gz)?)$', r'\1\2_probability\3'),
            (r'(.*)(/([a-z]+)_deformation_field/)i?y_(sub-.*)(\.nii(\.gz)?)$',
             r'\1/normalized_space/\4_target-Ixi549Space_transformation-\3_deformation\5'),
            (r'(.*)(/t1_mni/)w(sub-.*)_T1w(\.nii(\.gz)?)$', r'\1/normalized_space/\3_space-Ixi549Space_T1w\4'),
            (r'(.*)(/modulated_normalized/)mw(sub-.*)(\.nii(\.gz)?)$',
             r'\1/normalized_space/\3_space-Ixi549Space_modulated-on_probability\4'),
            (r'(.*)(/normalized/)w(sub-.*)(\.nii(\.gz)?)$',
             r'\1/normalized_space/\3_space-Ixi549Space_modulated-off_probability\4'),
            (r'(.*/dartel_input/)r(sub-.*)(\.nii(\.gz)?)$', r'\1\2_dartelinput\3'),
            (r'trait_added', r'')
        ]

        self.connect([
            # Writing CAPS
            (self.output_node, write_node, datasink_connections)
        ])

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipelines.
        """

        import os
        import platform
        import nipype.interfaces.spm as spm
        import nipype.interfaces.matlab as mlab
        import nipype.pipeline.engine as npe
        import nipype.interfaces.utility as nutil
        import clinica.pipelines.t1_volume_tissue_segmentation.t1_volume_tissue_segmentation_utils as utils
        from clinica.utils.io import unzip_nii

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
                    if platform.system() == 'Darwin':
                        tissue_map = os.path.join(str(spm_home), 'spm12.app/Contents/MacOS/spm12_mcr/spm12/spm12/tpm/TPM.nii')
                    else:
                        tissue_map = os.path.join(str(spm_home), 'spm12_mcr/spm/spm12/tpm/TPM.nii')
                else:
                    raise RuntimeError('SPM standalone version not supported. Please upgrade SPM standalone.')
        else:
            raise RuntimeError('SPM could not be found. Please verify your SPM_HOME environment variable.')

        # Unzipping
        # ===============================
        unzip_node = npe.MapNode(nutil.Function(input_names=['in_file'],
                                                output_names=['out_file'],
                                                function=unzip_nii),
                                 name='unzip_node', iterfield=['in_file'])

        # Unified Segmentation
        # ===============================
        new_segment = npe.MapNode(spm.NewSegment(),
                                  name='new_segment',
                                  iterfield=['channel_files'])

        if self.parameters['affine_regularization'] is not None:
            new_segment.inputs.affine_regularization = self.parameters['affine_regularization']
        if self.parameters['channel_info'] is not None:
            new_segment.inputs.channel_info = self.parameters['channel_info']
        if self.parameters['sampling_distance'] is not None:
            new_segment.inputs.sampling_distance = self.parameters['sampling_distance']
        if self.parameters['warping_regularization'] is not None:
            new_segment.inputs.warping_regularization = self.parameters['warping_regularization']

        # Check if we need to save the forward transformation for registering the T1 to the MNI space
        if self.parameters['save_t1_mni'] is not None and self.parameters['save_t1_mni']:
            if self.parameters['write_deformation_fields'] is not None:
                self.parameters['write_deformation_fields'][1] = True
            else:
                self.parameters['write_deformation_fields'] = [False, True]

        if self.parameters['write_deformation_fields'] is not None:
            new_segment.inputs.write_deformation_fields = self.parameters['write_deformation_fields']

        if self.parameters['tpm'] is not None:
            tissue_map = self.parameters['tpm']

        new_segment.inputs.tissues = utils.get_tissue_tuples(tissue_map,
                                                             self.parameters['tissue_classes'],
                                                             self.parameters['dartel_tissues'],
                                                             self.parameters['save_warped_unmodulated'],
                                                             self.parameters['save_warped_modulated'])

        # Apply segmentation deformation to T1 (into MNI space)
        # ========================================================
        if self.parameters['save_t1_mni'] is not None and self.parameters['save_t1_mni']:

            t1_to_mni = npe.MapNode(utils.ApplySegmentationDeformation(),
                                    name='t1_to_mni',
                                    iterfield=['deformation_field', 'in_files'])
            self.connect([
                (unzip_node, t1_to_mni, [('out_file', 'in_files')]),
                (new_segment, t1_to_mni, [('forward_deformation_field', 'deformation_field')]),
                (t1_to_mni, self.output_node, [('out_files', 't1_mni')])
            ])

        # Connection
        # ==========
        self.connect([
            (self.input_node, unzip_node, [('input_images', 'in_file')]),
            (unzip_node, new_segment, [('out_file', 'channel_files')]),
            (new_segment, self.output_node, [('bias_corrected_images', 'bias_corrected_images'),
                                             ('bias_field_images', 'bias_field_images'),
                                             ('dartel_input_images', 'dartel_input_images'),
                                             ('forward_deformation_field', 'forward_deformation_field'),
                                             ('inverse_deformation_field', 'inverse_deformation_field'),
                                             ('modulated_class_images', 'modulated_class_images'),
                                             ('native_class_images', 'native_class_images'),
                                             ('normalized_class_images', 'normalized_class_images'),
                                             ('transformation_mat', 'transformation_mat')])
        ])
