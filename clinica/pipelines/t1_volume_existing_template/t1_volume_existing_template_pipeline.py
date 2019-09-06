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


# Use hash instead of parameters for iterables folder names
# Otherwise path will be too long and generate OSError
from nipype import config
cfg = dict(execution={'parameterize_dirs': False})
config.update_config(cfg)


class T1VolumeExistingTemplate(cpe.Pipeline):
    """T1VolumeExistingTemplate

    Args:
        input_dir: A BIDS directory.
        output_dir: An empty output directory where CAPS structured data will be written.
        subjects_sessions_list: The Subjects-Sessions list file (in .tsv format).

    Returns:
        A clinica pipeline object containing the T1VolumeExistingTemplate pipeline.
    """
    def __init__(self, bids_directory=None, caps_directory=None, tsv_file=None, name=None, group_id='default'):
        from os.path import exists, join, abspath
        from os import listdir

        super(T1VolumeExistingTemplate, self).__init__(bids_directory, caps_directory, tsv_file, name)

        # Check that group already exists
        if not exists(join(abspath(caps_directory), 'groups', 'group-' + group_id)):
            error_message = 'group_id : ' + group_id + ' does not exists, ' \
                            + 'please choose another one. Groups that exist' \
                            + ' in your CAPS directory are: \n'
            list_groups = listdir(join(abspath(caps_directory), 'groups'))
            has_one_group = False
            for e in list_groups:
                if e.startswith('group-'):
                    error_message += e + ' \n'
                    has_one_group = True
            if not has_one_group:
                error_message = error_message + 'No group found ! ' \
                                + 'Use t1-volume pipeline if you do not ' \
                                + 'have a template yet!'
            raise ValueError(error_message)

        self._group_id = group_id

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
                            'save_t1_mni': False,
                            'iteration_parameters': None,
                            'optimization_parameters': None,
                            'regularization_form': None,
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

        return ['input_images', 'dartel_iteration_templates', 'dartel_final_template']

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
                't1_mni',
                'dartel_flow_fields',
                'normalized_files',
                'smoothed_normalized_files',
                'atlas_statistics']

    def build_input_node(self):
        """Build and connect an input node to the pipelines."""
        import nipype.pipeline.engine as npe
        import nipype.interfaces.utility as nutil
        import clinica.pipelines.t1_volume_tissue_segmentation.t1_volume_tissue_segmentation_utils as seg_utils
        from os.path import join
        import glob

        # Reading T1w
        # ============
        t1w_images = seg_utils.select_bids_images(self.subjects,
                                                  self.sessions,
                                                  'T1w',
                                                  self.bids_layout)

        # Dartel Iterations Templates
        # ============================
        g_id = self._group_id
        pattern_iter_dartel = join(self.caps_directory, 'groups', 'group-' + g_id, 't1', 'group-' + g_id + '_iteration-*_template.nii*')
        iter_template = glob.glob(pattern_iter_dartel)

        # Dartel Template
        # ================
        pattern_final_dartel = join(self.caps_directory, 'groups', 'group-' + g_id, 't1', 'group-' + g_id + '_template.nii*')
        final_template = glob.glob(pattern_final_dartel)
        if len(final_template) != 1:
            raise FileNotFoundError('Could find template !')
        else:
            final_template = final_template[0]

        read_node = npe.Node(name="read_node",
                             interface=nutil.IdentityInterface(fields=['t1w',
                                                                       'templates_iter',
                                                                       'final_template'],
                                                               mandatory_inputs=True),
                             iterables=[('t1w', t1w_images)],
                             synchronize=True)

        read_node.inputs.templates_iter = iter_template
        read_node.inputs.final_template = final_template

        self.connect([
            (read_node, self.input_node, [('t1w', 'input_images')]),
            (read_node, self.input_node, [('templates_iter', 'dartel_iteration_templates')]),
            (read_node, self.input_node, [('final_template', 'dartel_final_template')])
        ])

    def build_output_node(self):
        """Build and connect an output node to the pipelines.
        """
        import re
        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio
        import nipype.interfaces.utility as nutil
        from clinica.utils.io import zip_nii
        from ..t1_volume_tissue_segmentation import t1_volume_tissue_segmentation_utils as seg_utils
        from ..t1_volume_existing_template import t1_volume_existing_template_utils as existing_template_utils

        # Writing Segmentation output into CAPS
        # =====================================
        datasink_infields = ['native_space', 'dartel_input']

        datasink_connections = [(('native_class_images', seg_utils.group_nested_images_by_subject, True),
                                 'native_space'),
                                (('dartel_input_images', seg_utils.group_nested_images_by_subject, True),
                                 'dartel_input')]

        if self.parameters['save_warped_unmodulated']:
            datasink_connections.append(
                (('normalized_class_images', seg_utils.group_nested_images_by_subject, True), 'normalized'))
            datasink_infields.append('normalized')

        if self.parameters['save_warped_modulated']:
            datasink_connections.append(
                (('modulated_class_images', seg_utils.group_nested_images_by_subject, True), 'modulated_normalized'))
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

        write_segmentation_node = npe.MapNode(name='write_segmentation_node',
                                              iterfield=datasink_infields,
                                              interface=nio.DataSink(infields=datasink_infields))
        write_segmentation_node.inputs.base_directory = self.caps_directory
        write_segmentation_node.inputs.parameterization = False
        write_segmentation_node.inputs.regexp_substitutions = [
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

        extract_container_node_write_segmentation = npe.Node(nutil.Function(input_names=['filename'],
                                                                            output_names=['container_path'],
                                                                            function=existing_template_utils.container_name_for_write_segmentation),
                                                             name='extract_container_node_write_segmentation')

        self.connect([
            (self.output_node, write_segmentation_node, datasink_connections),
            (self.input_node, extract_container_node_write_segmentation, [('input_images', 'filename')]),
            (extract_container_node_write_segmentation, write_segmentation_node, [('container_path', 'container')])
        ])

        # Writing flowfields into CAPS
        # ============================
        write_flowfields_node = npe.MapNode(name='write_flowfields_node',
                                            iterfield=['flow_fields'],
                                            interface=nio.DataSink(infields=['flow_fields']))
        write_flowfields_node.inputs.base_directory = self.caps_directory
        write_flowfields_node.inputs.parameterization = False

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
        extract_container_node_write_flowfields = npe.Node(nutil.Function(input_names=['filename', 'group_id'],
                                                                          output_names=['container_path'],
                                                                          function=existing_template_utils.container_name_for_write_normalized),
                                                           name='extract_container_node_write_flowfields')
        extract_container_node_write_flowfields.inputs.group_id = self._group_id
        self.connect([
            (self.output_node, write_flowfields_node, [(('dartel_flow_fields', zip_nii, True), 'flow_fields')]),
            (self.input_node, extract_container_node_write_flowfields, [('input_images', 'filename')]),
            (extract_container_node_write_flowfields, write_flowfields_node, [('container_path', 'container')])
        ])

        # Writing normalized images (and smoothed) into CAPS
        # ==================================================
        write_normalized_node = npe.MapNode(name='write_normalized_node',
                                            iterfield=['normalized_files', 'smoothed_normalized_files'],
                                            interface=nio.DataSink(infields=['normalized_files',
                                                                             'smoothed_normalized_files']))
        write_normalized_node.inputs.base_directory = self.caps_directory
        write_normalized_node.inputs.parameterization = False

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

        extract_container_node_write_normalized = npe.Node(nutil.Function(input_names=['filename', 'group_id'],
                                                                          output_names=['container_path'],
                                                                          function=existing_template_utils.container_name_for_write_normalized),
                                                           name='extract_container_node_write_normalized')
        extract_container_node_write_normalized.inputs.group_id = self._group_id

        # Writing atlases statistics into CAPS
        # ==================================================
        write_atlas_node = npe.MapNode(name='write_atlas_node',
                                            iterfield=['atlas_statistics'],
                                            interface=nio.DataSink(infields=['atlas_statistics']))
        write_atlas_node.inputs.base_directory = self.caps_directory
        write_atlas_node.inputs.parameterization = False
        write_atlas_node.inputs.regexp_substitutions = [
            (r'(.*atlas_statistics)/atlas_statistics/mwc1(sub-.*_T1w).*(_space-.*_map-graymatter_statistics\.tsv)$',
             r'\1/\2\3'),
            (r'(.*atlas_statistics)/atlas_statistics/(m?w)?(sub-.*_T1w).*(_space-.*_map-graymatter_statistics).*(\.tsv)$',
             r'\1/\3\4\5'),
            (r'trait_added', r'')
        ]

        extract_container_node_atlas = npe.Node(nutil.Function(input_names=['filename', 'group_id'],
                                                               output_names=['container_path'],
                                                               function=existing_template_utils.container_name_for_atlas),
                                                name='extract_container_node_atlas')
        extract_container_node_atlas.inputs.group_id = self._group_id
        self.connect([
            (self.output_node, write_normalized_node, [(('normalized_files', zip_nii, True), 'normalized_files'),
                                                       (('smoothed_normalized_files', zip_nii, True),
                                                        'smoothed_normalized_files')]),
            (self.input_node, extract_container_node_write_normalized, [('input_images', 'filename')]),
            (extract_container_node_write_normalized, write_normalized_node, [('container_path', 'container')]),
            (self.output_node, write_atlas_node, [('atlas_statistics', 'atlas_statistics')]),
            (self.input_node, extract_container_node_atlas, [('input_images', 'filename')]),
            (extract_container_node_atlas, write_atlas_node, [('container_path', 'container')])
        ])

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipelines.
        """
        import nipype.interfaces.spm as spm
        import nipype.pipeline.engine as npe
        import nipype.interfaces.utility as nutil
        from ..t1_volume_tissue_segmentation import t1_volume_tissue_segmentation_utils as seg_utils
        from ..t1_volume_existing_dartel import t1_volume_existing_dartel_utils as existing_dartel_utils
        from ..t1_volume_dartel2mni import t1_volume_dartel2mni_utils as dartel2mni_utils
        from ..t1_volume_parcellation import t1_volume_parcellation_utils as parcellation_utils
        from clinica.utils.io import unzip_nii
        from clinica.utils.spm import get_tpm

        # Get Tissue Probability Map from SPM
        tissue_map = get_tpm()

        # Unzipping Images
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

        new_segment.inputs.tissues = seg_utils.get_tissue_tuples(tissue_map,
                                                                 self.parameters['tissue_classes'],
                                                                 self.parameters['dartel_tissues'],
                                                                 self.parameters['save_warped_unmodulated'],
                                                                 self.parameters['save_warped_modulated'])

        # Apply segmentation deformation to T1 (into MNI space)
        # ========================================================
        if self.parameters['save_t1_mni'] is not None and self.parameters['save_t1_mni']:

            t1_to_mni = npe.MapNode(seg_utils.ApplySegmentationDeformation(),
                                    name='t1_to_mni',
                                    iterfield=['deformation_field', 'in_files'])
            self.connect([
                (unzip_node, t1_to_mni, [('out_file', 'in_files')]),
                (new_segment, t1_to_mni, [('forward_deformation_field', 'deformation_field')]),
                (t1_to_mni, self.output_node, [('out_files', 't1_mni')])
            ])

        # Unzipping Templates
        # ============================
        unzip_templates_node = npe.Node(nutil.Function(input_names=['in_file'],
                                                       output_names=['out_file'],
                                                       function=unzip_nii),
                                        name='unzip_templates_node')

        # DARTEL with existing template
        # =============================
        dartel_existing_template = npe.MapNode(existing_dartel_utils.DARTELExistingTemplate(),
                                               name='dartel_existing_dartel',
                                               iterfield=['image_files'])
        if self.parameters['optimization_parameters'] is not None:
            dartel_existing_template.inputs.optimization_parameters = self.parameters['optimization_parameters']
        if self.parameters['regularization_form'] is not None:
            dartel_existing_template.inputs.regularization_form = self.parameters['regularization_form']

        # Unzipping Template
        # ============================
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
                                                                        function=dartel2mni_utils.join_smoothed_files),
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
                                                                   'atlas_list'],
                                                      output_names=['atlas_statistics'],
                                                      function=parcellation_utils.atlas_statistics),
                                       name='atlas_stats_node',
                                       iterfield=['in_image'])
        atlas_stats_node.inputs.atlas_list = self.parameters['atlas_list']

        # Connection
        # ==========
        self.connect([
            (self.input_node, unzip_node, [('input_images', 'in_file')]),
            (self.input_node, unzip_templates_node, [('dartel_iteration_templates', 'in_file')]),
            (self.input_node, unzip_template_node, [('dartel_final_template', 'in_file')]),
            (unzip_node, new_segment, [('out_file', 'channel_files')]),
            (new_segment, self.output_node, [('bias_corrected_images', 'bias_corrected_images'),
                                             ('bias_field_images', 'bias_field_images'),
                                             ('dartel_input_images', 'dartel_input_images'),
                                             ('forward_deformation_field', 'forward_deformation_field'),
                                             ('inverse_deformation_field', 'inverse_deformation_field'),
                                             ('modulated_class_images', 'modulated_class_images'),
                                             ('native_class_images', 'native_class_images'),
                                             ('normalized_class_images', 'normalized_class_images'),
                                             ('transformation_mat', 'transformation_mat')]),
            (new_segment, dartel_existing_template, [(('dartel_input_images',
                                                       existing_dartel_utils.prepare_images_from_segmentation),
                                                      'image_files')]),
            (unzip_templates_node, dartel_existing_template, [(('out_file',
                                                                existing_dartel_utils.create_iteration_parameters,
                                                                self.parameters['iteration_parameters']),
                                                               'iteration_parameters')]),
            (dartel_existing_template, self.output_node, [('dartel_flow_fields', 'dartel_flow_fields')]),
            (new_segment, dartel2mni_node, [(('native_class_images', seg_utils.group_nested_images_by_subject),
                                             'apply_to_files')]),
            (dartel_existing_template, dartel2mni_node, [(('dartel_flow_fields',
                                                           dartel2mni_utils.prepare_existing_dartel_flowfields,
                                                           self.parameters['tissue_classes']), 'flowfield_files')]),
            (unzip_template_node, dartel2mni_node, [('out_file', 'template_file')]),
            (dartel2mni_node, self.output_node, [('normalized_files', 'normalized_files')]),
            (dartel2mni_node, atlas_stats_node, [(('normalized_files', dartel2mni_utils.select_gm_images), 'in_image')]),
            (atlas_stats_node, self.output_node, [('atlas_statistics', 'atlas_statistics')])
        ])
