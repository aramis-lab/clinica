# coding: utf8

import clinica.pipelines.engine as cpe

# Use hash instead of parameters for iterables folder names
# Otherwise path will be too long and generate OSError
from nipype import config
cfg = dict(execution={'parameterize_dirs': False})
config.update_config(cfg)


class T1VolumeTissueSegmentation(cpe.Pipeline):
    """T1VolumeTissueSegmentation - Tissue segmentation, bias correction and
    spatial normalization to MNI space.

    Returns:
        A clinica pipeline object containing the T1VolumeTissueSegmentation pipeline.
    """
    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file."""
        pass

    def check_pipeline_parameters(self):
        """Check pipeline parameters."""
        from clinica.utils.spm import get_tpm

        if 'tissue_classes' not in self.parameters.keys():
            self.parameters['tissue_classes'] = [1, 2, 3]
        if 'dartel_tissues' not in self.parameters.keys():
            self.parameters['dartel_tissues'] = [1, 2, 3]
        if 'tissue_probability_maps' not in self.parameters.keys():
            self.parameters['tissue_probability_maps'] = None
        if 'save_warped_unmodulated' not in self.parameters.keys():
            self.parameters['save_warped_unmodulated'] = True
        if 'save_warped_modulated' not in self.parameters.keys():
            self.parameters['save_warped_modulated'] = False

        # Get Tissue Probability Map from SPM
        if self.parameters['tissue_probability_maps'] is None:
            self.parameters['tissue_probability_maps'] = get_tpm()

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipelines.

        Returns:
            A list of (string) input fields name.
        """
        return ['t1w']

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
        """Build and connect an input node to the pipeline.

        Raise:
            ClinicaBIDSError: If there are duplicated files or missing files for any subject
        """
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from clinica.utils.exceptions import ClinicaBIDSError, ClinicaException
        from clinica.utils.stream import cprint
        from clinica.iotools.utils.data_handling import check_volume_location_in_world_coordinate_system
        from clinica.utils.inputs import clinica_file_reader
        from clinica.utils.input_files import T1W_NII
        from clinica.utils.ux import print_images_to_process

        # Inputs from anat/ folder
        # ========================
        # T1w file:
        try:
            t1w_files = clinica_file_reader(self.subjects,
                                            self.sessions,
                                            self.bids_directory,
                                            T1W_NII)
        except ClinicaException as e:
            err = 'Clinica faced error(s) while trying to read files in your BIDS directory.\n' + str(e)
            raise ClinicaBIDSError(err)

        check_volume_location_in_world_coordinate_system(t1w_files, self.bids_directory)

        if len(self.subjects):
            print_images_to_process(self.subjects, self.sessions)
            cprint('The pipeline will last approximately 10 minutes per image.')

        read_node = npe.Node(name="ReadingFiles",
                             iterables=[
                                 ('t1w', t1w_files),
                             ],
                             synchronize=True,
                             interface=nutil.IdentityInterface(
                                 fields=self.get_input_fields())
                             )
        self.connect([
            (read_node, self.input_node, [('t1w', 't1w')]),
        ])

    def build_output_node(self):
        """Build and connect an output node to the pipeline."""
        pass

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline."""
        import nipype.pipeline.engine as npe
        import nipype.interfaces.utility as nutil
        import nipype.interfaces.io as nio
        import nipype.interfaces.spm as spm
        from ..t1_volume_tissue_segmentation import t1_volume_tissue_segmentation_utils as seg_utils
        from clinica.utils.filemanip import unzip_nii, zip_nii
        from clinica.utils.nipype import fix_join

        # Get <subject_id> (e.g. sub-CLNC01_ses-M00) from input_node
        # and print begin message
        # =======================
        init_node = npe.Node(interface=nutil.Function(
            input_names=self.get_input_fields(),
            output_names=['subject_id'] + self.get_input_fields(),
            function=seg_utils.init_input_node),
            name='0-InitNode')

        # Unzipping
        # =========
        unzip_node = npe.Node(nutil.Function(input_names=['in_file'],
                                             output_names=['out_file'],
                                             function=unzip_nii),
                              name='1-UnzipT1w')

        # Unified Segmentation
        # ====================
        new_segment = npe.Node(spm.NewSegment(),
                               name='2-SpmSegmentation')
        new_segment.inputs.write_deformation_fields = [True, True]
        new_segment.inputs.tissues = seg_utils.get_tissue_tuples(
            self.parameters['tissue_probability_maps'],
            self.parameters['tissue_classes'],
            self.parameters['dartel_tissues'],
            self.parameters['save_warped_unmodulated'],
            self.parameters['save_warped_modulated']
        )

        # Apply segmentation deformation to T1 (into MNI space)
        # =====================================================
        t1_to_mni = npe.Node(seg_utils.ApplySegmentationDeformation(),
                             name='3-T1wToMni')

        # Print end message
        # =================
        print_end_message = npe.Node(
            interface=nutil.Function(
                input_names=['subject_id', 'final_file'],
                function=seg_utils.print_end_pipeline),
            name='WriteEndMessage')

        # Connection
        # ==========
        self.connect([
            (self.input_node, init_node, [('t1w', 't1w')]),
            (init_node, unzip_node, [('t1w', 'in_file')]),
            (unzip_node, new_segment, [('out_file', 'channel_files')]),
            (init_node, print_end_message, [('subject_id', 'subject_id')]),
            (unzip_node, t1_to_mni, [('out_file', 'in_files')]),
            (new_segment, t1_to_mni, [('forward_deformation_field', 'deformation_field')]),
            (new_segment, self.output_node, [('bias_corrected_images', 'bias_corrected_images'),
                                             ('bias_field_images', 'bias_field_images'),
                                             ('dartel_input_images', 'dartel_input_images'),
                                             ('forward_deformation_field', 'forward_deformation_field'),
                                             ('inverse_deformation_field', 'inverse_deformation_field'),
                                             ('modulated_class_images', 'modulated_class_images'),
                                             ('native_class_images', 'native_class_images'),
                                             ('normalized_class_images', 'normalized_class_images'),
                                             ('transformation_mat', 'transformation_mat')]),
            (t1_to_mni, self.output_node, [('out_files', 't1_mni')]),
            (self.output_node, print_end_message, [('t1_mni', 'final_file')]),
        ])

        # Find container path from t1w filename
        # =====================================
        container_path = npe.Node(
            nutil.Function(input_names=['t1w_filename'],
                           output_names=['container'],
                           function=seg_utils.t1w_container_from_filename),
            name='ContainerPath')

        # Writing CAPS
        # ============
        write_node = npe.Node(name='WriteCAPS',
                              interface=nio.DataSink()
                              )
        write_node.inputs.base_directory = self.caps_directory
        write_node.inputs.parameterization = False
        write_node.inputs.regexp_substitutions = [
            (r'(.*)c1(sub-.*)(\.nii(\.gz)?)$',
             r'\1\2_segm-graymatter\3'),
            (r'(.*)c2(sub-.*)(\.nii(\.gz)?)$',
             r'\1\2_segm-whitematter\3'),
            (r'(.*)c3(sub-.*)(\.nii(\.gz)?)$',
             r'\1\2_segm-csf\3'),
            (r'(.*)c4(sub-.*)(\.nii(\.gz)?)$',
             r'\1\2_segm-bone\3'),
            (r'(.*)c5(sub-.*)(\.nii(\.gz)?)$',
             r'\1\2_segm-softtissue\3'),
            (r'(.*)c6(sub-.*)(\.nii(\.gz)?)$',
             r'\1\2_segm-background\3'),
            (r'(.*)(/native_space/sub-.*)(\.nii(\.gz)?)$',
             r'\1\2_probability\3'),
            (r'(.*)(/([a-z]+)_deformation_field/)i?y_(sub-.*)(\.nii(\.gz)?)$',
             r'\1/normalized_space/\4_target-Ixi549Space_transformation-\3_deformation\5'),
            (r'(.*)(/t1_mni/)w(sub-.*)_T1w(\.nii(\.gz)?)$', r'\1/normalized_space/\3_space-Ixi549Space_T1w\4'),
            (r'(.*)(/modulated_normalized/)mw(sub-.*)(\.nii(\.gz)?)$',
             r'\1/normalized_space/\3_space-Ixi549Space_modulated-on_probability\4'),
            (r'(.*)(/normalized/)w(sub-.*)(\.nii(\.gz)?)$',
             r'\1/normalized_space/\3_space-Ixi549Space_modulated-off_probability\4'),
            (r'(.*/dartel_input/)r(sub-.*)(\.nii(\.gz)?)$',
             r'\1\2_dartelinput\3'),
            # Will remove trait_added empty folder
            (r'trait_added', r'')
        ]

        self.connect([
            (self.input_node, container_path, [('t1w', 't1w_filename')]),
            (container_path, write_node, [(('container', fix_join, ''), 'container')]),
            (self.output_node, write_node, [(('native_class_images', seg_utils.zip_list_files, True), 'native_space'),
                                            (('dartel_input_images', seg_utils.zip_list_files, True), 'dartel_input')]),
            (self.output_node, write_node, [
                (('inverse_deformation_field', zip_nii, True), 'inverse_deformation_field')
            ]),
            (self.output_node, write_node, [
                (('forward_deformation_field', zip_nii, True), 'forward_deformation_field')
            ]),
            (self.output_node, write_node, [(('t1_mni', zip_nii, True), 't1_mni')]),
        ])
        if self.parameters['save_warped_unmodulated']:
            self.connect([
                (self.output_node, write_node, [
                    (('normalized_class_images', seg_utils.zip_list_files, True), 'normalized')
                ]),
            ])
        if self.parameters['save_warped_modulated']:
            self.connect([
                (self.output_node, write_node, [
                    (('modulated_class_images', seg_utils.zip_list_files, True), 'modulated_normalized')
                ]),
            ])
