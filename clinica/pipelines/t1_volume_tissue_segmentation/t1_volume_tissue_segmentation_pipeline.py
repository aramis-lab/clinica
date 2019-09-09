# coding: utf8

import clinica.pipelines.engine as cpe

__author__ = "Jorge Samper-Gonzalez"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Jorge Samper-Gonzalez", "Alexandre Routier"]
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


class T1VolumeTissueSegmentation(cpe.Pipeline):
    """T1VolumeTissueSegmentation - Tissue segmentation, bias correction and
    spatial normalization to MNI space.

    Args:
        bids_directory: A BIDS directory.
        caps_directory: An empty output directory where CAPS structured data will be written.
        subjects_sessions_list: The Subjects-Sessions list file (in .tsv format).

    Returns:
        A clinica pipeline object containing the T1VolumeTissueSegmentation pipeline.
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
        """Build and connect an input node to the pipelines.
        """
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from colorama import Fore
        from clinica.utils.exceptions import ClinicaBIDSError, ClinicaException
        from clinica.utils.io import check_input_bids_files
        from clinica.utils.stream import cprint

        # Remove 'sub-' prefix from participant IDs
        participant_labels = '|'.join(sub[4:] for sub in self.subjects)
        # Remove 'ses-' prefix from session IDs
        session_labels = '|'.join(ses[4:] for ses in self.sessions)

        error_message = ""
        # Inputs from anat/ folder
        # ========================
        # T1w file:
        t1w_files = self.bids_layout.get(type='T1w', return_type='file', extensions=['.nii|.nii.gz'],
                                         subject=participant_labels, session=session_labels)
        error_message += check_input_bids_files(t1w_files, "T1W_NII",
                                                self.bids_directory, self.subjects, self.sessions)

        if error_message:
            raise ClinicaBIDSError(error_message)

        if len(t1w_files) == 0:
            import sys
            cprint(
                '%s\nEither all the images were already run by the pipeline or no image was found to run the pipeline. '
                'The program will now exit.%s' % (Fore.BLUE, Fore.RESET))
            sys.exit(0)
        else:
            images_to_process = ', '.join(self.subjects[i][4:] + '|' + self.sessions[i][4:]
                                          for i in range(len(self.subjects)))
            cprint('The pipeline will be run on the following subject(s): %s' % images_to_process)
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
        """Build and connect an output node to the pipelines.
        """
        pass

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipelines.
        """
        import nipype.pipeline.engine as npe
        import nipype.interfaces.utility as nutil
        import nipype.interfaces.io as nio
        import nipype.interfaces.spm as spm
        from ..t1_volume_tissue_segmentation import t1_volume_tissue_segmentation_utils as seg_utils
        from clinica.utils.io import unzip_nii
        from clinica.utils.io import (fix_join, zip_nii)
        from clinica.utils.spm import get_tpm

        # Get Tissue Probability Map from SPM
        tissue_map = get_tpm()

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

        new_segment.inputs.tissues = seg_utils.get_tissue_tuples(
            tissue_map,
            self.parameters['tissue_classes'],
            self.parameters['dartel_tissues'],
            self.parameters['save_warped_unmodulated'],
            self.parameters['save_warped_modulated'])

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

        # Apply segmentation deformation to T1 (into MNI space)
        # =====================================================
        if self.parameters['save_t1_mni'] is not None and self.parameters['save_t1_mni']:

            t1_to_mni = npe.Node(seg_utils.ApplySegmentationDeformation(),
                                 name='3-T1wToMni')
            self.connect([
                (unzip_node, t1_to_mni, [('out_file', 'in_files')]),
                (new_segment, t1_to_mni, [('forward_deformation_field', 'deformation_field')]),
                (t1_to_mni, self.output_node, [('out_files', 't1_mni')]),
                (self.output_node, print_end_message, [('t1_mni', 'final_file')]),
            ])
        else:
            self.connect([
                (self.output_node, print_end_message, [('transformation_mat', 'final_file')]),
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
            (r'(.*)c1(sub-.*)(\.nii(\.gz)?)$',                                 r'\1\2_segm-graymatter\3'),  # noqa
            (r'(.*)c2(sub-.*)(\.nii(\.gz)?)$',                                 r'\1\2_segm-whitematter\3'),  # noqa
            (r'(.*)c3(sub-.*)(\.nii(\.gz)?)$',                                 r'\1\2_segm-csf\3'),  # noqa
            (r'(.*)c4(sub-.*)(\.nii(\.gz)?)$',                                 r'\1\2_segm-bone\3'),  # noqa
            (r'(.*)c5(sub-.*)(\.nii(\.gz)?)$',                                 r'\1\2_segm-softtissue\3'),  # noqa
            (r'(.*)c6(sub-.*)(\.nii(\.gz)?)$',                                 r'\1\2_segm-background\3'),  # noqa
            (r'(.*)(/native_space/sub-.*)(\.nii(\.gz)?)$',                     r'\1\2_probability\3'),  # noqa
            (r'(.*)(/([a-z]+)_deformation_field/)i?y_(sub-.*)(\.nii(\.gz)?)$', r'\1/normalized_space/\4_target-Ixi549Space_transformation-\3_deformation\5'),  # noqa
            (r'(.*)(/t1_mni/)w(sub-.*)_T1w(\.nii(\.gz)?)$',                    r'\1/normalized_space/\3_space-Ixi549Space_T1w\4'),  # noqa
            (r'(.*)(/modulated_normalized/)mw(sub-.*)(\.nii(\.gz)?)$',         r'\1/normalized_space/\3_space-Ixi549Space_modulated-on_probability\4'),  # noqa
            (r'(.*)(/normalized/)w(sub-.*)(\.nii(\.gz)?)$',                    r'\1/normalized_space/\3_space-Ixi549Space_modulated-off_probability\4'),  # noqa
            (r'(.*/dartel_input/)r(sub-.*)(\.nii(\.gz)?)$',                    r'\1\2_dartelinput\3'),  # noqa
            # Will remove trait_added empty folder
            (r'trait_added', r'')
        ]

        self.connect([
            (self.input_node, container_path, [('t1w', 't1w_filename')]),  # noqa
            (container_path, write_node, [(('container', fix_join, ''), 'container')]),  # noqa
            (self.output_node, write_node, [(('native_class_images', seg_utils.zip_list_files, True), 'native_space'),  # noqa
                                            (('dartel_input_images', seg_utils.zip_list_files, True), 'dartel_input')]),  # noqa
        ])
        if self.parameters['save_warped_unmodulated']:
            self.connect([
                (self.output_node, write_node, [(('normalized_class_images', seg_utils.zip_list_files, True), 'normalized')]),  # noqa
            ])
        if self.parameters['save_warped_modulated']:
            self.connect([
                (self.output_node, write_node, [(('modulated_class_images', seg_utils.zip_list_files, True), 'modulated_normalized')]),  # noqa
            ])
        if self.parameters['write_deformation_fields'] is not None:
            if self.parameters['write_deformation_fields'][0]:
                self.connect([
                    (self.output_node, write_node, [(('inverse_deformation_field', zip_nii, True), 'inverse_deformation_field')]),  # noqa
                ])
            if self.parameters['write_deformation_fields'][1]:
                self.connect([
                    (self.output_node, write_node, [(('forward_deformation_field', zip_nii, True), 'forward_deformation_field')]),  # noqa
                ])
        if self.parameters['save_t1_mni'] is not None and self.parameters['save_t1_mni']:
            self.connect([
                (self.output_node, write_node, [(('t1_mni', zip_nii, True), 't1_mni')]),  # noqa
            ])
