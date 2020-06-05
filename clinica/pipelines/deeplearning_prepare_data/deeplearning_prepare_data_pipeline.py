# -*- coding: utf-8 -*-


import clinica.pipelines.engine as cpe


# Use hash instead of parameters for iterables folder names
# Otherwise path will be too long and generate OSError
from nipype import config
cfg = dict(execution={'parameterize_dirs': False})
config.update_config(cfg)


class DeepLearningPrepareData(cpe.Pipeline):
    """Deeplearning prepare data - MRI in nifti format are transformed into
    PyTorch tensors. The transformation is applied to: the whole volume, a
    selection of 3D patches, or slices extracted from the 3D volume. By default
    it uses the cropped version of the MRI (see option "--use_uncropper_image")


    Returns:
        A clinica pipeline object containing the Deeplearning prepare data pipeline.

    Raises:

    """

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file."""
        pass

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Returns:
            A list of (string) input fields name.
        """

        return ['input_nifti']

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """

        return ['image_id']  # Fill here the list

    def build_input_node(self):
        """Build and connect an input node to the pipeline."""
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from clinica.utils.exceptions import ClinicaBIDSError, ClinicaException
        from clinica.utils.stream import cprint
        from clinica.utils.inputs import clinica_file_reader
        from clinica.utils.input_files import T1W_LINEAR
        from clinica.utils.input_files import T1W_LINEAR_CROPPED
        from clinica.utils.ux import print_images_to_process

        if self.parameters.get('use_uncropped_image'):
            FILE_TYPE = T1W_LINEAR
        else:
            FILE_TYPE = T1W_LINEAR_CROPPED

        # T1w_Linear file:
        try:
            t1w_files = clinica_file_reader(self.subjects,
                                            self.sessions,
                                            self.caps_directory,
                                            FILE_TYPE)
        except ClinicaException as e:
            err = 'Clinica faced error(s) while trying to read files in your CAPS directory.\n' + str(e)
            raise ClinicaBIDSError(err)

        if len(self.subjects):
            print_images_to_process(self.subjects, self.sessions)
            cprint('The pipeline will last approximately 30 seconds per image.')  # Replace by adequate computational time.

        if self.parameters.get('extract_method') == 'slice':
            self.slice_direction = self.parameters.get('slice_direction')
            self.slice_mode = self.parameters.get('slice_mode')
        else:
            self.slice_direction = 'axial'
            self.slice_mode = 'rgb'

        if self.parameters.get('extract_method') == 'patch':
            self.patch_size = self.parameters.get('patch_size')
            self.stride_size = self.parameters.get('stride_size')
        else:
            self.patch_size = 50
            self.stride_size = 50

        # The reading node
        # -------------------------
        read_node = npe.Node(name="ReadingFiles",
                             iterables=[
                                 ('input_nifti', t1w_files),
                             ],
                             synchronize=True,
                             interface=nutil.IdentityInterface(
                                 fields=self.get_input_fields())
                             )

        self.connect([
            (read_node, self.input_node, [('input_nifti', 'input_nifti')]),
        ])

    def build_output_node(self):
        """Build and connect an output node to the pipeline."""
        import nipype.interfaces.utility as nutil
        from nipype.interfaces.io import DataSink
        import nipype.pipeline.engine as npe
        from clinica.utils.nipype import (fix_join, container_from_filename)
        from clinica.utils.filemanip import get_subject_id

        # Write node
        # ----------------------
        write_node = npe.Node(
                name="WriteCaps",
                interface=DataSink()
                )
        write_node.inputs.base_directory = self.caps_directory
        write_node.inputs.parameterization = False

        # Get subject ID node
        # ----------------------
        image_id_node = npe.Node(
                interface=nutil.Function(
                    input_names=['bids_or_caps_file'],
                    output_names=['image_id'],
                    function=get_subject_id),
                name='ImageID'
                )

        # Find container path from t1w filename
        # ----------------------
        container_path = npe.Node(
                nutil.Function(
                    input_names=['bids_or_caps_filename'],
                    output_names=['container'],
                    function=container_from_filename),
                name='ContainerPath')

        self.connect([
            (self.input_node, image_id_node, [('input_nifti', 'bids_or_caps_file')]),
            (self.input_node, container_path, [('input_nifti', 'bids_or_caps_filename')]),
            # (image_id_node, write_node, [('image_id', '@image_id')]),
            (image_id_node, write_node, [('image_id', '@image_id')]),
            ])

        subfolder = 'image_based'
        if self.parameters.get('extract_method') == 'slice':
            subfolder = 'slice_based'
            self.connect([
                (self.output_node, write_node, [('slices_rgb_T1', '@slices_rgb_T1')]),
                (self.output_node, write_node, [('slices_original_T1', '@slices_original_T1')])
                ])

        elif self.parameters.get('extract_method') == 'patch':
            subfolder = 'patch_based'
            self.connect([
                (self.output_node, write_node, [('patches_T1', '@patches_T1')])
                ])
        else:
            self.connect([
                (self.output_node, write_node, [('output_pt_file', '@output_pt_file')])
                ])

        self.connect([
            (container_path, write_node, [(
                (
                    'container', fix_join,
                    'deeplearning_prepare_data', subfolder, 't1_linear'
                    ),
                'container')]),
            ])

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline."""

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from .deeplearning_prepare_data_utils import (extract_slices,
                                                      extract_patches,
                                                      save_as_pt)
        # The processing nodes

        # Node to save MRI in nii.gz format into pytorch .pt format
        # ----------------------
        save_as_pt = npe.MapNode(
               name='save_as_pt',
               iterfield=['input_img'],
               interface=nutil.Function(
                   function=save_as_pt,
                   input_names=['input_img'],
                   output_names=['output_file']
                   )
               )

        # Extract slices node (options: 3 directions, mode)
        # ----------------------
        extract_slices = npe.MapNode(
                name='extract_slices',
                iterfield=['input_tensor'],
                interface=nutil.Function(
                    function=extract_slices,
                    input_names=[
                        'input_tensor', 'slice_direction',
                        'slice_mode'
                        ],
                    output_names=['output_file_rgb', 'output_file_original']
                    )
                )

        extract_slices.inputs.slice_direction = self.slice_direction
        extract_slices.inputs.slice_mode = self.slice_mode

        # Extract patches node (options, patch size and stride size)
        # ----------------------
        extract_patches = npe.MapNode(
                name='extract_patches',
                iterfield=['input_tensor'],
                interface=nutil.Function(
                    function=extract_patches,
                    input_names=['input_tensor', 'patch_size', 'stride_size'],
                    output_names=['output_patch']
                    )
                )

        extract_patches.inputs.patch_size = self.patch_size
        extract_patches.inputs.stride_size = self.stride_size

        # Connections
        # ----------------------
        self.connect([
            (self.input_node, save_as_pt, [('input_nifti', 'input_img')]),
            ])

        if self.parameters.get('extract_method') == 'slice':
            self.connect([
                (save_as_pt, extract_slices, [('output_file', 'input_tensor')]),
                (extract_slices, self.output_node, [('output_file_rgb', 'slices_rgb_T1')]),
                (extract_slices, self.output_node, [('output_file_original', 'slices_original_T1')])
                ])
        elif self.parameters.get('extract_method') == 'patch':
            self.connect([
                (save_as_pt, extract_patches, [('output_file', 'input_tensor')]),
                (extract_patches, self.output_node, [('output_patch', 'patches_T1')])
                ])
        else:
            self.connect([
                (save_as_pt, self.output_node, [('output_file', 'output_pt_file')]),
                ])
