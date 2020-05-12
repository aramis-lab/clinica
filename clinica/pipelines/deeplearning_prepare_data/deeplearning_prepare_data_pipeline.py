# coding: utf8

import clinica.pipelines.engine as cpe


# Use hash instead of parameters for iterables folder names
# Otherwise path will be too long and generate OSError
from nipype import config
cfg = dict(execution={'parameterize_dirs': False})
config.update_config(cfg)


class Deeplearningpreparedata(cpe.Pipeline):
    """Deeplearning prepare data - MRI in nifty format are transformed into
    Pytorch tensors. The transformation is applied to: the whole volume, a
    selection of 3D patches, or slices extracted from the 3D volume. 

    Warnings:
        - A warning.

    Todos:
        - [x] A filled todo item.
        - [ ] An ongoing todo item.

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

        return ['t1w']  # Fill here the list

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """

        return []  # Fill here the list

    def build_input_node(self):
        """Build and connect an input node to the pipeline."""
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from clinica.utils.exceptions import ClinicaBIDSError, ClinicaException
        from clinica.utils.stream import cprint
        from clinica.iotools.utils.data_handling import check_volume_location_in_world_coordinate_system
        from clinica.utils.inputs import clinica_file_reader
        from clinica.utils.input_files import T1W_LINEAR
        from clinica.utils.ux import print_images_to_process
        from clinica.utils.filemanip import get_subject_id

        # Inputs from anat/ folder
        # ========================
        # T1w_Linear file:
        try:
            t1w_files = clinica_file_reader(self.subjects,
                                            self.sessions,
                                            self.caps_directory,
                                            T1W_LINEAR)
        except ClinicaException as e:
            err = 'Clinica faced error(s) while trying to read files in your BIDS directory.\n' + str(e)
            raise ClinicaBIDSError(err)

        if len(self.subjects):
            print_images_to_process(self.subjects, self.sessions)
            cprint('The pipeline will last approximately 30 seconds per image.')  # Replace by adequate computational time.

        # The reading node
        # -------------------------
        read_node = npe.Node(name="ReadingFiles",
                             iterables=[
                                 ('t1w', t1w_files),
                             ],
                             synchronize=True,
                             interface=nutil.IdentityInterface(
                                 fields=self.get_input_fields())
                             )
        
        # Get subject ID node
        # ----------------------
        image_id_node = npe.Node(
                interface=nutil.Function(
                    input_names=['bids_or_caps_file'],
                    output_names=['image_id'],
                    function=get_subject_id),
                name='ImageID'
                )
       
        self.connect([
            (read_node, self.input_node, [('t1w', 't1w')]),
        ])

    def build_output_node(self):
        """Build and connect an output node to the pipeline."""


    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline."""

        import deeplearning_prepare_data_utils as utils
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from .T1_preparedl_utils import (extract_slices,
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
                iterfield=['preprocessed_T1'],
                interface=nutil.Function(
                    function=extract_slices,
                    input_names=[
                        'preprocessed_T1', 'slice_direction',
                        'slice_mode'
                        ],
                    output_names=['output_file_rgb', 'output_file_original']
                    )
                )

        extract_slices.inputs.slice_direction = slice_direction
        extract_slices.inputs.slice_mode = slice_mode

        # Extract patches node (options, patch size and stride size)
        # ----------------------
        extract_patches = npe.MapNode(
                name='extract_patches',
                iterfield=['preprocessed_T1'],
                interface=nutil.Function(
                    function=extract_patches,
                    input_names=['preprocessed_T1', 'patch_size', 'stride_size'],
                    output_names=['output_patch']
                    )
                )

        extract_patches.inputs.patch_size = patch_size
        extract_patches.inputs.stride_size = stride_size
        # Step 1
        # ======
        node1 = npe.Node(name="Step1",
                         interface=nutil.Function(
                             input_names=['t1w', 'in_hello_word'],
                             output_names=[],
                             function=utils.step1))
        node1.inputs.in_hello_word = self.parameters['hello_word']

        # Step 2
        # ======
        node2 = npe.Node(name="Step2",
                         interface=nutil.Function(
                             input_names=['t1w', 'in_advanced_arg'],
                             output_names=[],
                             function=utils.step2))
        node2.inputs.in_advanced_arg = self.parameters['advanced_argument']

        # Connection
        # ==========
        self.connect([
            # STEP 1
            (self.input_node,      node1,    [('t1w',    't1w')]),
            # STEP 2
            (self.input_node,      node2,    [('t1w',    't1w')])
        ])
