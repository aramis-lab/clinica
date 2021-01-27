# coding: utf8

"""pet_linear - Clinica Pipeline.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details: http://www.clinica.run/doc/InteractingWithClinica/.
"""

# WARNING: Don't put any import statement here except if it's absolutly
# necessary. Put it *inside* the different methods.
# Otherwise it will slow down the dynamic loading of the pipelines list by the
# command line tool.
import clinica.pipelines.engine as cpe


# Use hash instead of parameters for iterables folder names
# Otherwise path will be too long and generate OSError
from nipype import config
cfg = dict(execution={'parameterize_dirs': False})
config.update_config(cfg)


class pet_linear(cpe.Pipeline):
    """pet_linear SHORT DESCRIPTION.

    Warnings:
        - A warning.

    Todos:
        - [x] A filled todo item.
        - [ ] An ongoing todo item.

    Returns:
        A clinica pipeline object containing the pet_linear pipeline.

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
        from clinica.utils.input_files import T1W_NII
        from clinica.utils.ux import print_images_to_process

        # This node is supposedly used to load BIDS and/or CAPS inputs when this pipeline is
        # not already connected to the output of a previous Clinica pipeline.
        # For this example, we read T1w MRI data which are passed to a read_node with iterable.
        # This allows to parallelize the pipelines accross sessions
        # when connected to the `self.input_node`.

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

        if len(self.subjects):
            print_images_to_process(self.subjects, self.sessions)
            # Replace by adequate computational time in the line below.
            cprint('The pipeline will last approximately 42 minutes per image.')

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

        # In the same idea as the input node, this output node is supposedly
        # used to write the output fields in a CAPS. It should be executed only
        # if this pipeline output is not already connected to a next Clinica
        # pipeline.

        pass

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline."""

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from clinica.utils.filemanip import get_filename_no_ext
        from nipype.interfaces import ants
        from .pet_linear_utils import suvr_normalization, crop_nifti, print_end_pipeline

        image_id_node = npe.Node(
                interface=nutil.Function(
                    input_names=['filename'],
                    output_names=['image_id'],
                    function=get_filename_no_ext),
                name='ImageID'
                )

        # The core (processing) nodes
        # =====================================

        # 1. `RegistrationSynQuick` by *ANTS*. It uses nipype interface.
        ants_registration_node = npe.Node(
            name='antsRegistration',
            interface=ants.RegistrationSynQuick()
        )
        # Fixed image is the MRI image and moving image is pet
        ants_registration_node.inputs.dimension = 3
        ants_registration_node.inputs.transform_type = 'r'

        # 2. `ApplyTransforms` by *ANTS*. It uses nipype interface. PET to MRI
        ants_applytransform1_node = npe.Node(
            name='antsApplyTransform1',
            interface=ants.ApplyTransforms()
        )
        # Reference image is the MRI image and input image is pet
        ants_applytransform1_node.inputs.dimension = 3

        # 3. `ApplyTransforms` by *ANTS*. It uses nipype interface. PET to MNI
        ants_applytransform2_node = npe.Node(
            name='antsApplyTransform2',
            interface=ants.ApplyTransforms()
        )
        # Input image is pet output from previous node
        ants_applytransform1_node.inputs.dimension = 3
        ants_applytransform2_node.inputs.reference_image = self.ref_template

        # 4. Normalize the image (using nifti). It uses custom interface, from utils file
        normalize_intensity_node = npe.Node(
            name='intensityNormalization',
            interface=nutil.Function(
                function=suvr_normalization,
                input_names=['input_img', 'ref_mask'],
                output_names=['output_img', 'mask_template']
            )
        )
        normalize_intensity_node.inputs.ref_mask = self.ref_mask

        # 5. Crop image (using nifti). It uses custom interface, from utils file
        crop_nifti_node = npe.Node(
                name='cropnifti',
                interface=nutil.Function(
                    function=crop_nifti,
                    input_names=['input_img', 'ref_crop'],
                    output_names=['output_img', 'crop_template']
                    )
                )
        crop_nifti_node.inputs.ref_crop = self.ref_crop

        # 4. Print end message
        print_end_message = npe.Node(
            interface=nutil.Function(
                input_names=['pet', 'final_file'],
                function=print_end_pipeline),
            name='WriteEndMessage')

        # Connection
        # ==========
        self.connect([
            (self.input_node, image_id_node, [('t1w', 'filename')]),
            # STEP 1
            (self.input_node, ants_registration_node, [('t1w', 'fixed_image')]),
            (self.input_node, ants_registration_node, [('pet', 'moving_image')]),
            (ants_registration_node, ants_applytransform1_node, [('out_matrix', 'transforms')]),
            (self.input_node, ants_applytransform1_node, [('pet', 'input_node')]),
            (self.input_node, ants_applytransform1_node, [('t1w', 'reference_image')]),
            # STEP 2
            (self.input_node, ants_applytransform2_node, [('mni_trans', 'transforms')]),
            (ants_applytransform1_node, ants_applytransform2_node, [('output_image', 'input_node')]),
            # STEP 3
            (ants_applytransform2_node, normalize_intensity_node, [('output_image', 'input_img')])

            # Connect to DataSink
            (image_id_node, self.output_node, [('image_id', 'image_id')]),
            (ants_registration_node, self.output_node, [('out_matrix', 'affine_mat')]),  # to change
            (ants_applytransform1_node, self.output_node, [('output_image', 'pet_in_mri')]),
            (ants_applytransform2_node, self.output_node, [('output_image', 'pet_in_mni')]),
            (normalize_intensity_node, self.output_node, [('output_image', 'suvr_pet')]),
            (normalize_intensity_node, self.output_node, [('mask_template', 'outfile_mask')]),
            (self.input_node, print_end_message, [('pet', 'pet')]),
        ])
        # STEP 4
        if not (self.parameters.get('uncropped_image')):
            self.connect([
                (normalize_intensity_node, crop_nifti_node, [('output_img', 'input_img')])
                (cropnifti, self.output_node, [('output_img', 'outfile_crop')]),
                (cropnifti, print_end_message, [('output_img', 'final_file')]),
                ])
        else:
            self.connect([
                (normalize_intensity_node, print_end_message, [('output_img', 'final_file')]),
            ])
