# Use hash instead of parameters for iterables folder names
# Otherwise path will be too long and generate OSError
from nipype import config

import clinica.pipelines.engine as cpe

cfg = dict(execution={"parameterize_dirs": False})
config.update_config(cfg)


class PETLinear(cpe.PETPipeline):
    """PET Linear - Affine registration of PET images to standard space.
    This preprocessing pipeline uses T1w MRI transformation into standard
    space computed in clinica t1-linear pipeline.
    It includes globally four steps:
    1) PET image registration into associated T1w MRI image space
        with RegistrationSynQuick from ANTs.
    2) PET image registration to MNI152NLin2009cSym template with
       RegistrationSynQuick from ANTs.
    3) SUVR voxel intensity normalisation using Pons Cerebellum masks.
    4) Crop the background (in order to save computational power).

    Returns:
        A clinica pipeline object containing the pet_linear pipeline.
    """

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file."""
        pass

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Returns:
            A list of (string) input fields name.
        """
        return ["pet", "t1w", "t1w_to_mni"]

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """
        return [
            "registered_pet",
            "transform_mat",
            "registered_pet_in_t1w",
        ]  # Fill here the list

    def build_input_node(self):
        """Build and connect an input node to the pipeline."""
        from os import pardir
        from os.path import abspath, dirname, exists, join

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from clinica.utils.exceptions import (
            ClinicaBIDSError,
            ClinicaCAPSError,
            ClinicaException,
        )
        from clinica.utils.input_files import T1W_NII, T1W_TO_MNI_TRANSFORM
        from clinica.utils.inputs import (
            RemoteFileStructure,
            clinica_file_reader,
            fetch_file,
        )
        from clinica.utils.pet import get_suvr_mask
        from clinica.utils.stream import cprint
        from clinica.utils.ux import print_images_to_process

        # Import references files
        root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
        path_to_mask = join(root, "resources", "masks")
        url_aramis = "https://aramislab.paris.inria.fr/files/data/img_t1_linear/"
        FILE1 = RemoteFileStructure(
            filename="mni_icbm152_t1_tal_nlin_sym_09c.nii",
            url=url_aramis,
            checksum="93359ab97c1c027376397612a9b6c30e95406c15bf8695bd4a8efcb2064eaa34",
        )
        FILE2 = RemoteFileStructure(
            filename="ref_cropped_template.nii.gz",
            url=url_aramis,
            checksum="67e1e7861805a8fd35f7fcf2bdf9d2a39d7bcb2fd5a201016c4d2acdd715f5b3",
        )

        self.ref_template = join(path_to_mask, FILE1.filename)
        self.ref_crop = join(path_to_mask, FILE2.filename)
        self.ref_mask = get_suvr_mask(self.parameters["suvr_reference_region"])

        if not (exists(self.ref_template)):
            try:
                fetch_file(FILE1, path_to_mask)
            except IOError as err:
                cprint(
                    msg=f"Unable to download required template (mni_icbm152) for processing: {err}",
                    lvl="error",
                )
        if not (exists(self.ref_crop)):
            try:
                fetch_file(FILE2, path_to_mask)
            except IOError as err:
                cprint(
                    msg=f"Unable to download required template (ref_crop) for processing: {err}",
                    lvl="error",
                )

        # Inputs from BIDS directory
        try:
            pet_files, _ = clinica_file_reader(
                self.subjects,
                self.sessions,
                self.bids_directory,
                self._get_pet_scans_query(),
            )
        except ClinicaException as e:
            err = (
                "Clinica faced error(s) while trying to read pet files in your BIDS directory.\n"
                + str(e)
            )
            raise ClinicaBIDSError(err)

        # T1w file:
        try:
            t1w_files, _ = clinica_file_reader(
                self.subjects, self.sessions, self.bids_directory, T1W_NII
            )
        except ClinicaException as e:
            err = (
                "Clinica faced error(s) while trying to read t1w files in your BIDS directory.\n"
                + str(e)
            )
            raise ClinicaBIDSError(err)

        # Inputs from t1-linear pipeline
        # Transformation files from T1w files to MNI:
        try:
            t1w_to_mni_transformation_files, _ = clinica_file_reader(
                self.subjects, self.sessions, self.caps_directory, T1W_TO_MNI_TRANSFORM
            )
        except ClinicaException as e:
            err = (
                "Clinica faced error(s) while trying to read transformation files in your CAPS directory.\n"
                + str(e)
            )
            raise ClinicaCAPSError(err)

        if len(self.subjects):
            print_images_to_process(self.subjects, self.sessions)
            cprint("The pipeline will last approximately 3 minutes per image.")

        read_input_node = npe.Node(
            name="LoadingCLIArguments",
            iterables=[
                ("t1w", t1w_files),
                ("pet", pet_files),
                ("t1w_to_mni", t1w_to_mni_transformation_files),
            ],
            synchronize=True,
            interface=nutil.IdentityInterface(fields=self.get_input_fields()),
        )
        # fmt: off
        self.connect(
            [
                (read_input_node, self.input_node, [("t1w", "t1w")]),
                (read_input_node, self.input_node, [("pet", "pet")]),
                (read_input_node, self.input_node, [("t1w_to_mni", "t1w_to_mni")]),
            ]
        )
        # fmt: on

    def build_output_node(self):
        """Build and connect an output node to the pipeline."""

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from nipype.interfaces.io import DataSink

        from clinica.utils.nipype import container_from_filename, fix_join

        from .pet_linear_utils import rename_into_caps

        # Writing node
        write_node = npe.Node(name="WriteCaps", interface=DataSink())
        write_node.inputs.base_directory = self.caps_directory
        write_node.inputs.parameterization = False

        # Other nodes
        rename_file_node_inputs = [
            "pet_filename_bids",
            "pet_filename_raw",
            "transformation_filename_raw",
            "suvr_reference_region",
            "uncropped_image",
        ]
        if self.parameters.get("save_PETinT1w"):
            rename_file_node_inputs.append("pet_filename_in_t1w_raw")
        rename_files = npe.Node(
            interface=nutil.Function(
                input_names=rename_file_node_inputs,
                output_names=[
                    "pet_filename_caps",
                    "transformation_filename_caps",
                    "pet_filename_in_t1w_caps",
                ],
                function=rename_into_caps,
            ),
            name="renameFileCAPS",
        )
        rename_files.interface.inputs.suvr_reference_region = self.parameters.get(
            "suvr_reference_region"
        )
        rename_files.interface.inputs.uncropped_image = self.parameters.get(
            "uncropped_image"
        )
        container_path = npe.Node(
            interface=nutil.Function(
                input_names=["bids_or_caps_filename"],
                output_names=["container"],
                function=container_from_filename,
            ),
            name="containerPath",
        )
        # fmt: off
        self.connect(
            [
                (self.input_node, container_path, [("pet", "bids_or_caps_filename")]),
                (container_path, write_node, [(("container", fix_join, "pet_linear"), "container")]),
                (self.input_node, rename_files, [("pet", "pet_filename_bids")]),
                (self.output_node, rename_files, [("affine_mat", "transformation_filename_raw")]),
                (rename_files, write_node, [("transformation_filename_caps", "@transform_mat")]),
            ]
        )
        if not (self.parameters.get("uncropped_image")):
            self.connect(
                [
                    (self.output_node, rename_files, [("outfile_crop", "pet_filename_raw")]),
                    (rename_files, write_node, [("pet_filename_caps", "@registered_pet")]),
                ]
            )
        else:
            self.connect(
                [
                    (self.output_node, rename_files, [("suvr_pet", "pet_filename_raw")]),
                    (rename_files, write_node, [("pet_filename_caps", "@registered_pet")]),
                ]
            )
        if self.parameters.get("save_PETinT1w"):
            self.connect(
                [
                    (self.output_node, rename_files, [("PETinT1w", "pet_filename_in_t1w_raw")]),
                    (rename_files, write_node, [("pet_filename_in_t1w_caps", "@registered_pet_in_t1w")]),
                ]
            )
        # fmt: on

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline."""

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from nipype.interfaces import ants

        import clinica.pipelines.pet_linear.pet_linear_utils as utils

        # Utilitary nodes
        init_node = npe.Node(
            interface=nutil.Function(
                input_names=["pet"],
                output_names=["pet"],
                function=utils.init_input_node,
            ),
            name="initPipeline",
        )
        concatenate_node = npe.Node(
            interface=nutil.Function(
                input_names=["pet_to_t1w_transform", "t1w_to_mni_transform"],
                output_names=["transforms_list"],
                function=utils.concatenate_transforms,
            ),
            name="concatenateTransforms",
        )

        # The core (processing) nodes

        # 1. `RegistrationSynQuick` by *ANTS*. It uses nipype interface.
        ants_registration_node = npe.Node(
            name="antsRegistration", interface=ants.RegistrationSynQuick()
        )
        ants_registration_node.inputs.dimension = 3
        ants_registration_node.inputs.transform_type = "r"

        # 2. `ApplyTransforms` by *ANTS*. It uses nipype interface. PET to MRI
        ants_applytransform_node = npe.Node(
            name="antsApplyTransformPET2MNI", interface=ants.ApplyTransforms()
        )
        ants_applytransform_node.inputs.dimension = 3
        ants_applytransform_node.inputs.reference_image = self.ref_template

        # 3. Normalize the image (using nifti). It uses custom interface, from utils file
        ants_registration_nonlinear_node = npe.Node(
            name="antsRegistrationT1W2MNI", interface=ants.Registration()
        )
        ants_registration_nonlinear_node.inputs.fixed_image = self.ref_template
        ants_registration_nonlinear_node.inputs.metric = ["MI"]
        ants_registration_nonlinear_node.inputs.metric_weight = [1.0]
        ants_registration_nonlinear_node.inputs.transforms = ["SyN"]
        ants_registration_nonlinear_node.inputs.transform_parameters = [(0.1, 3, 0)]
        ants_registration_nonlinear_node.inputs.dimension = 3
        ants_registration_nonlinear_node.inputs.shrink_factors = [[8, 4, 2]]
        ants_registration_nonlinear_node.inputs.smoothing_sigmas = [[3, 2, 1]]
        ants_registration_nonlinear_node.inputs.sigma_units = ["vox"]
        ants_registration_nonlinear_node.inputs.number_of_iterations = [[200, 50, 10]]
        ants_registration_nonlinear_node.inputs.convergence_threshold = [1e-05]
        ants_registration_nonlinear_node.inputs.convergence_window_size = [10]
        ants_registration_nonlinear_node.inputs.radius_or_number_of_bins = [32]
        ants_registration_nonlinear_node.inputs.winsorize_lower_quantile = 0.005
        ants_registration_nonlinear_node.inputs.winsorize_upper_quantile = 0.995
        ants_registration_nonlinear_node.inputs.collapse_output_transforms = True
        ants_registration_nonlinear_node.inputs.use_histogram_matching = False
        ants_registration_nonlinear_node.inputs.verbose = True

        ants_applytransform_nonlinear_node = npe.Node(
            name="antsApplyTransformNonLinear", interface=ants.ApplyTransforms()
        )
        ants_applytransform_nonlinear_node.inputs.dimension = 3
        ants_applytransform_nonlinear_node.inputs.reference_image = self.ref_template

        if random_seed := self.parameters.get("random_seed", None):
            ants_registration_nonlinear_node.inputs.random_seed = random_seed

        normalize_intensity_node = npe.Node(
            name="intensityNormalization",
            interface=nutil.Function(
                function=utils.suvr_normalization,
                input_names=["input_img", "norm_img", "ref_mask"],
                output_names=["output_img"],
            ),
        )
        normalize_intensity_node.inputs.ref_mask = self.ref_mask

        # 4. Crop image (using nifti). It uses custom interface, from utils file
        crop_nifti_node = npe.Node(
            name="cropNifti",
            interface=nutil.Function(
                function=utils.crop_nifti,
                input_names=["input_img", "ref_img"],
                output_names=["output_img"],
            ),
        )
        crop_nifti_node.inputs.ref_img = self.ref_crop

        # 5. Print end message
        print_end_message = npe.Node(
            interface=nutil.Function(
                input_names=["pet", "final_file"], function=utils.print_end_pipeline
            ),
            name="WriteEndMessage",
        )

        # 6. Optional node: compute PET image in T1w
        ants_applytransform_optional_node = npe.Node(
            name="antsApplyTransformPET2T1w", interface=ants.ApplyTransforms()
        )
        ants_applytransform_optional_node.inputs.dimension = 3

        # Connection
        # ==========
        # fmt: off
        self.connect(
            [
                (self.input_node, init_node, [("pet", "pet")]),
                # STEP 1
                (self.input_node, ants_registration_node, [("t1w", "fixed_image")]),
                (init_node, ants_registration_node, [("pet", "moving_image")]),
                # STEP 2
                (ants_registration_node, concatenate_node, [("out_matrix", "pet_to_t1w_transform")]),
                (self.input_node, concatenate_node, [("t1w_to_mni", "t1w_to_mni_transform")]),
                (self.input_node, ants_applytransform_node, [("pet", "input_image")]),
                (concatenate_node, ants_applytransform_node, [("transforms_list", "transforms")]),
                # STEP 3
                (self.input_node, ants_registration_nonlinear_node, [("t1w", "moving_image")]),
                (ants_registration_nonlinear_node, ants_applytransform_nonlinear_node,
                 [("reverse_forward_transforms", "transforms")]),
                (ants_applytransform_node, ants_applytransform_nonlinear_node, [("output_image", "input_image")]),
                (ants_applytransform_node, normalize_intensity_node, [("output_image", "input_img")]),
                (ants_applytransform_nonlinear_node, normalize_intensity_node, [("output_image", "norm_img")]),
                # Connect to DataSink
                (ants_registration_node, self.output_node, [("out_matrix", "affine_mat")]),
                (normalize_intensity_node, self.output_node, [("output_img", "suvr_pet")]),
                (self.input_node, print_end_message, [("pet", "pet")]),
            ]
        )
        # STEP 4
        if not (self.parameters.get("uncropped_image")):
            self.connect(
                [
                    (normalize_intensity_node, crop_nifti_node, [("output_img", "input_img")]),
                    (crop_nifti_node, self.output_node, [("output_img", "outfile_crop")]),
                    (crop_nifti_node, print_end_message, [("output_img", "final_file")]),
                ]
            )
        else:
            self.connect(
                [
                    (normalize_intensity_node, print_end_message, [("output_img", "final_file")]),
                ]
            )
        # STEP 6: Optional argument
        if self.parameters.get("save_PETinT1w"):
            self.connect(
                [
                    (self.input_node, ants_applytransform_optional_node, [("pet", "input_image"),
                                                                          ("t1w", "reference_image")]),
                    (ants_registration_node, ants_applytransform_optional_node, [("out_matrix", "transforms")]),
                    (ants_applytransform_optional_node, self.output_node, [("output_image", "PETinT1w")]),
                ]
            )
        # fmt: on
