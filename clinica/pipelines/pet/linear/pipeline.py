# Use hash instead of parameters for iterables folder names
# Otherwise path will be too long and generate OSError
from typing import List, Optional

import nipype.pipeline.engine as npe
from nipype import config

from clinica.pipelines.pet.engine import PETPipeline

cfg = dict(execution={"parameterize_dirs": False})
config.update_config(cfg)


class PETLinear(PETPipeline):
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

    def __init__(
        self,
        bids_directory: Optional[str] = None,
        caps_directory: Optional[str] = None,
        tsv_file: Optional[str] = None,
        overwrite_caps: Optional[bool] = False,
        base_dir: Optional[str] = None,
        parameters: Optional[dict] = None,
        name: Optional[str] = None,
        ignore_dependencies: Optional[List[str]] = None,
        use_antspy: bool = False,
    ):
        from clinica.utils.stream import cprint

        super().__init__(
            bids_directory=bids_directory,
            caps_directory=caps_directory,
            tsv_file=tsv_file,
            overwrite_caps=overwrite_caps,
            base_dir=base_dir,
            parameters=parameters,
            ignore_dependencies=ignore_dependencies,
            name=name,
        )
        self.use_antspy = use_antspy
        if self.use_antspy:
            self._ignore_dependencies.append("ants")
            cprint(
                (
                    "The PETLinear pipeline has been configured to use ANTsPy instead of ANTs.\n"
                    "This means that no installation of ANTs is required, but the antspyx Python "
                    "package must be installed in your environment.\nThis functionality has been "
                    "introduced in Clinica 0.9.0 and is considered experimental.\n"
                    "Please report any issue or unexpected results to the Clinica developer team."
                ),
                lvl="warning",
            )

    def _check_custom_dependencies(self) -> None:
        """Check dependencies that can not be listed in the `info.json` file."""
        pass

    def get_input_fields(self) -> List[str]:
        """Specify the list of possible inputs of this pipeline.

        Returns
        -------
        list of str :
            A list of (string) input fields name.
        """
        return ["pet", "t1w", "t1w_to_mni", "t1w_linear"]

    def get_output_fields(self) -> List[str]:
        """Specify the list of possible outputs of this pipeline.

        Returns
        -------
        list of str :
            A list of (string) output fields name.
        """
        return [
            "registered_pet",
            "transform_mat",
            "registered_pet_in_t1w",
        ]

    def _build_input_node(self):
        """Build and connect an input node to the pipeline."""
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from clinica.pipelines.pet.utils import get_suvr_mask
        from clinica.utils.exceptions import (
            ClinicaBIDSError,
            ClinicaCAPSError,
            ClinicaException,
        )
        from clinica.utils.image import get_mni_template
        from clinica.utils.input_files import (
            T1W_LINEAR,
            T1W_LINEAR_CROPPED,
            T1W_NII,
            T1W_TO_MNI_TRANSFORM,
        )
        from clinica.utils.inputs import (
            clinica_file_reader,
            format_clinica_file_reader_errors,
        )
        from clinica.utils.stream import cprint
        from clinica.utils.ux import print_images_to_process

        self.ref_template = get_mni_template("t1")
        self.ref_mask = get_suvr_mask(self.parameters["suvr_reference_region"])

        # Inputs from BIDS directory
        pet_files, pet_errors = clinica_file_reader(
            self.subjects,
            self.sessions,
            self.bids_directory,
            self._get_pet_scans_query(),
        )
        if pet_errors:
            raise ClinicaBIDSError(
                format_clinica_file_reader_errors(
                    pet_errors, self._get_pet_scans_query()
                )
            )

        # T1w file:
        t1w_files, t1w_errors = clinica_file_reader(
            self.subjects, self.sessions, self.bids_directory, T1W_NII
        )
        if t1w_errors:
            raise ClinicaBIDSError(
                format_clinica_file_reader_errors(t1w_errors, T1W_NII)
            )

        # Inputs from t1-linear pipeline
        # T1w images registered
        t1w_linear_file_pattern = (
            T1W_LINEAR
            if self.parameters.get("uncropped_image", False)
            else T1W_LINEAR_CROPPED
        )
        t1w_linear_files, t1w_linear_errors = clinica_file_reader(
            self.subjects, self.sessions, self.caps_directory, t1w_linear_file_pattern
        )
        if t1w_linear_errors:
            raise ClinicaCAPSError(
                format_clinica_file_reader_errors(
                    t1w_linear_errors, t1w_linear_file_pattern
                )
            )
        # Transformation files from T1w files to MNI:
        t1w_to_mni_transformation_files, t1w_to_mni_errors = clinica_file_reader(
            self.subjects, self.sessions, self.caps_directory, T1W_TO_MNI_TRANSFORM
        )
        if t1w_to_mni_errors:
            raise ClinicaCAPSError(
                format_clinica_file_reader_errors(
                    t1w_to_mni_errors, T1W_TO_MNI_TRANSFORM
                )
            )

        if len(self.subjects):
            print_images_to_process(self.subjects, self.sessions)
            cprint("The pipeline will last approximately 3 minutes per image.")

        read_input_node = npe.Node(
            name="LoadingCLIArguments",
            iterables=[
                ("t1w", t1w_files),
                ("pet", pet_files),
                ("t1w_to_mni", t1w_to_mni_transformation_files),
                ("t1w_linear", t1w_linear_files),
            ],
            synchronize=True,
            interface=nutil.IdentityInterface(fields=self.get_input_fields()),
        )
        self.connect(
            [
                (read_input_node, self.input_node, [("t1w", "t1w")]),
                (read_input_node, self.input_node, [("pet", "pet")]),
                (read_input_node, self.input_node, [("t1w_to_mni", "t1w_to_mni")]),
                (read_input_node, self.input_node, [("t1w_linear", "t1w_linear")]),
            ]
        )

    def _build_output_node(self):
        """Build and connect an output node to the pipeline."""
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from nipype.interfaces.io import DataSink

        from clinica.utils.nipype import container_from_filename, fix_join

        from .tasks import rename_into_caps_task

        write_node = npe.Node(name="WriteCaps", interface=DataSink())
        write_node.inputs.base_directory = str(self.caps_directory)
        write_node.inputs.parameterization = False

        rename_file_node_inputs = [
            "pet_bids_image_filename",
            "pet_preprocessed_image_filename",
            "pet_to_mri_transformation_filename",
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
                function=rename_into_caps_task,
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
        self.connect(
            [
                (self.input_node, container_path, [("pet", "bids_or_caps_filename")]),
                (
                    container_path,
                    write_node,
                    [(("container", fix_join, "pet_linear"), "container")],
                ),
                (self.input_node, rename_files, [("pet", "pet_bids_image_filename")]),
                (
                    self.output_node,
                    rename_files,
                    [("affine_mat", "pet_to_mri_transformation_filename")],
                ),
                (
                    rename_files,
                    write_node,
                    [("transformation_filename_caps", "@transform_mat")],
                ),
            ]
        )
        if not (self.parameters.get("uncropped_image")):
            node_out_name = "outfile_crop"
        else:
            node_out_name = "suvr_pet"
        self.connect(
            [
                (
                    self.output_node,
                    rename_files,
                    [(node_out_name, "pet_preprocessed_image_filename")],
                ),
                (
                    rename_files,
                    write_node,
                    [("pet_filename_caps", "@registered_pet")],
                ),
            ]
        )
        if self.parameters.get("save_PETinT1w"):
            self.connect(
                [
                    (
                        self.output_node,
                        rename_files,
                        [("PETinT1w", "pet_filename_in_t1w_raw")],
                    ),
                    (
                        rename_files,
                        write_node,
                        [("pet_filename_in_t1w_caps", "@registered_pet_in_t1w")],
                    ),
                ]
            )

    def _build_core_nodes(self):
        """Build and connect the core nodes of the pipeline."""
        import nipype.interfaces.utility as nutil
        from nipype.interfaces import ants

        from clinica.pipelines.tasks import (
            crop_nifti_task,
            get_filename_no_ext_task,
            run_ants_apply_transforms_task,
        )

        from .tasks import (
            clip_task,
            perform_suvr_normalization_task,
        )
        from .utils import concatenate_transforms, init_input_node, print_end_pipeline

        init_node = npe.Node(
            interface=nutil.Function(
                input_names=["pet"],
                output_names=["pet"],
                function=init_input_node,
            ),
            name="initPipeline",
        )
        # image_id_node = npe.Node(
        #    interface=nutil.Function(
        #        input_names=["filename"],
        #        output_names=["image_id"],
        #        function=get_filename_no_ext_task,
        #    ),
        #    name="ImageID",
        # )
        concatenate_node = npe.Node(
            interface=nutil.Function(
                input_names=["pet_to_t1w_transform", "t1w_to_mni_transform"],
                output_names=["transforms_list"],
                function=concatenate_transforms,
            ),
            name="concatenateTransforms",
        )

        # 1. Clipping node
        clipping_node = npe.Node(
            name="clipping",
            interface=nutil.Function(
                function=clip_task,
                input_names=["input_pet", "output_dir"],
                output_names=["output_image"],
            ),
        )
        clipping_node.inputs.output_dir = self.base_dir

        ants_registration_node = self._build_ants_registration_node()
        ants_applytransform_node = npe.Node(
            name="antsApplyTransformPET2MNI",
            interface=(
                nutil.Function(
                    function=run_ants_apply_transforms_task,
                    input_names=[
                        "reference_image",
                        "input_image",
                        "transforms",
                        "output_dir",
                    ],
                    output_names=["output_image"],
                )
                if self.use_antspy
                else ants.ApplyTransforms()
            ),
        )
        if not self.use_antspy:
            ants_applytransform_node.inputs.dimension = 3
        ants_applytransform_node.inputs.reference_image = self.ref_template

        ants_registration_nonlinear_node = (
            self._build_ants_registration_nonlinear_node()
        )

        ants_applytransform_nonlinear_node = npe.Node(
            name="antsApplyTransformNonLinear",
            interface=(
                nutil.Function(
                    function=run_ants_apply_transforms_task,
                    input_names=[
                        "reference_image",
                        "input_image",
                        "transforms",
                        "output_dir",
                    ],
                    output_names=["output_image"],
                )
                if self.use_antspy
                else ants.ApplyTransforms()
            ),
        )
        if not self.use_antspy:
            ants_applytransform_nonlinear_node.inputs.dimension = 3
        ants_applytransform_nonlinear_node.inputs.reference_image = self.ref_template

        normalize_intensity_node = npe.Node(
            name="intensityNormalization",
            interface=nutil.Function(
                function=perform_suvr_normalization_task,
                input_names=[
                    "pet_image_path",
                    "normalizing_image_path",
                    "reference_mask_path",
                ],
                output_names=["output_image"],
            ),
        )
        normalize_intensity_node.inputs.reference_mask_path = self.ref_mask

        # 5. Crop image (using nifti). It uses custom interface, from utils file
        crop_nifti_node = npe.Node(
            name="cropNifti",
            interface=nutil.Function(
                function=crop_nifti_task,
                input_names=["input_image", "output_path"],
                output_names=["output_image"],
            ),
        )
        crop_nifti_node.inputs.output_path = self.base_dir

        # 6. Print end message
        print_end_message = npe.Node(
            interface=nutil.Function(
                input_names=["pet", "final_file"], function=print_end_pipeline
            ),
            name="WriteEndMessage",
        )

        # 7. Optional node: compute PET image in T1w
        ants_applytransform_optional_node = npe.Node(
            name="antsApplyTransformPET2T1w",
            interface=(
                nutil.Function(
                    function=run_ants_apply_transforms_task,
                    input_names=[
                        "reference_image",
                        "input_image",
                        "transforms",
                        "output_dir",
                    ],
                    output_names=["output_image"],
                )
                if self.use_antspy
                else ants.ApplyTransforms()
            ),
        )
        if not self.use_antspy:
            ants_applytransform_optional_node.inputs.dimension = 3

        self.connect(
            [
                (self.input_node, init_node, [("pet", "pet")]),
                # STEP 1:
                (init_node, clipping_node, [("pet", "input_pet")]),
                # STEP 2
                (
                    clipping_node,
                    ants_registration_node,
                    [("output_image", "moving_image")],
                ),
                (self.input_node, ants_registration_node, [("t1w", "fixed_image")]),
                # STEP 3
                (
                    ants_registration_node,
                    concatenate_node,
                    [("out_matrix", "pet_to_t1w_transform")],
                ),
                (
                    self.input_node,
                    concatenate_node,
                    [("t1w_to_mni", "t1w_to_mni_transform")],
                ),
                (
                    clipping_node,
                    ants_applytransform_node,
                    [("output_image", "input_image")],
                ),
                (
                    concatenate_node,
                    ants_applytransform_node,
                    [("transforms_list", "transforms")],
                ),
                # STEP 4
                (
                    self.input_node,
                    ants_registration_nonlinear_node,
                    [("t1w_linear", "moving_image")],
                ),
                (
                    ants_registration_nonlinear_node,
                    ants_applytransform_nonlinear_node,
                    [("forward_transforms", "transforms")],
                ),
                (
                    ants_applytransform_node,
                    ants_applytransform_nonlinear_node,
                    [("output_image", "input_image")],
                ),
                (
                    ants_applytransform_node,
                    normalize_intensity_node,
                    [("output_image", "pet_image_path")],
                ),
                (
                    ants_applytransform_nonlinear_node,
                    normalize_intensity_node,
                    [("output_image", "normalizing_image_path")],
                ),
                # Connect to DataSink
                (
                    ants_registration_node,
                    self.output_node,
                    [("out_matrix", "affine_mat")],
                ),
                (
                    normalize_intensity_node,
                    self.output_node,
                    [("output_image", "suvr_pet")],
                ),
                (self.input_node, print_end_message, [("pet", "pet")]),
            ]
        )
        # STEP 5
        # Case 1: crop the image
        if not (self.parameters.get("uncropped_image")):
            self.connect(
                [
                    (
                        normalize_intensity_node,
                        crop_nifti_node,
                        [("output_image", "input_image")],
                    ),
                    (
                        crop_nifti_node,
                        self.output_node,
                        [("output_image", "outfile_crop")],
                    ),
                ]
            )
            last_node = crop_nifti_node
        # Case 2:  don't crop the image
        else:
            last_node = normalize_intensity_node
        self.connect(
            [
                (
                    last_node,
                    print_end_message,
                    [("output_image", "final_file")],
                ),
            ]
        )

        # STEP 6: Optional argument
        if self.parameters.get("save_PETinT1w"):
            self.connect(
                [
                    (
                        clipping_node,
                        ants_applytransform_optional_node,
                        [("output_image", "input_image"), ("t1w", "reference_image")],
                    ),
                    (
                        ants_registration_node,
                        ants_applytransform_optional_node,
                        [("out_matrix", "transforms")],
                    ),
                    (
                        ants_applytransform_optional_node,
                        self.output_node,
                        [("output_image", "PETinT1w")],
                    ),
                ]
            )

    def _build_ants_registration_node(self) -> npe.Node:
        import nipype.interfaces.utility as nutil
        from nipype.interfaces import ants

        from clinica.pipelines.tasks import run_ants_registration_synquick_task
        from clinica.pipelines.utils import AntsRegistrationSynQuickTransformType

        ants_registration_node = npe.Node(
            name="antsRegistration",
            interface=(
                nutil.Function(
                    function=run_ants_registration_synquick_task,
                    input_names=[
                        "fixed_image",
                        "moving_image",
                        "random_seed",
                        "transform_type",
                        "output_prefix",
                        "output_dir",
                    ],
                    output_names=["warped_image", "out_matrix"],
                )
                if self.use_antspy
                else ants.RegistrationSynQuick()
            ),
        )
        ants_registration_node.inputs.fixed_image = self.ref_template
        if self.use_antspy:
            ants_registration_node.inputs.output_dir = str(self.base_dir)
            ants_registration_node.inputs.transform_type = (
                AntsRegistrationSynQuickTransformType.RIGID
            )
        else:
            ants_registration_node.inputs.transform_type = "r"
            ants_registration_node.inputs.dimension = 3
        ants_registration_node.inputs.random_seed = (
            self.parameters.get("random_seed", None) or 0
        )

        return ants_registration_node

    def _build_ants_registration_nonlinear_node(self) -> npe.Node:
        import nipype.interfaces.utility as nutil
        from nipype.interfaces import ants

        from clinica.pipelines.tasks import run_ants_registration_task
        from clinica.pipelines.utils import AntsRegistrationTransformType

        ants_registration_nonlinear_node = npe.Node(
            name="antsRegistrationT1W2MNI",
            interface=(
                nutil.Function(
                    function=run_ants_registration_task,
                    input_names=[
                        "fixed_image",
                        "moving_image",
                        "random_seed",
                        "transform_type",
                        "output_prefix",
                        "output_dir",
                        "shrink_factors",
                        "smoothing_sigmas",
                        "number_of_iterations",
                        "return_inverse_transform",
                    ],
                    output_names=[
                        "warped_image",
                        "out_matrix",
                        "reverse_forward_transforms",
                    ],
                )
                if self.use_antspy
                else ants.Registration()
            ),
        )
        ants_registration_nonlinear_node.inputs.fixed_image = self.ref_template
        if self.use_antspy:
            ants_registration_nonlinear_node.inputs.transform_type = (
                AntsRegistrationTransformType.SYN
            )
            ants_registration_nonlinear_node.inputs.shrink_factors = (8, 4, 2)
            ants_registration_nonlinear_node.inputs.smoothing_sigmas = (3, 2, 1)
            ants_registration_nonlinear_node.inputs.number_of_iterations = (200, 50, 10)
            ants_registration_nonlinear_node.inputs.return_inverse_transform = True
        else:
            ants_registration_nonlinear_node.inputs.metric = ["MI"]
            ants_registration_nonlinear_node.inputs.metric_weight = [1.0]
            ants_registration_nonlinear_node.inputs.transforms = ["SyN"]
            ants_registration_nonlinear_node.inputs.dimension = 3
            ants_registration_nonlinear_node.inputs.shrink_factors = [[8, 4, 2]]
            ants_registration_nonlinear_node.inputs.smoothing_sigmas = [[3, 2, 1]]
            ants_registration_nonlinear_node.inputs.sigma_units = ["vox"]
            ants_registration_nonlinear_node.inputs.number_of_iterations = [
                [200, 50, 10]
            ]
            ants_registration_nonlinear_node.inputs.radius_or_number_of_bins = [32]

            ants_registration_nonlinear_node.inputs.transform_parameters = [(0.1, 3, 0)]
            ants_registration_nonlinear_node.inputs.convergence_threshold = [1e-05]
            ants_registration_nonlinear_node.inputs.convergence_window_size = [10]
            ants_registration_nonlinear_node.inputs.winsorize_lower_quantile = 0.005
            ants_registration_nonlinear_node.inputs.winsorize_upper_quantile = 0.995
            ants_registration_nonlinear_node.inputs.collapse_output_transforms = True
            ants_registration_nonlinear_node.inputs.use_histogram_matching = False
        ants_registration_nonlinear_node.inputs.verbose = True
        ants_registration_nonlinear_node.inputs.random_seed = (
            self.parameters.get("random_seed", None) or 0
        )

        return ants_registration_nonlinear_node
