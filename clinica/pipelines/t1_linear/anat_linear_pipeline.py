# Use hash instead of parameters for iterables folder names
# Otherwise path will be too long and generate OSError
from pathlib import Path
from typing import List, Optional

from nipype import config

from clinica.pipelines.engine import Pipeline
from clinica.utils.bids import Visit
from clinica.utils.check_dependency import ThirdPartySoftware
from clinica.utils.stream import log_and_warn

cfg = dict(execution={"parameterize_dirs": False})
config.update_config(cfg)


class AnatLinear(Pipeline):
    """Anat Linear - Affine registration of anat (t1w or flair) images to standard space.

    This preprocessing pipeline includes globally three steps:
    1) Bias correction with N4 algorithm from ANTs.
    2) Linear registration to MNI152NLin2009cSym template with
       RegistrationSynQuick from ANTs.
    3) Crop the background (in order to save computational power).

    Returns:
        A clinica pipeline object containing the  AnatLinear pipeline.
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
        caps_name: Optional[str] = None,
    ):
        self.use_antspy = use_antspy
        if self.use_antspy:
            if ignore_dependencies:
                ignore_dependencies.append(ThirdPartySoftware.ANTS)
            else:
                ignore_dependencies = [ThirdPartySoftware.ANTS]
            log_and_warn(
                (
                    "The AnatLinear pipeline has been configured to use ANTsPy instead of ANTs.\n"
                    "This means that no installation of ANTs is required, but the antspyx Python "
                    "package must be installed in your environment.\nThis functionality has been "
                    "introduced in Clinica 0.9.0 and is considered experimental.\n"
                    "Please report any issue or unexpected results to the Clinica developer team."
                ),
                UserWarning,
            )
        super().__init__(
            bids_directory=bids_directory,
            caps_directory=caps_directory,
            tsv_file=tsv_file,
            overwrite_caps=overwrite_caps,
            base_dir=base_dir,
            parameters=parameters,
            ignore_dependencies=ignore_dependencies,
            name=name,
            caps_name=caps_name,
        )

    def get_processed_visits(self) -> list[Visit]:
        from clinica.utils.filemanip import extract_visits
        from clinica.utils.input_files import T1W_LINEAR, T1W_LINEAR_CROPPED
        from clinica.utils.inputs import clinica_file_reader

        processed_visits: list[Visit] = []
        if self.caps_directory.is_dir():
            cropped_files, _ = clinica_file_reader(
                self.subjects,
                self.sessions,
                self.caps_directory,
                T1W_LINEAR
                if self.parameters.get("uncropped_image", False)
                else T1W_LINEAR_CROPPED,
            )
            processed_visits.extend(extract_visits(cropped_files))
        return processed_visits

    def _check_custom_dependencies(self) -> None:
        """Check dependencies that can not be listed in the `info.json` file."""
        pass

    def _check_pipeline_parameters(self) -> None:
        """Check pipeline parameters."""
        pass

    def get_input_fields(self) -> List[str]:
        """Specify the list of possible inputs of this pipeline.

        Returns
        -------
        list of str :
            A list of (string) input fields name.
        """
        return ["anat"]

    def get_output_fields(self) -> List[str]:
        """Specify the list of possible outputs of this pipeline.

        Returns
        -------
        list of str:
            A list of (string) output fields name.
        """
        return ["image_id"]

    def _build_input_node(self):
        """Build and connect an input node to the pipeline."""
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from clinica.utils.image import get_mni_template
        from clinica.utils.input_files import T1W_NII, Flair_T2W_NII
        from clinica.utils.inputs import clinica_file_filter
        from clinica.utils.stream import cprint
        from clinica.utils.ux import print_images_to_process

        self.ref_template = get_mni_template(
            "t1" if self.name == "t1-linear" else "flair"
        )

        # Inputs from anat/ folder
        # ========================
        # anat image file:
        query = T1W_NII if self.name == "t1-linear" else Flair_T2W_NII

        anat_files, filtered_subjects, filtered_sessions = clinica_file_filter(
            self.subjects, self.sessions, self.bids_directory, query
        )
        self.subjects = filtered_subjects
        self.sessions = filtered_sessions

        if len(self.subjects):
            print_images_to_process(self.subjects, self.sessions)
            cprint("The pipeline will last approximately 6 minutes per image.")

        read_node = npe.Node(
            name="ReadingFiles",
            iterables=[
                ("anat", anat_files),
            ],
            synchronize=True,
            interface=nutil.IdentityInterface(fields=self.get_input_fields()),
        )
        self.connect(
            [
                (read_node, self.input_node, [("anat", "anat")]),
            ]
        )

    def _build_output_node(self):
        """Build and connect an output node to the pipeline."""
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from nipype.interfaces.io import DataSink

        from clinica.utils.nipype import container_from_filename, fix_join

        if self.name == "flair-linear":
            from .anat_linear_utils import (
                get_substitutions_datasink_flair as get_substitutions,
            )
        else:
            from .anat_linear_utils import (
                get_substitutions_datasink_t1_linear as get_substitutions,
            )

        # Writing node
        write_node = npe.Node(name="WriteCaps", interface=DataSink())
        write_node.inputs.base_directory = str(self.caps_directory)
        write_node.inputs.parameterization = False

        # Other nodes
        # =====================================
        # Get substitutions to rename files
        get_ids = npe.Node(
            interface=nutil.Function(
                input_names=["bids_image_id"],
                output_names=["substitutions"],
                function=get_substitutions,
            ),
            name="GetIDs",
        )
        # Find container path from t1w filename
        container_path = npe.Node(
            nutil.Function(
                input_names=["bids_or_caps_filename"],
                output_names=["container"],
                function=container_from_filename,
            ),
            name="ContainerPath",
        )
        self.connect(
            [
                (self.input_node, container_path, [("anat", "bids_or_caps_filename")]),
                (self.output_node, get_ids, [("image_id", "bids_image_id")]),
                (
                    container_path,
                    write_node,
                    [
                        (
                            ("container", fix_join, self.name.replace("-", "_")),
                            "container",
                        )
                    ],
                ),
                (get_ids, write_node, [("substitutions", "substitutions")]),
                (self.output_node, write_node, [("image_id", "@image_id")]),
                (self.output_node, write_node, [("outfile_reg", "@outfile_reg")]),
                (self.output_node, write_node, [("affine_mat", "@affine_mat")]),
            ]
        )

        if not (self.parameters.get("uncropped_image")):
            self.connect(
                [
                    (self.output_node, write_node, [("outfile_crop", "@outfile_crop")]),
                ]
            )

    def _build_core_nodes(self):
        """Build and connect the core nodes of the pipeline."""
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from nipype.interfaces import ants

        from clinica.pipelines.t1_linear.tasks import (
            run_ants_registration_task,
            run_n4biasfieldcorrection_task,
        )
        from clinica.pipelines.tasks import crop_nifti_task, get_filename_no_ext_task

        from .anat_linear_utils import print_end_pipeline

        image_id_node = npe.Node(
            interface=nutil.Function(
                input_names=["filename"],
                output_names=["image_id"],
                function=get_filename_no_ext_task,
            ),
            name="ImageID",
        )

        # 1. N4biascorrection by ANTS. It uses nipype interface.
        n4biascorrection = npe.Node(
            name="n4biascorrection",
            interface=(
                nutil.Function(
                    function=run_n4biasfieldcorrection_task,
                    input_names=[
                        "input_image",
                        "bspline_fitting_distance",
                        "output_prefix",
                        "output_dir",
                        "save_bias",
                        "verbose",
                    ],
                    output_names=["output_image"],
                )
                if self.use_antspy
                else ants.N4BiasFieldCorrection(dimension=3)
            ),
        )
        n4biascorrection.inputs.save_bias = True
        if self.use_antspy:
            n4biascorrection.inputs.output_dir = str(self.base_dir)
            n4biascorrection.inputs.verbose = True
        if self.name == "t1-linear":
            n4biascorrection.inputs.bspline_fitting_distance = 600
        else:
            n4biascorrection.inputs.bspline_fitting_distance = 100

        # 2. `RegistrationSynQuick` by *ANTS*. It uses nipype interface.
        ants_registration_node = npe.Node(
            name="antsRegistrationSynQuick",
            interface=(
                nutil.Function(
                    function=run_ants_registration_task,
                    input_names=[
                        "fixed_image",
                        "moving_image",
                        "random_seed",
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
        if not self.use_antspy:
            ants_registration_node.inputs.transform_type = "a"
            ants_registration_node.inputs.dimension = 3

        random_seed = self.parameters.get("random_seed", None)
        ants_registration_node.inputs.random_seed = random_seed or 0

        # 3. Crop image (using nifti). It uses custom interface, from utils file

        cropnifti = npe.Node(
            name="cropnifti",
            interface=nutil.Function(
                function=crop_nifti_task,
                input_names=["input_image", "output_path"],
                output_names=["output_img"],
            ),
        )
        cropnifti.inputs.output_path = self.base_dir

        # 4. Print end message
        print_end_message = npe.Node(
            interface=nutil.Function(
                input_names=["anat", "final_file"], function=print_end_pipeline
            ),
            name="WriteEndMessage",
        )
        self.connect(
            [
                (self.input_node, image_id_node, [("anat", "filename")]),
                (self.input_node, n4biascorrection, [("anat", "input_image")]),
                (
                    n4biascorrection,
                    ants_registration_node,
                    [("output_image", "moving_image")],
                ),
                (
                    image_id_node,
                    ants_registration_node,
                    [("image_id", "output_prefix")],
                ),
                # Connect to DataSink
                (image_id_node, self.output_node, [("image_id", "image_id")]),
                (
                    ants_registration_node,
                    self.output_node,
                    [("out_matrix", "affine_mat")],
                ),
                (
                    ants_registration_node,
                    self.output_node,
                    [("warped_image", "outfile_reg")],
                ),
                (self.input_node, print_end_message, [("anat", "anat")]),
            ]
        )
        if self.use_antspy:
            self.connect(
                [
                    (
                        image_id_node,
                        n4biascorrection,
                        [("image_id", "output_prefix")],
                    ),
                ]
            )
        if not (self.parameters.get("uncropped_image")):
            self.connect(
                [
                    (
                        ants_registration_node,
                        cropnifti,
                        [("warped_image", "input_image")],
                    ),
                    (cropnifti, self.output_node, [("output_img", "outfile_crop")]),
                    (cropnifti, print_end_message, [("output_img", "final_file")]),
                ]
            )
        else:
            self.connect(
                [
                    (
                        ants_registration_node,
                        print_end_message,
                        [("warped_image", "final_file")],
                    ),
                ]
            )
