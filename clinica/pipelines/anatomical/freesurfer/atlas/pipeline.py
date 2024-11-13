from pathlib import Path
from typing import List, Optional

from nipype import config

from clinica.pipelines.engine import Pipeline

cfg = dict(execution={"parameterize_dirs": False})
config.update_config(cfg)


class T1FreeSurferAtlas(Pipeline):
    """Projection of the results of t1-freesurfer on another atlass

    Returns:
        A clinica pipeline object containing the T1FreeSurfer atlas projection pipeline.
    """

    def __init__(
        self,
        caps_directory: str,
        atlas_path: Optional[str] = None,
    ):
        if atlas_path:
            self.atlas_path = Path(atlas_path)
        else:
            self.atlas_path = None
        super().__init__(caps_directory=caps_directory)

    def _check_custom_dependencies(self) -> None:
        """Check dependencies provided by the developer."""
        pass

    def _check_pipeline_parameters(self) -> None:
        """Check pipeline parameters."""
        pass

    @staticmethod
    def get_to_process_with_atlases(
        caps_directory: Path,
        subjects: List[str],
        sessions: List[str],
        atlas_dir_path: Path,
    ) -> List[str]:
        import itertools

        from clinica.pipelines.t1_freesurfer_longitudinal.longitudinal_utils import (
            grab_image_ids_from_caps_directory,
        )
        from clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_longitudinal_correction_utils import (
            get_processed_images,
        )
        from clinica.utils.filemanip import extract_image_ids
        from clinica.utils.image import HemiSphere
        from clinica.utils.input_files import (
            Parcellation,
            get_t1_freesurfer_segmentation,
            get_t1_freesurfer_statistics,
        )
        from clinica.utils.inputs import clinica_file_reader

        part_ids, sess_ids, list_long_id = grab_image_ids_from_caps_directory(
            caps_directory
        )

        initial_list_to_process = []
        atlas_list = []
        for path in atlas_dir_path.rglob("*rh*6p0.gcs"):
            atlas_name = path.name.split(".")[1].split("_")[0]
            atlas_list.append(atlas_name)

        if caps_directory.is_dir():
            for atlas in atlas_list:
                t1_freesurfer_longitudinal_output = get_processed_images(
                    caps_directory, part_ids, sess_ids, list_long_id
                )
                t1_freesurfer_longitudinal_output_atlas = get_processed_images(
                    caps_directory, part_ids, sess_ids, list_long_id, atlas
                )
                to_process_long = list(
                    set(t1_freesurfer_longitudinal_output)
                    - set(t1_freesurfer_longitudinal_output_atlas)
                )
                t1_freesurfer_output, _ = clinica_file_reader(
                    subjects,
                    sessions,
                    caps_directory,
                    get_t1_freesurfer_segmentation(Parcellation.DESTRIEUX),
                )
                t1_freesurfer_files, _ = clinica_file_reader(
                    subjects,
                    sessions,
                    caps_directory,
                    get_t1_freesurfer_statistics(atlas, HemiSphere.RIGHT),
                )
                image_ids = extract_image_ids(t1_freesurfer_files)
                image_ids_2 = extract_image_ids(t1_freesurfer_output)
                to_process = list(set(image_ids_2) - set(image_ids)) + to_process_long
                initial_list_to_process.append(([atlas], to_process))

        list_to_process = []
        for i in initial_list_to_process:
            if list(itertools.product(i[0], i[1])) != []:
                list_to_process = list_to_process + list(itertools.product(i[0], i[1]))
        return list_to_process

    def get_input_fields(self) -> List[str]:
        """Specify the list of possible inputs of this pipeline.

        Notes
        -----
        The list of inputs of the T1FreeSurferAtlas pipeline is:
            * to_process_with_atlases (str): list of the tuples (atlas, sub-ses) to process

        Returns
        -------
        list of str :
            A list of (string) input fields name.
        """
        return ["to_process_with_atlases"]

    def get_output_fields(self) -> List[str]:
        """Specify the list of possible outputs of this pipeline.

        Notes
        -----
        The list of outputs of the T1FreeSurfer pipeline is:
            * image_id (str): Image ID (e.g. sub-CLNC01_ses-M000)

        Returns
        -------
        list of str :
            A list of (string) output fields name.
        """
        return ["image_id"]

    def _build_input_node(self):
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        to_process_with_atlases = self.get_to_process_with_atlases(
            self.caps_directory,
            self.subjects,
            self.sessions,
            atlas_dir_path=self.atlas_path,
        )
        read_node = npe.Node(
            name="ReadingFiles",
            iterables=[
                ("to_process_with_atlases", to_process_with_atlases),
            ],
            synchronize=True,
            interface=nutil.IdentityInterface(fields=self.get_input_fields()),
        )
        self.connect(
            [
                (
                    read_node,
                    self.input_node,
                    [("to_process_with_atlases", "to_process_with_atlases")],
                )
            ]
        )

    def _build_output_node(self):
        return super().build_output_node()

    def _build_core_nodes(self):
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from .tasks import compute_atlas_task, write_tsv_files_task

        # Run an additional two Freesurfer commands if there is an atlas_path specified
        compute_other_atlases = npe.Node(
            interface=nutil.Function(
                input_names=[
                    "caps_directory",
                    "to_process_with_atlases",
                    "path_to_atlas",
                ],
                output_names=["subject_dir", "image_id", "atlas"],
                function=compute_atlas_task,
            ),
            name="1-ComputeOtherAtlases",
        )
        compute_other_atlases.inputs.path_to_atlas = str(self.atlas_path)
        compute_other_atlases.inputs.caps_directory = str(self.caps_directory)

        create_tsv = npe.Node(
            interface=nutil.Function(
                input_names=["subject_dir", "image_id", "atlas"],
                output_names=["image_id"],
                function=write_tsv_files_task,
            ),
            name="2-CreateTsvFiles",
        )

        self.connect(
            [
                (
                    self.input_node,
                    compute_other_atlases,
                    [("to_process_with_atlases", "to_process_with_atlases")],
                ),
                (compute_other_atlases, create_tsv, [("subject_dir", "subject_dir")]),
                (compute_other_atlases, create_tsv, [("image_id", "image_id")]),
                (compute_other_atlases, create_tsv, [("atlas", "atlas")]),
            ]
        )
