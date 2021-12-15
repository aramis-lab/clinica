# coding: utf8

# Use hash instead of parameters for iterables folder names
# Otherwise path will be too long and generate OSError
from typing import Optional

from nipype import config

import clinica.pipelines.engine as cpe
from clinica.pipelines import t1_freesurfer_longitudinal
from clinica.pipelines.t1_freesurfer_longitudinal import (
    t1_freesurfer_longitudinal_correction_utils,
)
from clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_longitudinal_correction_utils import (
    get_processed_images,
)

cfg = dict(execution={"parameterize_dirs": False})
config.update_config(cfg)


class T1FreeSurferAtlas(cpe.Pipeline):
    """Projection of the results of t1-freesurfer on another atlass

    Returns:
        A clinica pipeline object containing the T1FreeSurfer atlas projection pipeline.
    """

    def __init__(
        self,
        caps_directory: str,
        atlas_path: Optional[str] = None,
    ):
        self.atlas_path = atlas_path
        super().__init__(
            caps_directory=caps_directory,
        )

    @staticmethod
    def get_to_process_with_atlases(
        caps_directory: str, subjects: list, sessions: list, atlas_dir_path: str
    ) -> list:
        import itertools
        import os
        from pathlib import Path

        from clinica.pipelines.t1_freesurfer_longitudinal.longitudinal_utils import (
            grab_image_ids_from_caps_directory,
        )
        from clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_longitudinal_correction_utils import (
            get_processed_images,
        )
        from clinica.utils.filemanip import extract_image_ids
        from clinica.utils.input_files import T1_FS_DESTRIEUX
        from clinica.utils.inputs import clinica_file_reader

        part_ids, sess_ids, list_long_id = grab_image_ids_from_caps_directory(
            caps_directory
        )

        initial_list_to_process = []
        atlas_list = []
        for path in Path(atlas_dir_path).rglob("*rh*6p0.gcs"):
            atlas_name = path.name.split(".")[1].split("_")[0]
            atlas_list.append(atlas_name)

        if os.path.isdir(caps_directory):
            for atlas in atlas_list:

                atlas_info = dict(
                    {
                        "pattern": "t1/freesurfer_cross_sectional/sub-*_ses-*/stats/rh."
                        + atlas
                        + ".stats",
                        "description": atlas + "-based segmentation",
                        "needed_pipeline": "t1-freesurfer",
                    }
                )
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
                    subjects, sessions, caps_directory, T1_FS_DESTRIEUX, False
                )
                t1_freesurfer_files, _ = clinica_file_reader(
                    subjects, sessions, caps_directory, atlas_info, False
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

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Note:
            The list of inputs of the T1FreeSurferAtlas pipeline is:
                * to_process_with_atlases (str): list of the tuples (atlas, sub-ses) to process

        Returns:
            A list of (string) input fields name.
        """
        return ["to_process_with_atlases"]

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Note:
            The list of outputs of the T1FreeSurfer pipeline is:
                * image_id (str): Image ID (e.g. sub-CLNC01_ses-M00)

        Returns:
            A list of (string) output fields name.
        """
        return ["image_id"]

    def build_input_node(self):
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
        # fmt: off
        self.connect(
            [(read_node, self.input_node, [("to_process_with_atlases", "to_process_with_atlases")])]
        )
        # fmt: on

    def build_output_node(self):
        return super().build_output_node()

    def build_core_nodes(self):
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from .t1_freesurfer_atlas_utils import compute_atlas, write_tsv_files

        # Run an additional two Freesurfer commands if there is an atlas_path specified
        compute_other_atlases = npe.Node(
            interface=nutil.Function(
                input_names=[
                    "caps_directory",
                    "to_process_with_atlases",
                    "path_to_atlas",
                ],
                output_names=["subject_dir", "image_id", "atlas"],
                function=compute_atlas,
            ),
            name="1-ComputeOtherAtlases",
        )
        compute_other_atlases.inputs.path_to_atlas = self.atlas_path
        compute_other_atlases.inputs.caps_directory = self.caps_directory

        create_tsv = npe.Node(
            interface=nutil.Function(
                input_names=["subject_dir", "image_id", "atlas"],
                output_names=["image_id"],
                function=write_tsv_files,
            ),
            name="2-CreateTsvFiles",
        )

        self.connect(
            [
                # Run compute_atlases command
                (
                    self.input_node,
                    compute_other_atlases,
                    [("to_process_with_atlases", "to_process_with_atlases")],
                ),
                # Generate TSV files
                (compute_other_atlases, create_tsv, [("subject_dir", "subject_dir")]),
                (compute_other_atlases, create_tsv, [("image_id", "image_id")]),
                (compute_other_atlases, create_tsv, [("atlas", "atlas")]),
            ]
        )
