# Use hash instead of parameters for iterables folder names
# Otherwise path will be too long and generate OSError
from nipype import config

import clinica.pipelines.engine as cpe

cfg = dict(execution={"parameterize_dirs": False})
config.update_config(cfg)


class AnatLinear(cpe.Pipeline):
    """Anat Linear - Affine registration of anat (t1w or flair) images to standard space.

    This preprocessing pipeline includes globally three steps:
    1) Bias correction with N4 algorithm from ANTs.
    2) Linear registration to MNI152NLin2009cSym template with
       RegistrationSynQuick from ANTs.
    3) Crop the background (in order to save computational power).

    Returns:
        A clinica pipeline object containing the  AnatLinear pipeline.
    """

    @staticmethod
    def get_processed_images(caps_directory, subjects, sessions):
        import os

        from clinica.utils.filemanip import extract_image_ids
        from clinica.utils.input_files import T1W_LINEAR_CROPPED
        from clinica.utils.inputs import clinica_file_reader

        image_ids = []
        if os.path.isdir(caps_directory):
            cropped_files, _ = clinica_file_reader(
                subjects, sessions, caps_directory, T1W_LINEAR_CROPPED, False
            )
            image_ids = extract_image_ids(cropped_files)
        return image_ids

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file."""

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Returns:
            A list of (string) input fields name.
        """
        return ["anat"]

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """
        return ["image_id"]

    def build_input_node(self):
        """Build and connect an input node to the pipeline."""
        from os import pardir
        from os.path import abspath, dirname, exists, join

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from clinica.utils.exceptions import ClinicaBIDSError, ClinicaException
        from clinica.utils.filemanip import extract_subjects_sessions_from_filename
        from clinica.utils.input_files import T1W_NII, Flair_T2W_NII
        from clinica.utils.inputs import (
            RemoteFileStructure,
            clinica_file_reader,
            fetch_file,
        )
        from clinica.utils.stream import cprint
        from clinica.utils.ux import print_images_to_process

        root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
        path_to_mask = join(root, "resources", "masks")
        if self.name == "t1-linear":
            url_aramis = "https://aramislab.paris.inria.fr/files/data/img_t1_linear/"
            FILE2 = RemoteFileStructure(
                filename="mni_icbm152_t1_tal_nlin_sym_09c.nii",
                url=url_aramis,
                checksum="93359ab97c1c027376397612a9b6c30e95406c15bf8695bd4a8efcb2064eaa34",
            )
        else:
            url_aramis = "https://aramislab.paris.inria.fr/files/data/img_flair_linear/"
            FILE2 = RemoteFileStructure(
                filename="GG-853-FLAIR-1.0mm.nii.gz",
                url=url_aramis,
                checksum="b1d2d359a4c3671685227bb14014ce50ac232012b628335a4c049e2911c64ce1",
            )

        FILE1 = RemoteFileStructure(
            filename="ref_cropped_template.nii.gz",
            url=url_aramis,
            checksum="67e1e7861805a8fd35f7fcf2bdf9d2a39d7bcb2fd5a201016c4d2acdd715f5b3",
        )

        self.ref_template = join(path_to_mask, FILE2.filename)
        self.ref_crop = join(path_to_mask, FILE1.filename)

        if not (exists(self.ref_template)):
            try:
                fetch_file(FILE2, path_to_mask)
            except IOError as err:
                cprint(
                    msg=f"Unable to download required template (mni_icbm152) for processing: {err}",
                    lvl="error",
                )

        if not (exists(self.ref_crop)):
            try:
                fetch_file(FILE1, path_to_mask)
            except IOError as err:
                cprint(
                    msg=f"Unable to download required template (ref_crop) for processing: {err}",
                    lvl="error",
                )

        # Display image(s) already present in CAPS folder
        # ===============================================
        processed_ids = self.get_processed_images(
            self.caps_directory, self.subjects, self.sessions
        )
        if len(processed_ids) > 0:
            cprint(
                msg=f"Clinica found {len(processed_ids)} image(s) already processed in CAPS directory:",
                lvl="warning",
            )
            for image_id in processed_ids:
                cprint(msg=f"{image_id.replace('_', ' | ')}", lvl="warning")
            cprint(msg=f"Image(s) will be ignored by Clinica.", lvl="warning")
            input_ids = [
                p_id + "_" + s_id for p_id, s_id in zip(self.subjects, self.sessions)
            ]
            to_process_ids = list(set(input_ids) - set(processed_ids))
            self.subjects, self.sessions = extract_subjects_sessions_from_filename(
                to_process_ids
            )

        # Inputs from anat/ folder
        # ========================
        # anat image file:
        try:
            file = T1W_NII if self.name == "t1-linear" else Flair_T2W_NII
            anat_files, _ = clinica_file_reader(
                self.subjects, self.sessions, self.bids_directory, file
            )
        except ClinicaException as e:
            err = (
                "Clinica faced error(s) while trying to read files in your BIDS directory.\n"
                + str(e)
            )
            raise ClinicaBIDSError(err)

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

    def build_output_node(self):
        """Build and connect an output node to the pipeline."""
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from nipype.interfaces.io import DataSink

        from clinica.utils.nipype import container_from_filename, fix_join

        from .anat_linear_utils import get_substitutions_datasink

        # Writing node
        write_node = npe.Node(name="WriteCaps", interface=DataSink())
        write_node.inputs.base_directory = self.caps_directory
        write_node.inputs.parameterization = False

        # Other nodes
        # =====================================
        # Get substitutions to rename files
        get_ids = npe.Node(
            interface=nutil.Function(
                input_names=["bids_file", "pipeline_name"],
                output_names=["image_id_out", "subst_ls"],
                function=get_substitutions_datasink,
            ),
            name="GetIDs",
        )
        get_ids.inputs.pipeline_name = self.name
        # Find container path from t1w filename
        container_path = npe.Node(
            nutil.Function(
                input_names=["bids_or_caps_filename"],
                output_names=["container"],
                function=container_from_filename,
            ),
            name="ContainerPath",
        )
        # fmt: off
        self.connect(
            [
                (self.input_node, container_path, [("anat", "bids_or_caps_filename")]),
                (self.output_node, get_ids, [("image_id", "bids_file")]),
                (container_path, write_node, [(("container", fix_join, self.name.replace("-", "_")), "container")]),
                (get_ids, write_node, [("subst_ls", "substitutions")]),
                (get_ids, write_node, [("image_id_out", "@image_id")]),
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
        # fmt: on

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline."""
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from nipype.interfaces import ants

        from clinica.utils.filemanip import get_filename_no_ext

        from .anat_linear_utils import crop_nifti, print_end_pipeline

        image_id_node = npe.Node(
            interface=nutil.Function(
                input_names=["filename"],
                output_names=["image_id"],
                function=get_filename_no_ext,
            ),
            name="ImageID",
        )

        # The core (processing) nodes
        # =====================================

        # 1. N4biascorrection by ANTS. It uses nipype interface.
        n4biascorrection = npe.Node(
            name="n4biascorrection",
            interface=ants.N4BiasFieldCorrection(dimension=3, save_bias=True),
        )

        if self.name == "t1-linear":
            n4biascorrection.inputs.bspline_fitting_distance = 600
        else:
            n4biascorrection.inputs.bspline_fitting_distance = 100

        # 2. `RegistrationSynQuick` by *ANTS*. It uses nipype interface.
        ants_registration_node = npe.Node(
            name="antsRegistrationSynQuick", interface=ants.RegistrationSynQuick()
        )
        ants_registration_node.inputs.fixed_image = self.ref_template
        ants_registration_node.inputs.transform_type = "a"
        ants_registration_node.inputs.dimension = 3

        # 3. Crop image (using nifti). It uses custom interface, from utils file

        cropnifti = npe.Node(
            name="cropnifti",
            interface=nutil.Function(
                function=crop_nifti,
                input_names=["input_img", "ref_crop"],
                output_names=["output_img", "crop_template"],
            ),
        )
        cropnifti.inputs.ref_crop = self.ref_crop

        # 4. Print end message
        print_end_message = npe.Node(
            interface=nutil.Function(
                input_names=["anat", "final_file"], function=print_end_pipeline
            ),
            name="WriteEndMessage",
        )

        # Connection
        # ==========
        # fmt: off
        self.connect(
            [
                (self.input_node, image_id_node, [("anat", "filename")]),
                (self.input_node, n4biascorrection, [("anat", "input_image")]),
                (n4biascorrection, ants_registration_node, [("output_image", "moving_image")]),
                (image_id_node, ants_registration_node, [("image_id", "output_prefix")]),
                # Connect to DataSink
                (image_id_node, self.output_node, [("image_id", "image_id")]),
                (ants_registration_node, self.output_node, [("out_matrix", "affine_mat")]),
                (ants_registration_node, self.output_node, [("warped_image", "outfile_reg")]),
                (self.input_node, print_end_message, [("anat", "anat")]),
            ]
        )
        if not (self.parameters.get("uncropped_image")):
            self.connect(
                [
                    (ants_registration_node, cropnifti, [("warped_image", "input_img")]),
                    (cropnifti, self.output_node, [("output_img", "outfile_crop")]),
                    (cropnifti, print_end_message, [("output_img", "final_file")]),
                ]
            )
        else:
            self.connect(
                [
                    (ants_registration_node, print_end_message, [("warped_image", "final_file")]),
                ]
            )
        # fmt: on
