# Use hash instead of parameters for iterables folder names
from typing import List

from nipype import config

from clinica.pipelines.engine import Pipeline

cfg = dict(execution={"parameterize_dirs": False})
config.update_config(cfg)


class DwiConnectome(Pipeline):
    """Connectome-based processing of corrected DWI datasets.

    Returns:
        A clinica pipeline object containing the DwiConnectome pipeline.
    """

    def _check_pipeline_parameters(self) -> None:
        """Check pipeline parameters."""
        if self.parameters["n_tracks"] < 0:
            raise ValueError(
                f"The n_tracks parameter ({self.parameters['n_tracks']}) should be positive."
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
        return [
            "t1_brain_file",
            "wm_mask_file",
            "dwi_file",
            "dwi_brainmask_file",
            "grad_fsl",
            "atlas_files",
        ]

    def get_output_fields(self) -> List[str]:
        """Specify the list of possible outputs of this pipeline.

        Returns
        -------
        list of str :
            A list of (string) output fields name.
        """
        return ["response", "fod", "tracts", "nodes", "connectomes"]

    def _build_input_node(self):
        """Build and connect an input node to the pipeline."""
        import re

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        import clinica.utils.input_files as input_files
        from clinica.utils.exceptions import ClinicaCAPSError
        from clinica.utils.filemanip import save_participants_sessions
        from clinica.utils.inputs import clinica_list_of_files_reader
        from clinica.utils.stream import cprint
        from clinica.utils.ux import print_images_to_process

        # Read CAPS files
        list_caps_files = clinica_list_of_files_reader(
            self.subjects,
            self.sessions,
            self.caps_directory,
            [
                # Inputs from t1-freesurfer pipeline
                input_files.T1_FS_WM,  # list_caps_files[0]
                input_files.T1_FS_DESIKAN,  # list_caps_files[1]
                input_files.T1_FS_DESTRIEUX,  # list_caps_files[2]
                input_files.T1_FS_BRAIN,  # list_caps_files[3]
                # Inputs from dwi-preprocessing pipeline
                input_files.DWI_PREPROC_NII,  # list_caps_files[4]
                input_files.DWI_PREPROC_BRAINMASK,  # list_caps_files[5]
                input_files.DWI_PREPROC_BVEC,  # list_caps_files[6]
                input_files.DWI_PREPROC_BVAL,  # list_caps_files[7]
            ],
            raise_exception=True,
        )

        # Check space of DWI dataset
        dwi_file_spaces = [
            re.search(
                ".*_space-(.*)_desc-preproc_dwi.nii.*", file, re.IGNORECASE
            ).group(1)
            for file in list_caps_files[4]
        ]

        # Return an error if all the DWI files are not in the same space
        if any(a != dwi_file_spaces[0] for a in dwi_file_spaces):
            raise ClinicaCAPSError(
                "Preprocessed DWI files are not all in the same space. "
                "Please process them separately using the appropriate subjects/sessions `.tsv` file (-tsv option)."
            )
        list_atlas_files = [
            [aparc_aseg, aparc_aseg_a2009]
            for aparc_aseg, aparc_aseg_a2009 in zip(
                list_caps_files[1], list_caps_files[2]
            )
        ]

        list_grad_fsl = [
            (bvec, bval) for bvec, bval in zip(list_caps_files[6], list_caps_files[7])
        ]

        # Save subjects to process in <WD>/<Pipeline.name>/participants.tsv
        save_participants_sessions(
            self.subjects, self.sessions, self.base_dir / self.name
        )

        if len(self.subjects):
            print_images_to_process(self.subjects, self.sessions)
            cprint(
                "Computational time will depend of the number of volumes in your DWI dataset and "
                "the number of streamlines you selected."
            )

        if dwi_file_spaces[0] == "b0":
            self.parameters["dwi_space"] = "b0"
            read_node = npe.Node(
                name="ReadingFiles",
                iterables=[
                    ("wm_mask_file", list_caps_files[0]),
                    ("t1_brain_file", list_caps_files[3]),
                    ("dwi_file", list_caps_files[4]),
                    ("dwi_brainmask_file", list_caps_files[5]),
                    ("grad_fsl", list_grad_fsl),
                    ("atlas_files", list_atlas_files),
                ],
                synchronize=True,
                interface=nutil.IdentityInterface(fields=self.get_input_fields()),
            )
            self.connect(
                [
                    (read_node, self.input_node, [("t1_brain_file", "t1_brain_file")]),
                    (read_node, self.input_node, [("wm_mask_file", "wm_mask_file")]),
                    (read_node, self.input_node, [("dwi_file", "dwi_file")]),
                    (
                        read_node,
                        self.input_node,
                        [("dwi_brainmask_file", "dwi_brainmask_file")],
                    ),
                    (read_node, self.input_node, [("grad_fsl", "grad_fsl")]),
                    (read_node, self.input_node, [("atlas_files", "atlas_files")]),
                ]
            )
        elif dwi_file_spaces[0] == "T1w":
            self.parameters["dwi_space"] = "T1w"
            read_node = npe.Node(
                name="ReadingFiles",
                iterables=[
                    ("wm_mask_file", list_caps_files[0]),
                    ("dwi_file", list_caps_files[4]),
                    ("dwi_brainmask_file", list_caps_files[5]),
                    ("grad_fsl", list_grad_fsl),
                    ("atlas_files", list_atlas_files),
                ],
                synchronize=True,
                interface=nutil.IdentityInterface(fields=self.get_input_fields()),
            )
            self.connect(
                [
                    (read_node, self.input_node, [("wm_mask_file", "wm_mask_file")]),
                    (read_node, self.input_node, [("dwi_file", "dwi_file")]),
                    (
                        read_node,
                        self.input_node,
                        [("dwi_brainmask_file", "dwi_brainmask_file")],
                    ),
                    (read_node, self.input_node, [("grad_fsl", "grad_fsl")]),
                    (read_node, self.input_node, [("atlas_files", "atlas_files")]),
                ]
            )
        else:
            raise ClinicaCAPSError(
                "Bad preprocessed DWI space. Please check your CAPS folder."
            )

    def _build_output_node(self):
        """Build and connect an output node to the pipeline."""
        import nipype.interfaces.io as nio
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from .utils import get_containers

        join_node = npe.JoinNode(
            name="JoinOutputs",
            joinsource="ReadingFiles",
            interface=nutil.IdentityInterface(fields=self.get_output_fields()),
        )

        write_node = npe.MapNode(
            name="WritingCAPS",
            iterfield=["container"]
            + [f"connectome_based_processing.@{o}" for o in self.get_output_fields()],
            interface=nio.DataSink(
                infields=[
                    f"connectome_based_processing.@{o}"
                    for o in self.get_output_fields()
                ]
            ),
        )
        write_node.inputs.base_directory = str(self.caps_directory)
        write_node.inputs.container = get_containers(self.subjects, self.sessions)
        write_node.inputs.substitutions = [("trait_added", "")]
        write_node.inputs.parameterization = False

        self.connect(
            [
                # Writing CAPS
                (self.output_node, join_node, [("response", "response")]),
                (self.output_node, join_node, [("fod", "fod")]),
                (self.output_node, join_node, [("tracts", "tracts")]),
                (self.output_node, join_node, [("nodes", "nodes")]),
                (self.output_node, join_node, [("connectomes", "connectomes")]),
                (
                    join_node,
                    write_node,
                    [("response", "connectome_based_processing.@response")],
                ),
                (join_node, write_node, [("fod", "connectome_based_processing.@fod")]),
                (
                    join_node,
                    write_node,
                    [("tracts", "connectome_based_processing.@tracts")],
                ),
                (
                    join_node,
                    write_node,
                    [("nodes", "connectome_based_processing.@nodes")],
                ),
                (
                    join_node,
                    write_node,
                    [("connectomes", "connectome_based_processing.@connectomes")],
                ),
            ]
        )

    def _build_core_nodes(self):
        """Build and connect the core nodes of the pipeline.

        Notes:
            - If `FSLOUTPUTTYPE` environment variable is not set, `nipype` takes
            NIFTI by default.

        Todo:
            - [x] Detect space automatically.
            - [ ] Allow for custom parcellations (See TODOs in utils).

        """
        import nipype.interfaces.freesurfer as fs
        import nipype.interfaces.fsl as fsl
        import nipype.interfaces.mrtrix3 as mrtrix3
        import nipype.interfaces.utility as niu
        import nipype.pipeline.engine as npe
        from nipype.interfaces.mrtrix.preprocess import MRTransform
        from nipype.interfaces.mrtrix3 import (
            ConstrainedSphericalDeconvolution,
            Tractography,
        )

        from clinica.utils.exceptions import ClinicaCAPSError
        from clinica.utils.mri_registration import (
            convert_flirt_transformation_to_mrtrix_transformation,
        )

        from .utils import (
            get_caps_filenames,
            get_conversion_luts,
            get_luts,
            print_begin_pipeline,
            print_end_pipeline,
        )

        # Nodes
        # =====
        # B0 Extraction (only if space=b0)
        # -------------
        split_node = npe.Node(name="Reg-0-DWI-B0Extraction", interface=fsl.Split())
        split_node.inputs.output_type = "NIFTI_GZ"
        split_node.inputs.dimension = "t"
        select_node = npe.Node(name="Reg-0-DWI-B0Selection", interface=niu.Select())
        select_node.inputs.index = 0

        # B0 Brain Extraction (only if space=b0)
        # -------------------
        mask_node = npe.Node(name="Reg-0-DWI-BrainMasking", interface=fsl.ApplyMask())
        mask_node.inputs.output_type = "NIFTI_GZ"

        # T1-to-B0 Registration (only if space=b0)
        # ---------------------
        t12b0_reg_node = npe.Node(
            name="Reg-1-T12B0Registration",
            interface=fsl.FLIRT(
                dof=6,
                interp="spline",
                cost="normmi",
                cost_func="normmi",
            ),
        )
        t12b0_reg_node.inputs.output_type = "NIFTI_GZ"

        # MGZ File Conversion (only if space=b0)
        # -------------------
        t1_brain_conv_node = npe.Node(
            name="Reg-0-T1-T1BrainConvertion", interface=fs.MRIConvert()
        )
        wm_mask_conv_node = npe.Node(
            name="Reg-0-T1-WMMaskConvertion", interface=fs.MRIConvert()
        )

        # WM Transformation (only if space=b0)
        # -----------------
        wm_transform_node = npe.Node(
            name="Reg-2-WMTransformation", interface=fsl.ApplyXFM()
        )
        wm_transform_node.inputs.apply_xfm = True

        # Nodes Generation
        # ----------------
        label_convert_node = npe.MapNode(
            name="0-LabelsConversion",
            iterfield=["in_file", "in_config", "in_lut", "out_file"],
            interface=mrtrix3.LabelConvert(),
        )
        label_convert_node.inputs.in_config = get_conversion_luts()
        label_convert_node.inputs.in_lut = get_luts()

        # FSL flirt matrix to MRtrix matrix Conversion (only if space=b0)
        # --------------------------------------------
        fsl2mrtrix_conv_node = npe.Node(
            name="Reg-2-FSL2MrtrixConversion",
            interface=niu.Function(
                input_names=[
                    "in_source_image",
                    "in_reference_image",
                    "in_flirt_matrix",
                    "name_output_matrix",
                ],
                output_names=["out_mrtrix_matrix"],
                function=convert_flirt_transformation_to_mrtrix_transformation,
            ),
        )

        # Parc. Transformation (only if space=b0)
        # --------------------
        parc_transform_node = npe.MapNode(
            name="Reg-2-ParcTransformation",
            iterfield=["in_files", "out_filename"],
            interface=MRTransform(),
        )

        # Response Estimation
        # -------------------
        resp_estim_node = npe.Node(
            name="1a-ResponseEstimation", interface=mrtrix3.ResponseSD()
        )
        resp_estim_node.inputs.algorithm = "tournier"

        # FOD Estimation
        # --------------
        fod_estim_node = npe.Node(
            name="1b-FODEstimation",
            interface=ConstrainedSphericalDeconvolution(),
        )
        fod_estim_node.inputs.algorithm = "csd"

        # Tracts Generation
        # -----------------
        tck_gen_node = npe.Node(name="2-TractsGeneration", interface=Tractography())
        tck_gen_node.inputs.select = self.parameters["n_tracks"]
        tck_gen_node.inputs.algorithm = "iFOD2"

        # Connectome Generation
        # ---------------------
        # only the parcellation and output filename should be iterable, the tck
        # file stays the same.
        conn_gen_node = npe.MapNode(
            name="3-ConnectomeGeneration",
            iterfield=["in_parc", "out_file"],
            interface=mrtrix3.BuildConnectome(),
        )

        # Print begin message
        # -------------------
        print_begin_message = npe.MapNode(
            interface=niu.Function(
                input_names=["in_bids_or_caps_file"],
                function=print_begin_pipeline,
            ),
            iterfield="in_bids_or_caps_file",
            name="WriteBeginMessage",
        )

        # Print end message
        # -----------------
        print_end_message = npe.MapNode(
            interface=niu.Function(
                input_names=["in_bids_or_caps_file", "final_file"],
                function=print_end_pipeline,
            ),
            iterfield=["in_bids_or_caps_file"],
            name="WriteEndMessage",
        )

        # CAPS File names Generation
        # --------------------------
        caps_filenames_node = npe.Node(
            name="CAPSFilenamesGeneration",
            interface=niu.Function(
                input_names="dwi_file",
                output_names=self.get_output_fields(),
                function=get_caps_filenames,
            ),
        )
        self.connect(
            [
                (
                    self.input_node,
                    print_begin_message,
                    [("dwi_file", "in_bids_or_caps_file")],
                ),
                (self.input_node, caps_filenames_node, [("dwi_file", "dwi_file")]),
                # Response Estimation
                (
                    self.input_node,
                    resp_estim_node,
                    [("dwi_file", "in_file")],
                ),  # Preproc. DWI
                (
                    self.input_node,
                    resp_estim_node,
                    [("dwi_brainmask_file", "in_mask")],
                ),  # B0 brain mask
                (
                    self.input_node,
                    resp_estim_node,
                    [("grad_fsl", "grad_fsl")],
                ),  # bvecs and bvals
                (
                    caps_filenames_node,
                    resp_estim_node,
                    [("response", "wm_file")],
                ),  # output response filename
                # FOD Estimation
                (
                    self.input_node,
                    fod_estim_node,
                    [("dwi_file", "in_file")],
                ),  # Preproc. DWI
                (
                    resp_estim_node,
                    fod_estim_node,
                    [("wm_file", "wm_txt")],
                ),  # Response (txt file)
                (
                    self.input_node,
                    fod_estim_node,
                    [("dwi_brainmask_file", "mask_file")],
                ),  # B0 brain mask
                (
                    self.input_node,
                    fod_estim_node,
                    [("grad_fsl", "grad_fsl")],
                ),  # T1-to-B0 matrix file
                (
                    caps_filenames_node,
                    fod_estim_node,
                    [("fod", "wm_odf")],
                ),  # output odf filename
                # Tracts Generation
                (fod_estim_node, tck_gen_node, [("wm_odf", "in_file")]),  # ODF file
                (
                    caps_filenames_node,
                    tck_gen_node,
                    [("tracts", "out_file")],
                ),  # output tck filename
                # Label Conversion
                (
                    self.input_node,
                    label_convert_node,
                    [("atlas_files", "in_file")],
                ),  # atlas image files
                (
                    caps_filenames_node,
                    label_convert_node,
                    [("nodes", "out_file")],
                ),  # converted atlas image filenames
                # Connectomes Generation
                (tck_gen_node, conn_gen_node, [("out_file", "in_file")]),
                (caps_filenames_node, conn_gen_node, [("connectomes", "out_file")]),
            ]
        )
        # Registration T1-DWI (only if space=b0)
        # -------------------
        if self.parameters["dwi_space"] == "b0":
            self.connect(
                [
                    # MGZ Files Conversion
                    (
                        self.input_node,
                        t1_brain_conv_node,
                        [("t1_brain_file", "in_file")],
                    ),
                    (self.input_node, wm_mask_conv_node, [("wm_mask_file", "in_file")]),
                    # B0 Extraction
                    (self.input_node, split_node, [("dwi_file", "in_file")]),
                    (split_node, select_node, [("out_files", "inlist")]),
                    # Masking
                    (select_node, mask_node, [("out", "in_file")]),  # B0
                    (
                        self.input_node,
                        mask_node,
                        [("dwi_brainmask_file", "mask_file")],
                    ),  # Brain mask
                    # T1-to-B0 Registration
                    (
                        t1_brain_conv_node,
                        t12b0_reg_node,
                        [("out_file", "in_file")],
                    ),  # Brain
                    (
                        mask_node,
                        t12b0_reg_node,
                        [("out_file", "reference")],
                    ),  # B0 brain-masked
                    # WM Transformation
                    (
                        wm_mask_conv_node,
                        wm_transform_node,
                        [("out_file", "in_file")],
                    ),  # Brain mask
                    (
                        mask_node,
                        wm_transform_node,
                        [("out_file", "reference")],
                    ),  # BO brain-masked
                    (
                        t12b0_reg_node,
                        wm_transform_node,
                        [("out_matrix_file", "in_matrix_file")],
                    ),  # T1-to-B0 matrix file
                    # FSL flirt matrix to MRtrix matrix Conversion
                    (
                        t1_brain_conv_node,
                        fsl2mrtrix_conv_node,
                        [("out_file", "in_source_image")],
                    ),
                    (
                        mask_node,
                        fsl2mrtrix_conv_node,
                        [("out_file", "in_reference_image")],
                    ),
                    (
                        t12b0_reg_node,
                        fsl2mrtrix_conv_node,
                        [("out_matrix_file", "in_flirt_matrix")],
                    ),
                    # Apply registration without resampling on parcellations
                    (
                        label_convert_node,
                        parc_transform_node,
                        [("out_file", "in_files")],
                    ),
                    (
                        fsl2mrtrix_conv_node,
                        parc_transform_node,
                        [("out_mrtrix_matrix", "linear_transform")],
                    ),
                    (
                        caps_filenames_node,
                        parc_transform_node,
                        [("nodes", "out_filename")],
                    ),
                ]
            )
        # Special care for Parcellation & WM mask
        # ---------------------------------------
        if self.parameters["dwi_space"] == "b0":
            self.connect(
                [
                    (wm_transform_node, tck_gen_node, [("out_file", "seed_image")]),
                    (parc_transform_node, conn_gen_node, [("out_file", "in_parc")]),
                    (parc_transform_node, self.output_node, [("out_file", "nodes")]),
                ]
            )
        elif self.parameters["dwi_space"] == "T1w":
            self.connect(
                [
                    (self.input_node, tck_gen_node, [("wm_mask_file", "seed_image")]),
                    (label_convert_node, conn_gen_node, [("out_file", "in_parc")]),
                    (label_convert_node, self.output_node, [("out_file", "nodes")]),
                ]
            )
        else:
            raise ClinicaCAPSError(
                "Bad preprocessed DWI space. Please check your CAPS folder."
            )
        self.connect(
            [
                (resp_estim_node, self.output_node, [("wm_file", "response")]),
                (fod_estim_node, self.output_node, [("wm_odf", "fod")]),
                (tck_gen_node, self.output_node, [("out_file", "tracts")]),
                (conn_gen_node, self.output_node, [("out_file", "connectomes")]),
                (
                    self.input_node,
                    print_end_message,
                    [("dwi_file", "in_bids_or_caps_file")],
                ),
                (conn_gen_node, print_end_message, [("out_file", "final_file")]),
            ]
        )
