import clinica.pipelines.engine as cpe


class StatisticsVolumeCorrection(cpe.Pipeline):
    """StatisticsVolumeCorrection - Statistical correction of Statistical correction of StatisticsVolume pipeline.

    Returns:
        A clinica pipeline object containing the StatisticsVolumeCorrection pipeline.
    """

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file."""

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Returns:
            A list of (string) input fields name.
        """
        return ["t_map"]

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """
        return []

    def build_input_node(self):
        """Build and connect an input node to the pipeline."""
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from clinica.utils.inputs import clinica_group_reader

        t_map = clinica_group_reader(
            self.caps_directory,
            {
                "pattern": self.parameters["t_map"] + "*",
                "description": "statistics t map",
                "needed_pipeline": "statistics-volume",
            },
        )

        read_parameters_node = npe.Node(
            name="LoadingCLIArguments",
            interface=nutil.IdentityInterface(
                fields=self.get_input_fields(), mandatory_inputs=True
            ),
        )
        read_parameters_node.inputs.t_map = t_map

        self.connect([(read_parameters_node, self.input_node, [("t_map", "t_map")])])

    def build_output_node(self):
        """Build and connect an output node to the pipeline."""

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline."""
        from os.path import abspath, dirname, exists, join, pardir

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import numpy as np

        import clinica.pipelines.statistics_volume_correction.statistics_volume_correction_utils as utils
        from clinica.utils.inputs import RemoteFileStructure, fetch_file
        from clinica.utils.stream import cprint

        peak_correction_FWE = npe.Node(
            name="peak_correction_FWE",
            interface=nutil.Function(
                input_names=["t_map", "t_threshold"],
                output_names=["output"],
                function=utils.peak_correction,
            ),
        )
        peak_correction_FWE.inputs.t_threshold = self.parameters["FWEp"]

        peak_correction_FDR = peak_correction_FWE.clone(name="peak_correction_FDR")
        peak_correction_FDR.inputs.t_threshold = self.parameters["FDRp"]

        cluster_correction_FWE = npe.Node(
            name="cluster_correction_FWE",
            interface=nutil.Function(
                input_names=["t_map", "t_thresh", "c_thresh"],
                output_names=["output"],
                function=utils.cluster_correction,
            ),
        )
        cluster_correction_FWE.inputs.t_thresh = self.parameters["height_threshold"]
        cluster_correction_FWE.inputs.c_thresh = self.parameters["FWEc"]

        cluster_correction_FDR = cluster_correction_FWE.clone(
            name="cluster_correction_FDR"
        )
        cluster_correction_FDR.inputs.t_thresh = self.parameters["height_threshold"]
        cluster_correction_FDR.inputs.c_thresh = self.parameters["FDRc"]

        produce_fig_FWE_peak_correction = npe.Node(
            name="produce_figure_FWE_peak_correction",
            interface=nutil.Function(
                input_names=[
                    "nii_file",
                    "template",
                    "type_of_correction",
                    "t_thresh",
                    "c_thresh",
                    "n_cuts",
                ],
                output_names=["figs"],
                function=utils.produce_figures,
            ),
        )
        produce_fig_FWE_peak_correction.inputs.n_cuts = self.parameters["n_cuts"]

        root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
        path_to_mask = join(root, "resources", "masks")
        url_aramis = "https://aramislab.paris.inria.fr/files/data/img_t1_linear/"
        FILE1 = RemoteFileStructure(
            filename="mni_icbm152_t1_tal_nlin_sym_09a.nii.gz",
            url=url_aramis,
            checksum="3b244ee7e287319d36a25263744c468ef0ab2fe5a94b15a2138844db73b49adf",
        )

        produce_fig_FWE_peak_correction.inputs.template = join(
            path_to_mask, FILE1.filename
        )
        if not (exists(produce_fig_FWE_peak_correction.inputs.template)):
            try:
                fetch_file(FILE1, path_to_mask)
            except IOError as err:
                cprint(
                    msg=f"Unable to download required template (mni_icbm152) for processing: {err}",
                    lvl="error",
                )

        produce_fig_FDR_peak_correction = produce_fig_FWE_peak_correction.clone(
            name="produce_figure_FDR_peak_correction"
        )
        produce_fig_FWE_cluster_correction = produce_fig_FWE_peak_correction.clone(
            name="produce_figure_FWE_cluster_correction"
        )
        produce_fig_FDR_cluster_correction = produce_fig_FWE_peak_correction.clone(
            name="produce_figure_FDR_cluster_correction"
        )

        # fmt: off
        produce_fig_FWE_peak_correction.inputs.type_of_correction = "FWE"
        produce_fig_FDR_peak_correction.inputs.type_of_correction = "FDR"
        produce_fig_FWE_cluster_correction.inputs.type_of_correction = "FWE"
        produce_fig_FDR_cluster_correction.inputs.type_of_correction = "FDR"

        produce_fig_FWE_peak_correction.inputs.t_thresh = self.parameters["FWEp"]
        produce_fig_FDR_peak_correction.inputs.t_thresh = self.parameters["FDRp"]
        produce_fig_FWE_cluster_correction.inputs.t_thresh = self.parameters["height_threshold"]
        produce_fig_FDR_cluster_correction.inputs.t_thresh = self.parameters["height_threshold"]

        produce_fig_FWE_peak_correction.inputs.c_thresh = np.nan
        produce_fig_FDR_peak_correction.inputs.c_thresh = np.nan
        produce_fig_FWE_cluster_correction.inputs.c_thresh = self.parameters["FWEc"]
        produce_fig_FDR_cluster_correction.inputs.c_thresh = self.parameters["FDRc"]
        # fmt: on

        save_fig_peak_correction_FWE = npe.Node(
            name="save_figure_peak_correction_FWE",
            interface=nutil.Function(
                input_names=["t_map", "figs", "name"],
                output_names=[],
                function=utils.generate_output,
            ),
        )
        save_fig_peak_correction_FWE.inputs.name = "FWEp"

        save_fig_peak_correction_FDR = save_fig_peak_correction_FWE.clone(
            name="save_fig_peak_correction_FDR"
        )
        save_fig_peak_correction_FDR.inputs.name = "FDRp"

        save_fig_cluster_correction_FWE = save_fig_peak_correction_FWE.clone(
            name="save_fig_cluster_correction_FWE"
        )
        save_fig_cluster_correction_FWE.inputs.name = "FWEc"

        save_fig_cluster_correction_FDR = save_fig_peak_correction_FWE.clone(
            name="save_fig_cluster_correction_FDR"
        )
        save_fig_cluster_correction_FDR.inputs.name = "FDRc"

        # Connection
        # ==========
        # fmt: off
        self.connect(
            [
                (self.input_node, peak_correction_FWE, [("t_map", "t_map")]),
                (self.input_node, peak_correction_FDR, [("t_map", "t_map")]),
                (self.input_node, cluster_correction_FWE, [("t_map", "t_map")]),
                (self.input_node, cluster_correction_FDR, [("t_map", "t_map")]),
                (peak_correction_FWE, produce_fig_FWE_peak_correction, [("output", "nii_file")]),
                (peak_correction_FDR, produce_fig_FDR_peak_correction, [("output", "nii_file")]),
                (cluster_correction_FWE, produce_fig_FWE_cluster_correction, [("output", "nii_file")]),
                (cluster_correction_FDR, produce_fig_FDR_cluster_correction, [("output", "nii_file")]),
                (produce_fig_FWE_peak_correction, save_fig_peak_correction_FWE, [("figs", "figs")]),
                (produce_fig_FDR_peak_correction, save_fig_peak_correction_FDR, [("figs", "figs")]),
                (produce_fig_FWE_cluster_correction, save_fig_cluster_correction_FWE, [("figs", "figs")]),
                (produce_fig_FDR_cluster_correction, save_fig_cluster_correction_FDR, [("figs", "figs")]),
                (self.input_node, save_fig_peak_correction_FWE, [("t_map", "t_map")]),
                (self.input_node, save_fig_peak_correction_FDR, [("t_map", "t_map")]),
                (self.input_node, save_fig_cluster_correction_FWE, [("t_map", "t_map")]),
                (self.input_node, save_fig_cluster_correction_FDR, [("t_map", "t_map")]),
            ]
        )
