# coding: utf8

import clinica.engine as ce


class PETVolumeCLI(ce.CmdParser):
    def define_name(self):
        """Define the sub-command name to run this pipeline."""
        self._name = "pet-volume"

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = (
            "SPM-based pre-processing of PET images:\n"
            "https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/PET_Volume/"
        )

    def define_options(self):
        """Define the sub-command arguments."""
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES
        from clinica.utils.pet import LIST_SUVR_REFERENCE_REGIONS

        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_label)
        clinica_comp = self._args.add_argument_group(
            PIPELINE_CATEGORIES["CLINICA_COMPULSORY"]
        )
        clinica_comp.add_argument("bids_directory", help="Path to the BIDS directory.")
        clinica_comp.add_argument("caps_directory", help="Path to the CAPS directory.")
        clinica_comp.add_argument(
            "group_label",
            help="User-defined identifier for the provided group of subjects.",
        )
        clinica_comp.add_argument(
            "acq_label",
            type=str,
            help="Name of the label given to the PET acquisition, specifying the tracer used (acq-<acq_label>).",
        )
        clinica_comp.add_argument(
            "suvr_reference_region",
            choices=LIST_SUVR_REFERENCE_REGIONS,
            help="Intensity normalization using the average PET uptake in reference regions "
            "resulting in a standardized uptake value ratio (SUVR) map. It can be "
            "cerebellumPons (used for amyloid tracers) or pons (used for 18F-FDG tracers).",
        )
        # Optional arguments (e.g. FWHM)
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES["OPTIONAL"])
        optional.add_argument(
            "-psf",
            "--pvc_psf_tsv",
            help="TSV file containing for each PET image its point spread function (PSF) measured "
            "in mm at x, y & z coordinates. Columns must contain: "
            "participant_id, session_id, acq_label, psf_x, psf_y and psf_z.",
        )
        # Clinica standard arguments (e.g. --n_procs)
        self.add_clinica_standard_arguments()
        # Advanced arguments (i.e. tricky parameters)
        advanced = self._args.add_argument_group(PIPELINE_CATEGORIES["ADVANCED"])
        advanced.add_argument(
            "-mask",
            "--mask_tissues",
            nargs="+",
            type=int,
            default=[1, 2, 3],
            choices=range(1, 7),
            help="Tissue classes (1: gray matter (GM), 2: white matter (WM), 3: cerebrospinal "
            "fluid (CSF), 4: bone, 5: soft-tissue, 6: background) to use "
            "for masking the PET image (default: GM, WM and CSF i.e. --mask_tissues 1 2 3).",
        )
        advanced.add_argument(
            "-threshold",
            "--mask_threshold",
            type=float,
            default=0.3,
            help="Value used as threshold to binarize the tissue maps (default: --mask_threshold 0.3).",
        )
        advanced.add_argument(
            "-pvc_mask",
            "--pvc_mask_tissues",
            nargs="+",
            type=int,
            default=[1, 2, 3],
            choices=range(1, 7),
            help="Tissue classes (1: gray matter (GM), 2: white matter (WM), 3: cerebrospinal "
            "fluid (CSF), 4: bone, 5: soft-tissue, 6: background) to use "
            "as mask for PVC (default: GM, WM and CSF i.e. --pvc_mask_tissues 1 2 3).",
        )
        advanced.add_argument(
            "-smooth",
            "--smooth",
            nargs="+",
            type=int,
            default=[8],
            metavar="FWHM_N FWHM_M",
            help="A list of integers (e.g. --smooth 6 8) specifying the different isomorphic FWHM "
            "in millimeters to smooth the image (default: --smooth 8).",
        )

    def run_command(self, args):
        """Run the pipeline with defined args."""
        from networkx import Graph

        from clinica.pipelines.pet_volume.pet_volume_pipeline import PETVolume
        from clinica.utils.ux import print_crash_files_and_exit, print_end_pipeline

        parameters = {
            "group_label": args.group_label,
            "acq_label": args.acq_label,
            "suvr_reference_region": args.suvr_reference_region,
            "pvc_psf_tsv": self.absolute_path(args.pvc_psf_tsv),
            "mask_tissues": args.mask_tissues,
            "mask_threshold": args.mask_threshold,
            "pvc_mask_tissues": args.pvc_mask_tissues,
            "smooth": args.smooth,
        }
        pipeline = PETVolume(
            bids_directory=self.absolute_path(args.bids_directory),
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv),
            base_dir=self.absolute_path(args.working_directory),
            parameters=parameters,
            name=self.name,
        )

        if args.n_procs:
            exec_pipeline = pipeline.run(
                plugin="MultiProc", plugin_args={"n_procs": args.n_procs}
            )
        else:
            exec_pipeline = pipeline.run()

        if isinstance(exec_pipeline, Graph):
            print_end_pipeline(
                self.name, pipeline.base_dir, pipeline.base_dir_was_specified
            )
        else:
            print_crash_files_and_exit(args.logname, pipeline.base_dir)
