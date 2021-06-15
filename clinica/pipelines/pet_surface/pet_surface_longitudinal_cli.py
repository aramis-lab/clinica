# coding: utf-8

import clinica.engine as ce


class PetSurfaceLongitudinalCLI(ce.CmdParser):
    def define_name(self):
        """Define the sub-command name to run this pipeline."""
        self._name = "pet-surface-longitudinal"

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = (
            "Longitudinal surface-based processing of PET images:\n"
            "https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/PET_Surface_Longitudinal/"
        )

    def define_options(self):
        """Define the sub-command arguments."""
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES

        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_label)
        clinica_comp = self._args.add_argument_group(
            PIPELINE_CATEGORIES["CLINICA_COMPULSORY"]
        )
        clinica_comp.add_argument("bids_directory", help="Path to the BIDS directory.")
        clinica_comp.add_argument(
            "caps_directory",
            help="Path to the CAPS directory. (Filled with results from t1-freesurfer-longitudinal pipeline",
        )
        clinica_comp.add_argument(
            "acq_label",
            type=str,
            help="Name of the PET tracer label in the acquisition entity "
            "(acq-<acq_label>).",
        )
        clinica_comp.add_argument(
            "suvr_reference_region",
            choices=["cerebellumPons", "pons"],
            help="Intensity normalization using the average PET uptake in reference regions "
            "resulting in a standardized uptake value ratio (SUVR) map. It can be "
            "cerebellumPons (used for amyloid tracers) or pons (used for 18F-FDG tracers).",
        )
        clinica_comp.add_argument(
            "pvc_psf_tsv",
            help="TSV file containing for each PET image its point spread function (PSF) measured "
            "in mm at x, y & z coordinates. Columns must contain: "
            "participant_id, session_id, acq_label, psf_x, psf_y and psf_z.",
        )
        # Clinica standard arguments (e.g. --n_procs)
        self.add_clinica_standard_arguments()

    def run_command(self, args):
        """Run the pipeline with defined args."""
        from networkx import Graph

        from clinica.utils.ux import print_crash_files_and_exit, print_end_pipeline

        from .pet_surface_pipeline import PetSurface

        parameters = {
            "acq_label": args.acq_label,
            "suvr_reference_region": args.suvr_reference_region,
            "pvc_psf_tsv": self.absolute_path(args.pvc_psf_tsv),
            "longitudinal": True,
        }
        pipeline = PetSurface(
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
