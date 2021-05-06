# coding: utf8

import clinica.engine as ce


class AiblToBidsCLI(ce.CmdParser):
    def define_name(self):
        """Define the sub-command name to run this command."""
        self._name = "aibl-to-bids"

    def define_description(self):
        """Define a description of this command."""
        self._description = (
            "Convert AIBL (https://aibl.csiro.au/adni/index.html) into BIDS."
        )

    def define_options(self):
        """Define the sub-command arguments."""
        self._args.add_argument(
            "dataset_directory", help="Path to the AIBL images directory."
        )
        self._args.add_argument(
            "clinical_data_directory", help="Path to the AIBL clinical data directory."
        )
        self._args.add_argument("bids_directory", help="Path to the BIDS directory.")
        self._args.add_argument(
            "--overwrite",
            action="store_true",
            default=False,
            help="Overwrites previously written nifti and json files.",
        )
        self._args.add_argument(
            "-c",
            "--clinical_data_only",
            action="store_true",
            help="(Optional) Given the path to an already existing ADNI BIDS folder, convert only "
            "the clinical data. Mutually exclusive with --force_new_extraction",
        )

    def run_command(self, args):
        """Run the converter with defined args."""
        from os import makedirs

        from clinica.iotools.converters.aibl_to_bids.aibl_to_bids import (
            convert_clinical_data,
            convert_images,
        )
        from clinica.utils.check_dependency import (
            check_dcm2nii,
            check_dcm2niix,
            check_freesurfer,
        )

        check_dcm2nii()
        check_dcm2niix()
        check_freesurfer()

        makedirs(args.bids_directory, exist_ok=True)

        if not args.clinical_data_only:
            convert_images(
                args.dataset_directory,
                args.clinical_data_directory,
                args.bids_directory,
                args.overwrite,
            )
        convert_clinical_data(args.bids_directory, args.clinical_data_directory)
