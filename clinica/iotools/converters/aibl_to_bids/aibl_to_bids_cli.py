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

    def run_command(self, args):
        """Run the converter with defined args."""
        from os import makedirs
        from os.path import exists

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

        if not exists(args.bids_directory):
            makedirs(args.bids_directory)

        convert_images(
            args.dataset_directory,
            args.clinical_data_directory,
            args.bids_directory,
            args.overwrite,
        )
        convert_clinical_data(args.bids_directory, args.clinical_data_directory)
