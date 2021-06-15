# coding: utf8

import clinica.engine as ce


class NifdToBidsCLI(ce.CmdParser):
    def define_name(self):
        """Define the sub-command name to run this command."""
        self._name = "nifd-to-bids"

    def define_description(self):
        """Define a description of this command."""
        self._description = "Convert NIFD (http://4rtni-ftldni.ini.usc.edu/) into BIDS."

    def define_options(self):
        """Define the sub-command arguments."""
        self._args.add_argument(
            "dataset_directory", help="Path to the NIFD images directory."
        )
        self._args.add_argument(
            "clinical_data_directory",
            help="Path to the directory containing the NIFD clinical files. "
            "(NIFD_Clinical_Data_2017_final_updated.xlsx, "
            "DataDictionary_NIFD_2017.10.18.xlsx, idaSearch_1_17_2019_NIFD_all.csv)",
        )
        # self._args.add_argument("ida_file",
        #                         help='Path to the NIFD ida.tsv file, path/to/file/ida.tsv')
        self._args.add_argument("bids_directory", help="Path to the BIDS directory.")

    def run_command(self, args):
        """Run the converter with defined args."""

        from clinica.iotools.converters.nifd_to_bids.nifd_to_bids import (
            convert_clinical_data,
            convert_images,
        )
        from clinica.utils.check_dependency import check_dcm2niix
        from clinica.utils.stream import cprint

        check_dcm2niix()

        # to_convert = convert_images(args.dataset_directory, args.ida_file, args.bids_directory)
        # convert_clinical_data(args.bids_directory, args.ida_file, args.clinical_data_file, to_convert)
        to_convert = convert_images(
            args.dataset_directory, args.bids_directory, args.clinical_data_directory
        )
        convert_clinical_data(
            args.bids_directory, args.clinical_data_directory, to_convert
        )
        cprint(msg="Conversion to BIDS succeeded", lvl="info")
