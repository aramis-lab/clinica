# coding: utf8

import clinica.engine as ce


class OasisToBidsCLI(ce.CmdParser):
    def define_name(self):
        """Define the sub-command name to run this command."""
        self._name = 'oasis-to-bids'

    def define_description(self):
        """Define a description of this command."""
        self._description = 'Convert OASIS (http://oasis-brains.org/) into BIDS.'

    def define_options(self):
        """Define the sub-command arguments."""
        self._args.add_argument("dataset_directory",
                                help='Path to the OASIS images directory.')
        self._args.add_argument("clinical_data_directory",
                                help='Path to the OASIS clinical data directory.')
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')

    def run_command(self, args):
        """Run the converter with defined args."""
        from clinica.iotools.converters.oasis_to_bids.oasis_to_bids import OasisToBids

        oasis_to_bids = OasisToBids()
        oasis_to_bids.convert_images(args.dataset_directory, args.bids_directory)
        oasis_to_bids.convert_clinical_data(args.clinical_data_directory, args.bids_directory)
