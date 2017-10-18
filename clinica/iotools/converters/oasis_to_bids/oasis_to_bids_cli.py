import clinica.engine as ce

__author__ = "Sabrina Fontanella"
__copyright__ = "Copyright 2017, The Aramis Lab Team"
__credits__ = [""]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Simona Bottani"
__email__ = "simona.bottani@icm-institute.org"
__status__ = "Completed"

class OasisToBidsCLI(ce.CmdParser):
    def define_name(self):
        """Define the sub-command name to run this pipelines.
        """
        self._name = 'oasis-to-bids'

    def define_options(self):
        """Define the sub-command arguments
        """
        self._args.add_argument("dataset_directory",
                               help='Path to the OASIS images directory.')
        self._args.add_argument("clinical_data_directory",
                                help='Path to the OASIS clinical data directory.')
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')

    def run_pipeline(self, args):
        from clinica.iotools.converters.oasis_to_bids.oasis_to_bids import OasisToBids
        oasis_to_bids = OasisToBids()

        oasis_to_bids.convert_images(args.dataset_directory, args.bids_directory)
        oasis_to_bids.convert_clinical_data(args.clinical_data_directory, args.bids_directory)