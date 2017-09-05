import clinica.engine as ce


class AiblToBidsCLI(ce.CmdParser):
    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """
        self._name = 'convert-aibl-to-bids'

    def define_options(self):
        """Define the sub-command arguments
        """
        self._args.add_argument("dataset_directory",
                               help='Path to the AIBL images directory.')
        self._args.add_argument("clinical_data_directory",
                                help='Path to the OASIS clinical data directory.')
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')

    def run_pipeline(self, args):
        from clinica.iotools.converters.aibl_to_bids.aibl_to_bids import convert_clinical_data, convert_images

        convert_images(self.dataset_directory, self.clinical_data_directory, self.bids_directory)
        convert_clinical_data(self.bids_directory, self.clinical_data_directory)