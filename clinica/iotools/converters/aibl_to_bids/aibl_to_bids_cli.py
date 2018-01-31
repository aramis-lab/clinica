# coding: utf8

import clinica.engine as ce


class AiblToBidsCLI(ce.CmdParser):
    """
    todo:add description
    """

    def define_name(self):
        """Define the sub-command name to run this pipelines.
        """
        self._name = 'aibl-to-bids'

    def define_description(self):
        """Define a description of this pipeline.
        """
        self._description = 'Convert AIBL (http://adni.loni.usc.edu/) into BIDS'

    def define_options(self):
        """Define the sub-command arguments
        """
        self._args.add_argument("dataset_directory",
                               help='Path to the AIBL images directory.')
        self._args.add_argument("clinical_data_directory",
                                help='Path to the AIBL clinical data directory.')
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')

    def run_command(self, args):
        from clinica.iotools.converters.aibl_to_bids.aibl_to_bids import convert_clinical_data, convert_images

        convert_images(args.dataset_directory, args.clinical_data_directory, args.bids_directory)
        convert_clinical_data(args.bids_directory, args.clinical_data_directory)
