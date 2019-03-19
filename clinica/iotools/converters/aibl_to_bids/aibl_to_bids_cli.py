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
        self._description = 'Convert AIBL (https://aibl.csiro.au/adni/index.html”) into BIDS.'

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
        from clinica.utils.stream import cprint
        from os.path import exists
        from os import makedirs
        from colorama import Fore
        from clinica.iotools.converter_utils import check_bin
        import sys

        # Check existence of binaries dcm2nii, dcm2niix and mri_convert
        # If they are not found we warn the user and tell him that the bin
        # were not foud. He then has the possibility to run the converter anyway
        missing_bin = 0
        bin_to_test = ['dcm2nii', 'dcm2niix', 'mri_convert']
        for binary in bin_to_test:
            missing_bin = missing_bin + check_bin(binary)
        if missing_bin > 0:
            cprint(Fore.RED + str(missing_bin) + ' binary(es) is (are) missing. '
                   + 'Most important are : dcm2nii and dcm2niix.' + Fore.RESET)
            while True:
                cprint('Do you still want to run the converter ? (yes/no): ')
                answer = input('')
                if answer.lower() in ['yes', 'no']:
                    break
                else:
                    cprint('Possible answers are yes or no.\n')
            if answer.lower() == 'yes':
                cprint('Running the pipeline anyway.')
            if answer.lower() == 'no':
                cprint('Exiting clinica...')
                sys.exit()

        if not exists(args.bids_directory):
            makedirs(args.bids_directory)

        convert_images(args.dataset_directory, args.clinical_data_directory, args.bids_directory)
        convert_clinical_data(args.bids_directory, args.clinical_data_directory)
