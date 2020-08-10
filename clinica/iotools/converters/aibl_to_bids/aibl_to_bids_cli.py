# coding: utf8

import clinica.engine as ce


class AiblToBidsCLI(ce.CmdParser):
    def define_name(self):
        """Define the sub-command name to run this pipeline."""
        self._name = 'aibl-to-bids'

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = 'Convert AIBL (https://aibl.csiro.au/adni/index.html) into BIDS.'

    def define_options(self):
        """Define the sub-command arguments."""
        self._args.add_argument("dataset_directory",
                                help='Path to the AIBL images directory.')
        self._args.add_argument("clinical_data_directory",
                                help='Path to the AIBL clinical data directory.')
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')

    def run_command(self, args):
        from os.path import exists
        from os import makedirs
        from colorama import Fore
        from clinica.iotools.converters.aibl_to_bids.aibl_to_bids import convert_clinical_data, convert_images
        from clinica.utils.check_dependency import is_binary_present, check_freesurfer
        from clinica.utils.exceptions import ClinicaMissingDependencyError

        list_binaries = ['dcm2niix', 'dcm2nii']
        for binary in list_binaries:
            if not is_binary_present(binary):
                raise ClinicaMissingDependencyError(
                    '%s\n[Error] Clinica could not find %s software: it is not present in your '
                    'PATH environment.%s' % (Fore.RED, binary, Fore.RESET))
        check_freesurfer()

        if not exists(args.bids_directory):
            makedirs(args.bids_directory)

        convert_images(args.dataset_directory, args.clinical_data_directory, args.bids_directory)
        convert_clinical_data(args.bids_directory, args.clinical_data_directory)
