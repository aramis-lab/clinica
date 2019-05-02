# coding: utf-8

"""
 Command line for adni_to_bids converter
"""

import clinica.engine as ce

__author__ = "Jorge Samper-Gonzalez and Sabrina Fontanella"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper-Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


class AdniToBidsCLI(ce.CmdParser):
    def define_name(self):
        """Define the sub-command name to run this pipelines.
        """
        self._name = 'adni-to-bids'

    def define_description(self):
        """Define a description of this pipeline.
        """
        self._description = 'Convert ADNI (http://adni.loni.usc.edu/) into BIDS'

    def define_options(self):
        """Define the sub-command arguments
        """
        self._args.add_argument("dataset_directory",
                                help='Path to the ADNI images directory.')
        self._args.add_argument("clinical_data_directory",
                                help='Path to the ADNI clinical data directory.')
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')
        self._args.add_argument("-c", "--clinical_data_only", action='store_true',
                                help='(Optional) Given the path to an already existing ADNI BIDS folder, convert only '
                                     'the clinical data.')
        self._args.add_argument("-sl", "--subjects_list",
                                help='(Optional) A path to a .txt file containing a list of subject to convert '
                                     '(one for each row).')
        self._args.add_argument("-m", "--modalities", nargs='+', default=['T1', 'PET_FDG', 'PET_AV45', 'DWI', 'FLAIR', 'fMRI'],
                                choices=['T1', 'PET_FDG', 'PET_AV45', 'DWI', 'fMRI'],
                                help='(Optional) Convert only the list of selected modalities. '
                                     'By default all modalities are converted. Modalities available: '
                                      'T1, PET_FDG, PET_AV45, DWI, fMRI, FLAIR.')

    def run_command(self, args):
        from clinica.iotools.converters.adni_to_bids.adni_to_bids import AdniToBids
        from os.path import exists, join
        from clinica.utils.stream import cprint
        from colorama import Fore
        from time import sleep

        adni_to_bids = AdniToBids()

        # Check dcm2nii and dcm2niix dependencies
        adni_to_bids.check_adni_dependencies()

        if 'fMRI' in args.modalities and not exists(join(args.clinical_data_directory, 'IDA_MR_Metadata_Listing.csv')):
            cprint(Fore.RED + '[Warning] We could not detect file IDA_MR_Metadata_Listing.csv in your clinical data folder.'
                   + ' Therefore, fMRI conversion is disabled.' + Fore.GREEN
                   + '\nADNI conversion will start in a few seconds...' + Fore.RESET)
            # Enough time for the user to read the message
            sleep(7)
            args.modalities.remove('fMRI')

        if not args.clinical_data_only:
            adni_to_bids.convert_images(self.absolute_path(args.dataset_directory),
                                        self.absolute_path(args.clinical_data_directory),
                                        self.absolute_path(args.bids_directory),
                                        args.subjects_list,
                                        args.modalities)

        adni_to_bids.convert_clinical_data(self.absolute_path(args.clinical_data_directory),
                                           self.absolute_path(args.bids_directory))
