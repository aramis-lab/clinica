# coding: utf-8

"""
 Command line for adni_to_bids converter
"""

import clinica.engine as ce

__author__ = "Jorge Samper Gonzalez and Sabrina Fontanella"
__copyright__ = "Copyright 2017, The Aramis Lab Team"
__credits__ = [""]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


class AdniToBidsCLI(ce.CmdParser):
    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """

        self._name = 'adni-to-bids'

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
        self._args.add_argument("-m", "--modality",
                                help='(Optional) Convert only a selected modality. Modalities available: '
                                      'T1, PET_FDG, PET_AV45.')

    def run_pipeline(self, args):
        from clinica.iotools.converters.adni_to_bids.adni_to_bids import AdniToBids
        adni_to_bids = AdniToBids()

        # Check dcm2nii and dcm2niix dependencies
        adni_to_bids.check_adni_dependencies()

        if not args.clinical_data_only:
            adni_to_bids.convert_images(self.absolute_path(args.dataset_directory),
                                        self.absolute_path(args.clinical_data_directory),
                                        self.absolute_path(args.bids_directory),
                                        args.subjects_list,
                                        args.modality)

        adni_to_bids.convert_clinical_data(self.absolute_path(args.clinical_data_directory),
                                           self.absolute_path(args.bids_directory))
