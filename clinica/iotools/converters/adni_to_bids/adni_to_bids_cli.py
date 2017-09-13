import clinica.engine as ce

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
                                help='(Optional) Given an already existing ADNI BIDS folder, convert only the clinical '
                                     'data.')
        self._args.add_argument("-sl", "--subjects_list",
                                 help='(Optional) A path to a .txt file containing a list of subject to convert '
                                      '(one for each row).')
        self._args.add_argument("-m", "--modality",
                                 help='(Optional) Convert only a selected modality. Modalities available: '
                                      'T1, PET_FDG, PET_AV45.')


    def run_pipeline(self, args):
        from clinica.iotools.converters.adni_to_bids.adni_to_bids import AdniToBids
        adni_to_bids = AdniToBids()

        # adni_to_bids.convert_clinical_data(args.clinical_data_directory, args.bids_directory)
        adni_to_bids.convert_images(args.dataset_directory, args.clinical_data_directory, args.bids_directory, args.subjects_list, args.modality)

