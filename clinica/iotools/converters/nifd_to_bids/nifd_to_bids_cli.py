# coding: utf8

import clinica.engine as ce


class NifdToBidsCLI(ce.CmdParser):
    """
    todo:add description
    """

    def define_name(self):
        """Define the sub-command name to run this pipelines.
        """
        self._name = 'nifd-to-bids'

    def define_description(self):
        """Define a description of this pipeline.
        """
        self._description = 'Convert NIFD (http://4rtni-ftldni.ini.usc.edu/) into BIDS.'

    def define_options(self):
        """Define the sub-command arguments
        """
        self._args.add_argument("dataset_directory",
                                help='Path to the NIFD images directory.')
        self._args.add_argument("clinical_data_directory",
                                help='Path to the directory containing the NIFD clinical files. '
                                     '(NIFD_Clinical_Data_2017_final_updated.xlsx, '
                                     'DataDictionary_NIFD_2017.10.18.xlsx, idaSearch_1_17_2019_NIFD_all.csv)')
        # self._args.add_argument("ida_file",
        #                         help='Path to the NIFD ida.tsv file, path/to/file/ida.tsv')
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')

    def run_command(self, args):

        from clinica.iotools.converters.nifd_to_bids.nifd_to_bids import convert_images, convert_clinical_data
        from clinica.iotools.converter_utils import check_bin
        from clinica.utils.stream import cprint
        from colorama import Fore
        import sys

        missing_bin = check_bin('dcm2niix')
        if missing_bin:
            cprint(Fore.RED + 'dcm2niix is required.'
                   + ' Install it and re-run the converter. You can use the following link :'
                   + Fore.BLUE + '\nhttps://github.com/rordenlab/dcm2niix'
                   + Fore.RESET)
            cprint('Exiting clinica...')
            sys.exit()

        # to_convert = convert_images(args.dataset_directory, args.ida_file, args.bids_directory)
        # convert_clinical_data(args.bids_directory, args.ida_file, args.clinical_data_file, to_convert)
        to_convert = convert_images(args.dataset_directory, args.bids_directory, args.clinical_data_directory)
        convert_clinical_data(args.bids_directory, args.clinical_data_directory, to_convert)
        cprint(Fore.GREEN + 'Conversion to BIDS succeeded' + Fore.RESET)
