# coding: utf8

import clinica.engine as ce


class CmdParserSubjectsSessions(ce.CmdParser):

    def define_name(self):
        self._name = 'create-subjects-visits'

    def define_description(self):
        self._description = 'Create a TSV file containing participants with their sessions'

    def define_options(self):
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS dataset directory.')
        self._args.add_argument("out_tsv",
                                help='Output TSV file containing the participants with their sessions.')  # noqa

    def run_command(self, args):
        import os
        import errno
        from clinica.iotools.utils import data_handling as dt
        from clinica.utils.stream import cprint
        from clinica.utils.inputs import check_bids_folder

        check_bids_folder(args.bids_directory)
        output_directory = os.path.dirname(os.path.abspath(args.out_tsv))
        if not os.path.exists(output_directory):
            try:
                os.makedirs(output_directory)
            except OSError as exc:  # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise
        dt.create_subs_sess_list(args.bids_directory, output_directory, os.path.basename(args.out_tsv))
        cprint("The TSV file was saved to %s" % os.path.abspath(args.out_tsv))


class CmdParserMergeTsv(ce.CmdParser):

    def define_name(self):
        self._name = 'merge-tsv'

    def define_description(self):
        self._description = 'Merge TSV files containing clinical data of a BIDS dataset into a single TSV file.'

    def define_options(self):
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES
        # Clinica compulsory arguments (e.g. BIDS, out_tsv)
        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("bids_directory",
                                  help='Path to the BIDS dataset directory.')
        clinica_comp.add_argument("out_tsv",
                                  help='Path to the output file.')

        # Optional arguments
        iotools_options = self._args.add_argument_group(PIPELINE_CATEGORIES['IOTOOLS_OPTIONS'])
        iotools_options.add_argument('-caps', "--caps_directory", type=str, default=None,
                                     help='(Optional) path to a CAPS directory.')
        iotools_options.add_argument('-p', "--pipelines", nargs="*", type=str, default=None,
                                     help='(Optional) pipelines which will be merged to the TSV file. \n'
                                          'Currently, only t1-volume and pet-volume are supported.\n'
                                          'Default: all pipeline are merged')
        iotools_options.add_argument('-atlas', "--atlas_selection", nargs="*", type=str, default=None,
                                     help='(Optional) atlas that will be merged. \n'
                                          'Default: all atlases are merged')
        iotools_options.add_argument('-pvc', "--pvc_restriction", type=int, default=None,
                                     help='(Optional) indicates restriction on the label [_pvc-rbv]\n'
                                          'Default: all atlases are merged\n'
                                          '0: atlases without the label only are merged\n'
                                          '1: atlases with the label only are merged')
        iotools_options.add_argument('-group', '--group_selection', nargs="*", type=str, default=None,
                                     help='(Optional) groups that will be merged. \n'
                                          'Default: all groups are merged')
        iotools_options.add_argument("-tsv", "--subjects_sessions_tsv",
                                     help='TSV file containing the subjects with their sessions.')

    def run_command(self, args):
        from clinica.iotools.utils import data_handling as dt
        from clinica.utils.inputs import check_bids_folder

        check_bids_folder(args.bids_directory)
        dt.create_merge_file(args.bids_directory, args.out_tsv,
                             caps_dir=args.caps_directory, pipelines=args.pipelines,
                             atlas_selection=args.atlas_selection, pvc_restriction=args.pvc_restriction,
                             tsv_file=args.subjects_sessions_tsv, group_selection=args.group_selection)


class CmdParserMissingModalities(ce.CmdParser):

    def define_name(self):
        self._name = 'check-missing-modalities'

    def define_description(self):
        self._description = 'Check missing modalities in a BIDS directory'

    def define_options(self):
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS dataset directory.')
        self._args.add_argument("out_directory",
                                help='Path to the output directory.')  # noqa
        self._args.add_argument("-op", '--output_prefix',
                                type=str, default='',
                                help='Prefix for the name of output files (default: --output_prefix missing_mods).')  # noqa

    def run_command(self, args):
        from clinica.iotools.utils import data_handling as dt
        from clinica.utils.inputs import check_bids_folder

        check_bids_folder(args.bids_directory)
        dt.compute_missing_mods(args.bids_directory, args.out_directory, args.output_prefix)


class CmdParserCenterNifti(ce.CmdParser):

    def define_name(self):
        self._name = 'center-nifti'

    def define_description(self):
        self._description = 'Center NIFTI of a BIDS directory. Tool mainly used when SPM is not able \nto segment some'\
                            + ' T1w images because the centers of these volumes are not\naligned with the origin of the'\
                            + 'world coordinate system. By default, only\nproblematic images are converted. The rest' \
                            + ' of the images are also copied\nto the new BIDS directory, but left untouched.'

    def define_options(self):
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS dataset directory in which you want to center the NIfTI '
                                     + 'files/images.')
        self._args.add_argument("output_bids_directory",
                                help='Path to the output directory. This is where your BIDS folder with centered '
                                     + 'NIfTI files/images will appear.')
        self._args.add_argument("--modality", '-m',
                                help='List of modalities you want to center the NIfTI images (e.g.  --modality "t1w fdg_pet dwi"). '
                                     'All the files whose names contains one the keywords specified '
                                     ' in this list will be kept (default: --modality "t1w")',
                                default='t1w')
        self._args.add_argument('--center_all_files',
                                help='Use this flag without argument to force the conversion of all NIfTI specified by'
                                     + ' the --modality flag.',
                                action='store_true',
                                dest='center_all_files',
                                default=False)

    def run_command(self, args):
        from colorama import Fore
        from os.path import isdir, abspath, join, isfile
        from os import listdir
        from os import makedirs
        from clinica.iotools.utils.data_handling import center_all_nifti, write_list_of_files
        from clinica.utils.stream import cprint
        import sys
        import time

        # check that output_folder does not exist, or is an empty folder
        if isdir(args.output_bids_directory):
            file_list = [file for file in listdir(args.output_bids_directory) if not file.startswith('.')]
            if len(file_list) > 0:
                error_str = Fore.YELLOW + '[Warning] Some files or directory have been found in ' \
                            + abspath(args.output_bids_directory) + ': \n'
                for f in file_list:
                    error_str += '\t' + f + '\n'
                error_str += 'Do you wish to continue ? (If yes, this may overwrite the files mentioned above).'
                error_str += Fore.RESET
                cprint(error_str)
                while True:
                    cprint('Your answer [yes/no]:')
                    answer = input()
                    if answer.lower() in ['yes', 'no']:
                        break
                    else:
                        cprint(Fore.RED + 'You must answer yes or no' + Fore.RESET)
                if answer.lower() == 'no':
                    cprint(Fore.RED + 'Clinica will now exit...' + Fore.RESET)
                    sys.exit(0)
        else:
            makedirs(args.output_bids_directory)

        split_modality = args.modality.split(' ')
        # Remove empty str in list
        split_modality = [element for element in split_modality if element]

        cprint('Clinica is now centering the requested images.')

        centered_files = center_all_nifti(abspath(args.bids_directory),
                                          abspath(args.output_bids_directory),
                                          split_modality,
                                          center_all_files=args.center_all_files)

        # Write list of created files
        timestamp = time.strftime("%Y%m%d-%H%M%S", time.localtime(time.time()))
        log_file = abspath(join(args.output_bids_directory, 'centered_nifti_list_' + timestamp + '.txt'))
        # If an error happen while creating the file, the function returns Nan
        if not write_list_of_files(centered_files, log_file):
            cprint(Fore.YELLOW + '[Warning] Could not create log file' + Fore.RESET)

        # Final message
        cprint(Fore.GREEN + str(len(centered_files)) + ' NIfTI files/images of BIDS folder:\n\t ' + Fore.BLUE
               + abspath(args.bids_directory) + Fore.GREEN + '\n for the modalities ' + Fore.YELLOW + args.modality
               + Fore.GREEN + ' have been centered in output folder:\n\t' + Fore.BLUE
               + abspath(args.output_bids_directory) + Fore.RESET)
        if isfile(log_file):
            cprint(Fore.GREEN + 'The list of centered NIfTI files is available here : ' + Fore.BLUE + log_file
                   + Fore.RESET)
        cprint('Please note that the rest of the input BIDS folder has also been copied to the output folder.')
