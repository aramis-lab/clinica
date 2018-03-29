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
        self._args.add_argument("out_directory",
                                help='Path to the output directory.')  # noqa
        self._args.add_argument("-on", '--output_name',
                                type=str, default='',
                                help='(Optional) Name of the output file.')  # noqa

    def run_command(self, args):
        from clinica.iotools.utils import data_handling as dt
        dt.create_subs_sess_list(args.bids_directory, args.out_directory, args.output_name)


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
                                     help='(Optional) pipelines which will be merged to the .tsv file. \n'
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
                                help='Prefix for the name of output files.')  # noqa

    def run_command(self, args):
        from clinica.iotools.utils import data_handling as dt
        dt.compute_missing_mods(args.bids_directory, args.out_directory, args.output_prefix)
