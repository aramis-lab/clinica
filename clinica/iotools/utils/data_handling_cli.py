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
        output_directory = os.path.dirname(args.out_tsv)
        if not os.path.exists(output_directory):
            try:
                os.makedirs(output_directory)
            except OSError as exc:  # Guard against race condition
                if exc.errno != errno.EEXIST:
                    raise
        dt.create_subs_sess_list(args.bids_directory, output_directory, os.path.basename(args.out_tsv))


class CmdParserMergeTsv(ce.CmdParser):

    def define_name(self):
        self._name = 'merge-tsv'

    def define_description(self):
        self._description = 'Merge TSV files containing clinical data of a BIDS dataset into a single TSV file.'

    def define_options(self):
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS dataset directory.')
        self._args.add_argument("out_directory",
                                help='Path to the output directory.')
        self._args.add_argument("-tf", '--true_false_mode', type=bool, default=False,
                                help='(Optional) Convert all the field with binary meaning into True and False values.')

    def run_command(self, args):
        from clinica.iotools.utils import data_handling as dt
        dt.create_merge_file(args.bids_directory, args.out_directory, args.true_false_mode)


class CmdParserMissingModalities(ce.CmdParser):

    def define_name(self):
        self._name = 'check-missing-modalities'

    def define_description(self):
        self._description = 'Check missing modalities in a BIDS directory.'

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
