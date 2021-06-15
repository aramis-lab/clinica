# coding: utf8

import clinica.engine as ce


class DwiPreprocessingUsingT1Cli(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline."""
        self._name = 'dwi-preprocessing-using-t1'

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = ('Preprocessing of raw DWI datasets using a T1w image:\n'
                             'https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/DWI_Preprocessing/')

    def define_options(self):
        """Define the sub-command arguments."""
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES

        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_label)
        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("bids_directory",
                                  help='Path to the BIDS directory.')
        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')
        # Optional arguments (e.g. FWHM)
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES['OPTIONAL'])
        optional.add_argument("--low_bval",
                              metavar='N', type=int, default=5,
                              help='Define the b0 volumes as all volume bval <= low_bval '
                                   '(default: --low_bval %(default)s).')
        # Clinica standard arguments (e.g. --n_procs)
        self.add_clinica_standard_arguments()

    def run_command(self, args):
        """Run the pipeline with defined args."""
        from networkx import Graph

        from clinica.utils.ux import print_crash_files_and_exit, print_end_pipeline

        from .dwi_preprocessing_using_t1_pipeline import DwiPreprocessingUsingT1

        parameters = {
            'low_bval': args.low_bval
        }
        pipeline = DwiPreprocessingUsingT1(
            bids_directory=self.absolute_path(args.bids_directory),
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv),
            base_dir=self.absolute_path(args.working_directory),
            parameters=parameters,
            name=self.name
        )

        if args.n_procs:
            exec_pipeline = pipeline.run(plugin='MultiProc',
                                         plugin_args={'n_procs': args.n_procs})
        else:
            exec_pipeline = pipeline.run()

        if isinstance(exec_pipeline, Graph):
            print_end_pipeline(self.name, pipeline.base_dir, pipeline.base_dir_was_specified)
        else:
            print_crash_files_and_exit(args.logname, pipeline.base_dir)
