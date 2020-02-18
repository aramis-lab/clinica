# coding: utf8

"""Statistics_Volume_Correction - Clinica Command Line Interface.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details:
http://clinica.run/doc/InteractingWithClinica/
"""


import clinica.engine as ce


class StatisticsVolumeCorrectionCLI(ce.CmdParser):

    def define_name(self):
        self._name = 'statistics-volume-correction'

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = ('Brief description:\n'
                             'http://clinica.run/doc/Pipelines/StatisticsVolumeCorrection/')

    def define_options(self):
        """Define the sub-command arguments."""
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES

        # Clinica compulsory arguments

        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])

        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')

        clinica_comp.add_argument("t_map",
                                  help='t-statistics map')

        clinica_comp.add_argument("height_threshold",
                                  help='Height Threshold associated with uncorrected p-value')

        clinica_comp.add_argument("FWEp",
                                  help='FWE p value for peak correction')

        clinica_comp.add_argument("FDRp",
                                  help='FDR p value for peak correction')

        clinica_comp.add_argument("FWEc",
                                  help='FWE cluster minimum size value for cluster correction')

        clinica_comp.add_argument("FDRc",
                                  help='FDR cluster minimum size value for cluster correction')
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES['OPTIONAL'])
        optional.add_argument("-nc", "--n_cuts", default=8,
                              help='Number of cuts along each direction')
        # Clinica standard arguments (e.g. --n_procs)
        self.add_clinica_standard_arguments()

    def run_command(self, args):
        import os
        from networkx import Graph
        from .statistics_volume_correction_pipeline import StatisticsVolumeCorrection
        from clinica.utils.ux import print_end_pipeline, print_crash_files_and_exit

        # Most of the time, you will want to instantiate your pipeline with a
        # BIDS and/or CAPS directory as inputs. If the BIDS directory is not needed
        # for your pipeline, simply remove:
        # bids_directory=self.absolute_path(args.bids_directory),
        pipeline = StatisticsVolumeCorrection(
            caps_directory=self.absolute_path(args.caps_directory),
            base_dir=self.absolute_path(args.working_directory),
            name=self.name
        )
        pipeline.parameters = {
            't_map': args.t_map,
            'height_threshold': args.height_threshold,
            'FWEp': args.FWEp,
            'FDRp': args.FDRp,
            'FWEc': args.FWEp,
            'FDRc': args.FDRp,
            'n_cuts': args.n_cuts
        }

        if args.n_procs:
            exec_pipeline = pipeline.run(plugin='MultiProc',
                                         plugin_args={'n_procs': args.n_procs})
        else:
            exec_pipeline = pipeline.run()

        if isinstance(exec_pipeline, Graph):
            print_end_pipeline(self.name, os.path.join(pipeline.base_dir, self.name))
        else:
            print_crash_files_and_exit(args.logname, pipeline.base_dir)
