# coding: utf8

import clinica.engine as ce


class SpatialSVMCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline."""
        self._name = 'machinelearning-prepare-spatial-svm'

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = ('Prepare input data for SVM with spatial and anatomical regularization:\n'
                             'http://clinica.run/doc/MachineLeaning_PrepareSpatialSVM')

    def define_options(self):
        """Define the sub-command arguments."""
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES
        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_id)
        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')
        clinica_comp.add_argument("group_id",
                                  help='User-defined identifier for the provided group of subjects.')
        # Optional arguments
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES['OPTIONAL'])
        optional.add_argument("-it", "--image_type",
                              default='t1',
                              help='Imaging modality. Can be t1 or pet (default: --image_type %(default)s)')
        optional.add_argument("-pt", "--pet_tracer",
                              default='fdg',
                              help='PET tracer. Can be fdg or av45 (default: --pet_tracer %(default)s)')
        # Clinica standard arguments (e.g. --n_procs)
        self.add_clinica_standard_arguments()
        # Advanced arguments (i.e. tricky parameters)
        advanced = self._args.add_argument_group(PIPELINE_CATEGORIES['ADVANCED'])
        advanced.add_argument("-fwhm", "--full_width_half_maximum",
                              type=float, metavar='N', default=4,
                              help='Amount of regularization (in mm). In practice, we found the default value '
                                   '(--full_width_half_maximum %(default)s) to be optimal. We therefore '
                                   'do not recommend to change it unless you have a specific reason to do so.')
        advanced.add_argument("-no_pvc", "--no_pvc",
                              action='store_true', default=False,
                              help="Force the use of non PVC PET data (by default, PVC PET data are used)")

    def run_command(self, args):
        """Run the pipeline with defined args."""
        from networkx import Graph
        from .spatial_svm_pipeline import SpatialSVM
        from clinica.utils.ux import print_end_pipeline, print_crash_files_and_exit

        parameters = {
            'group_id': args.group_id,
            'fwhm': args.fwhm,
            'image_type': args.image_type,
            'pet_tracer': args.pet_tracer,
            'no_pvc': args.no_pvc
        }
        pipeline = SpatialSVM(
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
