# coding: utf8

import clinica.engine as ce


class SpatialSVMCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """
        self._name = 'machinelearning-prepare-spatial-svm'

    def define_description(self):
        self._description = 'Prepare input data for SVM with spatial and anatomical regularization:\nhttp://clinica.run/doc/MachineLeaning_PrepareSpatialSVM'

    def define_options(self):
        """Define the sub-command arguments
        """
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES
        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_id)
        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')
        clinica_comp.add_argument("group_id",
                                  help='User-defined identifier for the provided group of subjects.')
        # Optional arguments
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES['OPTIONAL'])
        optional.add_argument("-it", "--image_type", default='t1',
                              help='Imaging modality. Can be t1 or pet (default: --image_type t1)')
        optional.add_argument("-pt", "--pet_tracer", default='FDG',
                              help='PET tracer. Can be fdg or av45 (default: --pet_tracer fdg)')
        # Clinica standard arguments (e.g. --n_procs)
        clinica_opt = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_OPTIONAL'])
        clinica_opt.add_argument("-tsv", "--subjects_sessions_tsv",
                                 help='TSV file containing a list of subjects with their sessions.')
        clinica_opt.add_argument("-wd", "--working_directory",
                                 help='Temporary directory to store pipeline intermediate results')
        clinica_opt.add_argument("-np", "--n_procs", type=int,
                                 help='Number of cores used to run in parallel')
        # Advanced arguments (i.e. tricky parameters)
        advanced = self._args.add_argument_group(PIPELINE_CATEGORIES['ADVANCED'])
        advanced.add_argument("-fwhm", "--full_width_half_maximum", type=float, default=4,
                              help='Amount of regularization (in mm). In practice, we found the default value (--full_width_half_maximum 4) to be optimal. We therefore do not recommend to change it unless you have a specific reason to do so.')
        advanced.add_argument("-no_pvc", "--no_pvc",
                              action='store_true', default='False',
                              help="Force the use of non PVC PET data (by default, PVC PET data are used)")

    def run_command(self, args):
        """
        """
        from tempfile import mkdtemp
        from clinica.utils.stream import cprint
        from clinica.pipelines.machine_learning_spatial_svm.spatial_svm_pipeline import SpatialSVM

        pipeline = SpatialSVM(
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv))

        pipeline.parameters = {
            'group_id': args.group_id,
            'fwhm': args.full_width_half_maximum,
            'image_type': args.image_type,
            'pet_type': args.pet_tracer,
            'no_pvc': args.no_pvc
        }
        if args.working_directory is None:
            args.working_directory = mkdtemp()
        pipeline.base_dir = self.absolute_path(args.working_directory)
        if args.n_procs:
            pipeline.run(plugin='MultiProc', plugin_args={'n_procs': args.n_procs})
        else:
            pipeline.run()

        cprint("The " + self._name + " pipeline has completed. You can now delete the working directory (" + args.working_directory + ").")
