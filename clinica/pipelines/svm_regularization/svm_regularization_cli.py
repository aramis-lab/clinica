# coding: utf8

"""svm_regularization - Clinica Command Line Interface.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details:
http://clinica.run/doc/InteractingWithClinica/
"""


import clinica.engine as ce


class SVMRegularizationCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """

        self._name = 'machinelearning-prepare-regularized-svm'

    def define_description(self):
        self._description = 'Spatial and anatomical regularization for SVM\n'  # link to the doc

    def define_options(self):
        """Define the sub-command arguments
        """
        # @todo mettere i parametri di default nell'help

        from clinica.engine.cmdparser import PIPELINE_CATEGORIES
        # self._args.add_argument("bids_directory",
        #                         help='Path to the BIDS directory.')
        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])

        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')
        clinica_comp.add_argument("group_id",
                                  help='User-defined identifier for the provided group of subjects.')

        optional = self._args.add_argument_group(PIPELINE_CATEGORIES['OPTIONAL'])

        optional.add_argument("-fwhm", "--full_width_half_maximum", type=float, default=4,
                              help='fwhm value for regularization (in mm )')
        optional.add_argument("-vs", "--voxel_size", type=float, default=1.5,
                              help='voxel size (in mm)')
        optional.add_argument("-image_type", "--image_type", default='t1',
                              help='Possible values: t1/pet')
        optional.add_argument("-pet_type", "--pet_type", default='FDG',
                              help='Type of pet used, ex: fdg/av45')

        clinica_opt = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_OPTIONAL'])
        clinica_opt.add_argument("-tsv", "--subjects_sessions_tsv",
                                 help='TSV file containing a list of subjects with their sessions.')
        clinica_opt.add_argument("-wd", "--working_directory",
                                 help='Temporary directory to store pipeline intermediate results')
        clinica_opt.add_argument("-np", "--n_procs", type=int,
                                 help='Number of cores used to run in parallel')
        clinica_opt.add_argument("-sl", "--slurm", action='store_true',
                                 help='Run the pipeline using SLURM')

    def run_command(self, args):
        """
        """

        from tempfile import mkdtemp
        from clinica.pipelines.svm_regularization.svm_regularization_pipeline import SVMRegularization

        pipeline = SVMRegularization(
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv))

        pipeline.parameters = {
            'group_id': args.group_id,
            'h': args.voxel_size,
            'fwhm': args.full_width_half_maximum,
            'image_type': args.image_type,
            'pet_type': args.pet_type,
        }
        if args.working_directory is None:
            args.working_directory = mkdtemp()
        pipeline.base_dir = self.absolute_path(args.working_directory)
        if args.n_procs:
            pipeline.run(plugin='MultiProc', plugin_args={'n_procs': args.n_procs})
        elif args.slurm:
            pipeline.run(plugin='SLURM')
        else:
            pipeline.run()
