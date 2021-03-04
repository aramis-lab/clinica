# coding: utf8

"""pet_linear - Clinica Command Line Interface.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details:
http://clinica.run/doc/InteractingWithClinica/
"""


import clinica.engine as ce


class PETLinearCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline."""
        self._name = 'pet-linear'

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = ('Affine registration of PET images to the MNI standard space:\n'
                             'http://clinica.run/doc/Pipelines/PET_Linear/')

    def define_options(self):
        """Define the sub-command arguments."""
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES
        from clinica.utils.pet import LIST_SUVR_REFERENCE_REGIONS

        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_label...)
        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("bids_directory",
                                  help='Path to the BIDS directory.')
        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')
        clinica_comp.add_argument("acq_label", type=str,
                                  help='Name of the label given to the PET acquisition, specifying the tracer used (acq-<acq_label>).')
        clinica_comp.add_argument("suvr_reference_region",  choices=LIST_SUVR_REFERENCE_REGIONS,
                                  help='Intensity normalization using the average PET uptake in reference regions '
                                       'resulting in a standardized uptake value ratio (SUVR) map. It can be '
                                       'cerebellumPons (used for amyloid tracers) or pons (used for 18F-FDG tracers).')

        # Clinica optional arguments
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES['OPTIONAL'])
        optional.add_argument("-ui", "--uncropped_image",
                              help='''Do not crop the image with template
                              (cropped image are suggested for using with DL
                              models)''',
                              action='store_true',
                              default=False)

        # Clinica standard arguments (e.g. --n_procs)
        self.add_clinica_standard_arguments()


    def run_command(self, args):
        """Run the pipeline with defined args."""
        import os
        from networkx import Graph
        from colorama import Fore
        from .pet_linear_pipeline import PETLinear
        from clinica.utils.stream import cprint
        from clinica.utils.ux import print_end_pipeline, print_crash_files_and_exit

        parameters = {
            'acq_label': args.acq_label,
            'suvr_reference_region': args.suvr_reference_region
            }

        pipeline = PETLinear(
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
