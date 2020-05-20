# coding: utf8


import clinica.engine as ce


class DeepLearningPrepareDataCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline."""
        self._name = 'deeplearning-prepare-data'

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = ('Tensor extraction (Pytorch) from T1w images:\n'
                             'http://clinica.run/doc/Pipelines/DeepLearning_PrepareData/')

    def define_options(self):
        """Define the sub-command arguments."""
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES

        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_id...)
        # Most of the time, you will want to read your pipeline inputs into
        # a BIDS and/or CAPS directory. If your pipeline does not require BIDS input,
        # simply remove the two lines involving the BIDS directory.
        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')
        clinica_comp.add_argument("tsv_file",
                                  help='Path to the tsv file.')
        clinica_comp.add_argument("extract_method",
                                  help='''Method used to extract features. Three options:
                                       'slice' to get 2D slices from the MRI,
                                       'patch' to get 3D volumetric patches or
                                       'whole' to get the complete MRI.''',
                                  choices=['slice', 'patch', 'whole'], default='whole'
                                  )

        # group_id can be used by certain pipelines when some operations are performed at the group level
        # (for example, generation of a template in pipeline t1-volume)
        # clinica_comp.add_argument("group_id",
        #                           help='User-defined identifier for the provided group of subjects.')

        # Add your own pipeline command line arguments here to be used in the
        # method below. Example below:
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES['OPTIONAL'])

        optional.add_argument(
                '-psz', '--patch_size',
                help='''Patch size (only for 'patch' extraction) e.g: --patch_size 50''',
                type=int, default=50
                )
        optional.add_argument(
                '-ssz', '--stride_size',
                help='''Stride size (only for 'patch' extraction) e.g.: --stride_size 50''',
                type=int, default=50
                )
        optional.add_argument(
                '-sd', '--slice_direction',
                help='''Slice direction (only for 'slice' extraction). Three options:
                '0' -> Sagittal plane,
                '1' -> Coronal plane or
                '2' -> Axial plane''',
                type=int, default=0
                )
        optional.add_argument(
                '-sm', '--slice_mode',
                help='''Slice mode (only for 'slice' extraction). Two options:
                'original' to save one single channel (intensity),
                'rgb' to saves three channel (with same intensity).''',
                choices=['original', 'rgb'], default='rgb'
                )

        # Clinica standard arguments (e.g. --n_procs)
        self.add_clinica_standard_arguments()


    def run_command(self, args):
        """Run the pipeline with defined args."""
        from networkx import Graph
        from .deeplearning_prepare_data_pipeline import DeepLearningPrepareData
        from clinica.utils.ux import print_end_pipeline, print_crash_files_and_exit

        parameters = {
            # Add your own pipeline parameters here to use them inside your
            # pipeline. See the file `deeplearning_prepare_data_pipeline.py` to
            # see an example of use.
            'extract_method': args.extract_method,
            'patch_size': args.patch_size,
            'stride_size': args.stride_size,
            'slice_direction': args.slice_direction,
            'slice_mode': args.slice_mode,
        }

        # Most of the time, you will want to instantiate your pipeline with a
        # BIDS and/or CAPS directory as inputs. If the BIDS directory is not needed
        # for your pipeline, simply remove:
        # bids_directory=self.absolute_path(args.bids_directory),
        pipeline = DeepLearningPrepareData(
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
