# coding: utf8


import clinica.engine as ce
from colorama import Fore

class DeepLearningPrepareDataCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline."""
        self._name = 'deeplearning-prepare-data'

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = ('Prepare data generated Clinica for PyTorch with Tensor extraction:\n'
                             'http://clinica.run/doc/Pipelines/DeepLearning_PrepareData/')

    def define_options(self):
        """Define the sub-command arguments."""
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES

        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_label...)
        # Most of the time, you will want to read your pipeline inputs into
        # a BIDS and/or CAPS directory. If your pipeline does not require BIDS input,
        # simply remove the two lines involving the BIDS directory.
        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')
        clinica_comp.add_argument("modality",
                help='''For which modality the tensor will be extracted.
                't1-linear': images prepocessed with t1-linear pipeline.
                't1-extensive': images preprocessed with t1-extensive pipeline.
                'custom': find images with a custom suffix in their filename and 
                transform them to tensor format.''',
                choices=['t1-linear', 't1-extensive', 'custom'], default='t1-linear'
                )                          
        clinica_comp.add_argument("extract_method",
                help='''Format of the extracted features. Three options:
                'image' to convert to PyTorch tensor the complete 3D image,
                'patch' to extract 3D volumetric patches and 
                'slice' to extract 2D slices from the image.
                By default the features are extracted from the cropped image.''',
                choices=['image', 'slice', 'patch'], default='image'
                )

        optional = self._args.add_argument_group(PIPELINE_CATEGORIES['OPTIONAL'])
        
        optional.add_argument('-uui', '--use_uncropped_image',
                help='''Use the uncropped image instead of the
                cropped image generated by t1-linear (option only 
                valid fo t1-linear modality).''',
                default=False, action="store_true"
                )

        optional_patch = self._args.add_argument_group(
                "%sPipeline options if you chose ‘patch’ extraction%s" % (Fore.BLUE, Fore.RESET)
                )
        optional_patch.add_argument(
                '-ps', '--patch_size',
                help='''Patch size (default: --patch_size 50).''',
                type=int, default=50
                )
        optional_patch.add_argument(
                '-ss', '--stride_size',
                help='''Stride size (default: --stride_size 50).''',
                type=int, default=50
                )

        optional_slice = self._args.add_argument_group(
                "%sPipeline options if you chose ‘slice’ extraction%s" % (Fore.BLUE, Fore.RESET)
                )
        optional_slice.add_argument(
                '-sd', '--slice_direction',
                help='''Slice direction. Three options:
                '0' -> Sagittal plane,
                '1' -> Coronal plane or
                '2' -> Axial plane
                (default: sagittal plane i.e. --slice_direction 0)''',
                type=int, default=0
                )
        optional_slice.add_argument(
                '-sm', '--slice_mode',
                help='''Slice mode. Two options: 'rgb' to save the slice in
                three identical channels, ‘single’ to save the slice in a
                single channel (default: --slice_mode rgb).''',
                choices=['rgb', 'single'], default='rgb'
                )

        optional_custom = self._args.add_argument_group(
                "%sPipeline options if you chose ‘custom’ modality%s" % (Fore.BLUE, Fore.RESET)
                )
        optional_custom.add_argument(
                '-cn', '--custom_suffix',
                help='''Custom suffix filename, e.g.:
                'graymatter_space-Ixi549Space_modulated-off_probability.nii.gz', or
                'segm-whitematter_probability.nii.gz'
                ''',
                type=str, default=''
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
            'modality': args.modality,    
            'extract_method': args.extract_method,
            'patch_size': args.patch_size,
            'stride_size': args.stride_size,
            'slice_direction': args.slice_direction,
            'slice_mode': args.slice_mode,
            'use_uncropped_image': args.use_uncropped_image,
            'custom_suffix': args.custom_suffix,
        }

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
