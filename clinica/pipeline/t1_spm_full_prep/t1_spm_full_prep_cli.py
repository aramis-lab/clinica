"""T1 SPM Full Prep - Clinica Command Line Interface.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details: https://gitlab.icm-institute.org/aramislab/clinica/wikis/docs/InteractingWithClinica.
"""


import clinica.engine as ce


class T1SPMFullPrepCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """
        self._name = 't1-spm-full-prep'

    def define_options(self):
        """Define the sub-command arguments
        """
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')
        self._args.add_argument("caps_directory",
                                help='Path to the CAPS directory.')
        self._args.add_argument("group_id",
                                help='Current group name')
        self._args.add_argument("-tsv", "--subjects_sessions_tsv",
                                help='TSV file containing the subjects with their sessions.')
        self._args.add_argument("-ti", "--tissue_classes", nargs='+', type=int, default=[1, 2, 3], choices=range(1, 7),
                                help="Tissue classes (gray matter, GM; white matter, WM; cerebro-spinal fluid, CSF...) to save. Up to 6 tissue classes can be saved. Ex: 1 2 3 is GM, WM and CSF")
        self._args.add_argument("-dt", "--dartel_tissues", nargs='+', type=int, default=[1, 2, 3], choices=range(1, 7),
                                help='Tissues to use for DARTEL template calculation. Ex: 1 is only GM')
        self._args.add_argument("-swu", "--save_warped_unmodulated", action='store_true',
                                help="Save warped unmodulated images for tissues specified in --tissue_classes")
        self._args.add_argument("-swm", "--save_warped_modulated", action='store_true',
                                help="Save warped modulated images for tissues specified in --tissue_classes")
        self._args.add_argument("-fwhm", "--fwhm", nargs='+', type=int, default=[8],
                                help="A list of integers specifying the different isomorphic fwhm in milimeters to smooth the image")
        self._args.add_argument("-m", "--modulate", type=bool, default=True,
                                help='A boolean. Modulate output images - no modulation preserves concentrations')
        self._args.add_argument("-vs", "--voxel_size", nargs=3, type=float,
                                help="A list of 3 floats specifying voxel sizes for each dimension of output image")
        self._args.add_argument("-atlases", "--atlases", nargs='+', type=str,
                                default=['AAL2', 'LPBA40', 'Neuromorphometrics', 'AICHA', 'Hammers'],
                                choices=['AAL2', 'LPBA40', 'Neuromorphometrics', 'AICHA', 'Hammers'],
                                help='A list of atlases to use to calculate the mean GM concentration at each region')
        self._args.add_argument("-wd", "--working_directory",
                                help='Temporary directory to store pipeline intermediate results')
        self._args.add_argument("-np", "--n_procs", type=int,
                                help='Number of cores used to run in parallel')

    def run_pipeline(self, args):
        """
        """
        from t1_spm_full_prep_pipeline import T1SPMFullPrep

        pipeline = T1SPMFullPrep(bids_directory=self.absolute_path(args.bids_directory),
                                 caps_directory=self.absolute_path(args.caps_directory),
                                 tsv_file=self.absolute_path(args.subjects_sessions_tsv),
                                 group_id=args.group_id
                                 )

        pipeline.parameters.update({'tissue_classes': args.tissue_classes,
                                    'dartel_tissues': args.dartel_tissues,
                                    'save_warped_unmodulated': args.save_warped_unmodulated,
                                    'save_warped_modulated': args.save_warped_modulated,
                                    'write_deformation_fields': [True, True],  # args.write_deformation_fields
                                    'save_t1_mni': True,
                                    'voxel_size': tuple(args.voxel_size) if args.voxel_size is not None else None,
                                    'modulation': args.modulate,
                                    'fwhm': args.fwhm,
                                    'atlas_list': args.atlases
                                    })

        pipeline.base_dir = self.absolute_path(args.working_directory)

        if args.n_procs:
            pipeline.run(plugin='MultiProc', plugin_args={'n_procs': args.n_procs})
        else:
            pipeline.run()
