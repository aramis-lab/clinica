# coding: utf8

import clinica.engine as ce


class T1VolumeTissueSegmentationCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipelines.
        """
        self._name = 't1-volume-tissue-segmentation'

    def define_description(self):
        """Define a description of this pipeline.
        """
        self._description = 'Tissue segmentation, bias correction and spatial normalization to MNI space' \
                            + ' of T1w images with SPM:\nhttp://clinica.run/doc/Pipelines/T1_Volume/'

    def define_options(self):
        """Define the sub-command arguments
        """
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES
        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_id)
        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("bids_directory",
                                  help='Path to the BIDS directory.')
        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')
        # Optional arguments (e.g. FWHM)
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES['OPTIONAL'])
        optional.add_argument("-s", "--smooth",
                              nargs='+', type=int, default=[8],
                              help="A list of integers specifying the different isomorphic FWHM in millimeters "
                                   "to smooth the image (default: --smooth 8).")
        # Clinica standard arguments (e.g. --n_procs)
        clinica_opt = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_OPTIONAL'])
        clinica_opt.add_argument("-tsv", "--subjects_sessions_tsv",
                                 help='TSV file containing a list of subjects with their sessions.')
        clinica_opt.add_argument("-wd", "--working_directory",
                                 help='Temporary directory to store pipelines intermediate results')
        clinica_opt.add_argument("-np", "--n_procs",
                                 metavar=('N'), type=int,
                                 help='Number of cores used to run in parallel')
        # Advanced arguments (i.e. tricky parameters)
        advanced = self._args.add_argument_group(PIPELINE_CATEGORIES['ADVANCED'])
        advanced.add_argument("-t", "--tissue_classes",
                              metavar='', nargs='+', type=int, default=[1, 2, 3], choices=range(1, 7),
                              help="Tissue classes (1: gray matter (GM), 2: white matter (WM), "
                                   "3: cerebrospinal fluid (CSF), 4: bone, 5: soft-tissue, 6: background) to save "
                                   "(default: GM, WM and CSF i.e. --tissue_classes 1 2 3).")
        advanced.add_argument("-dt", "--dartel_tissues",
                              metavar='', nargs='+', type=int, default=[1, 2, 3], choices=range(1, 7),
                              help='Tissues to use for DARTEL template calculation '
                                   '(default: GM, WM and CSF i.e. --dartel_tissues 1 2 3).')
        advanced.add_argument("-tpm", "--tissue_probability_maps",
                              metavar=('TissueProbabilityMap.nii'),
                              help='Tissue probability maps to use for segmentation (default: TPM from SPM software).')
        advanced.add_argument("-swu", "--save_warped_unmodulated",
                              action='store_true', default=True,
                              help="Save warped unmodulated images for tissues specified in --tissue_classes flag.")
        advanced.add_argument("-swm", "--save_warped_modulated",
                              action='store_true',
                              help="Save warped modulated images for tissues specified in --tissue_classes flag.")
        # advanced.add_argument("-wdf", "--write_deformation_fields", nargs=2, type=bool,
        #                         help="Option to save the deformation fields from Unified Segmentation. Both inverse "
        #                              "and forward fields can be saved. Format: a list of 2 booleans. [Inverse, Forward]")

    def run_command(self, args):
        """
        """
        from clinica.utils.stream import cprint
        from .t1_volume_tissue_segmentation_pipeline import T1VolumeTissueSegmentation

        pipeline = T1VolumeTissueSegmentation(
            bids_directory=self.absolute_path(args.bids_directory),
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv),
            base_dir=self.absolute_path(args.working_directory),
        )

        pipeline.parameters.update({
            'tissue_classes': args.tissue_classes,
            'dartel_tissues': args.dartel_tissues,
            'tpm': args.tissue_probability_maps,
            'save_warped_unmodulated': args.save_warped_unmodulated,
            'save_warped_modulated': args.save_warped_modulated,
            'write_deformation_fields': [True, True],  # args.write_deformation_fields
            'save_t1_mni': True
            })

        if args.n_procs:
            pipeline.run(plugin='MultiProc', plugin_args={'n_procs': args.n_procs})
        else:
            pipeline.run()

        cprint("The " + self._name + " pipeline has completed. You can now delete the working directory (" + args.working_directory + ").")
