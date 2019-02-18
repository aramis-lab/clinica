# coding: utf8

import clinica.engine as ce


class PETVolumeCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipelines.
        """
        self._name = 'pet-volume'

    def define_description(self):
        """Define a description of this pipeline.
        """
        self._description = 'SPM-based pre-processing of PET images:\nhttp://clinica.run/doc/Pipelines/PET_Volume/'

    def define_options(self):
        """Define the sub-command arguments
        """
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES
        from colorama import Fore
        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_id)
        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("bids_directory",
                                  help='Path to the BIDS directory.')
        clinica_comp.add_argument("caps_directory",
                                  help='Path to the CAPS directory.')
        clinica_comp.add_argument("group_id",
                                  help='User-defined identifier for the provided group of subjects.')
        # Optional arguments (e.g. FWHM)
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES['OPTIONAL'])
        optional.add_argument("-psf", "--psf_tsv",
                              help='TSV file containing for each PET image its point spread function (PSF) measured at full-width and half maximum (fwhm_x, fwhm_y and fwhm_z).')
        optional.add_argument("-pet_tracer", "--pet_tracer", type=str, default='fdg', choices=['fdg', 'av45'],
                              help='PET tracer. Can be fdg or av45 (default: --pet_tracer fdg)')
        # Clinica standard arguments (e.g. --n_procs)
        clinica_opt = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_OPTIONAL'])
        clinica_opt.add_argument("-tsv", "--subjects_sessions_tsv",
                                 help='TSV file containing a list of subjects with their sessions.')
        clinica_opt.add_argument("-wd", "--working_directory",
                                 help='Temporary directory to store pipelines intermediate results')
        clinica_opt.add_argument("-np", "--n_procs",
                                 type=int,
                                 help='Number of cores used to run in parallel')
        # Advanced arguments (i.e. tricky parameters)
        advanced = self._args.add_argument_group(PIPELINE_CATEGORIES['ADVANCED'])
        advanced.add_argument("-mask", "--mask_tissues", nargs='+', type=int, default=[1, 2, 3], choices=range(1, 7),
                              help="Tissue classes (1: gray matter (GM), 2: white matter (WM), 3: cerebrospinal fluid (CSF), 4: bone, 5: soft-tissue, 6: background) to use for masking the PET image (default: GM, WM and CSF i.e. --mask_tissues 1 2 3).")
        advanced.add_argument("-threshold", "--mask_threshold", type=float, default=0.3,
                              help='Value used as threshold to binarize the tissue maps (default: --mask_threshold 0.3).')
        advanced.add_argument("-pvc_mask", "--pvc_mask_tissues",
                              nargs='+', type=int, default=[1, 2, 3], choices=range(1, 7),
                              help="Tissue classes (1: gray matter (GM), 2: white matter (WM), 3: cerebrospinal fluid (CSF), 4: bone, 5: soft-tissue, 6: background) to use as mask for PVC  (default: GM, WM and CSF i.e. --pvc_mask_tissues 1 2 3).")
        advanced.add_argument("-smooth", "--smooth", nargs='+', type=int, default=[8], metavar="FWHM_N FWHM_M",
                              help="A list of integers (e.g. --smooth 6 8) specifying the different isomorphic FWHM in millimeters to smooth the image (default: --smooth 8).")
        list_atlases = ['AAL2', 'LPBA40', 'Neuromorphometrics', 'AICHA', 'Hammers']
        advanced.add_argument("-atlases", "--atlases",
                              nargs='+', type=str, metavar="",
                              default=list_atlases, choices=list_atlases,
                              help='A list of atlases used to calculate regional mean SUVRs (default: all atlases i.e. --atlases AAL2 AICHA Hammers LPBA40 Neuromorphometrics).')

    def run_command(self, args):
        """
        """
        from tempfile import mkdtemp
        from clinica.utils.stream import cprint
        from clinica.pipelines.pet_volume.pet_volume_pipeline import PETVolume

        pipeline = PETVolume(
            bids_directory=self.absolute_path(args.bids_directory),
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv),
            group_id=args.group_id,
            fwhm_tsv=self.absolute_path(args.psf_tsv)
        )

        pipeline.parameters.update({
            'pet_type': args.pet_tracer,
            'mask_tissues': args.mask_tissues,
            'mask_threshold': args.mask_threshold,
            'pvc_mask_tissues': args.pvc_mask_tissues,
            'smooth': args.smooth,
            'atlas_list': args.atlases
        })

        if args.working_directory is None:
            args.working_directory = mkdtemp()
        pipeline.base_dir = self.absolute_path(args.working_directory)

        if args.n_procs:
            pipeline.run(plugin='MultiProc', plugin_args={'n_procs': args.n_procs})
        else:
            pipeline.run()

        cprint("The " + self._name + " pipeline has completed. You can now delete the working directory (" + args.working_directory + ").")
