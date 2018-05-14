# coding: utf8

"""T1 SPM Full Prep - Clinica Command Line Interface.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details: https://gitlab.icm-institute.org/aramislab/clinica/wikis/docs/InteractingWithClinica.
"""


import clinica.engine as ce

__author__ = "Jorge Samper Gonzalez"
__copyright__ = "Copyright 2016, The Aramis Lab Team"
__credits__ = ["Jorge Samper Gonzalez"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


class T1VolumeNewTemplateCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipelines.
        """
        self._name = 't1-volume-new-template'

    def define_description(self):
        """Define a description of this pipeline.
        """
        self._description = 'SPM-based pre-processing of T1w images:\nhttp://clinica.run/doc/Pipelines/T1_SPM_Full/'

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
        clinica_comp.add_argument("group_id",
                                  help='Current group name')
        # Optional arguments (e.g. FWHM)
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES['OPTIONAL'])
        optional.add_argument("-fwhm", "--fwhm",
                              nargs='+', type=int, default=[8],
                              help="A list of integers specifying the different isomorphic fwhm in millimeters to smooth the image")
        # Clinica standard arguments (e.g. --n_procs)
        clinica_opt = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_OPTIONAL'])
        clinica_opt.add_argument("-tsv", "--subjects_sessions_tsv",
                                 metavar=('participants.tsv'),
                                 help='TSV file containing the subjects with their sessions.')
        clinica_opt.add_argument("-wd", "--working_directory",
                                 metavar=('Working_Directory'),
                                 help='Temporary directory to store pipelines intermediate results')
        clinica_opt.add_argument("-np", "--n_procs",
                                 metavar=('N'),
                                 type=int,
                                 help='Number of cores used to run in parallel')
        clinica_opt.add_argument("-sl", "--slurm",
                                 action='store_true',
                                 help='Run the pipelines using SLURM')
        # Advanced arguments (i.e. tricky parameters)
        advanced = self._args.add_argument_group(PIPELINE_CATEGORIES['ADVANCED'])
        advanced.add_argument("-ti", "--tissue_classes",
                              metavar=('1 2 3 4 5 6'),
                              nargs='+', type=int, default=[1, 2, 3], choices=range(1, 7),
                              help="Tissue classes (gray matter, GM; white matter, WM; cerebro-spinal fluid, CSF...) to save. Up to 6 tissue classes can be saved. Ex: 1 2 3 is GM, WM and CSF")
        advanced.add_argument("-dt", "--dartel_tissues",
                              metavar=('1 2 3 4 5 6'),
                              nargs='+', type=int, default=[1, 2, 3], choices=range(1, 7),
                              help='Tissues to use for DARTEL template calculation. Ex: 1 is only GM')
        advanced.add_argument("-tpm", "--tissue_probability_maps",
                              metavar=('TissueProbabilityMap.nii'),
                              help='Tissue probability maps to use for segmentation.')
        advanced.add_argument("-swu", "--save_warped_unmodulated",
                              action='store_true', default=True,
                              help="Save warped unmodulated images for tissues specified in --tissue_classes")
        advanced.add_argument("-swm", "--save_warped_modulated",
                              action='store_true',
                              help="Save warped modulated images for tissues specified in --tissue_classes")
        advanced.add_argument("-m", "--modulate",
                              type=bool, default=True,
                              metavar=('True/False'),
                              help='A boolean. Modulate output images - no modulation preserves concentrations')
        advanced.add_argument("-vs", "--voxel_size",
                              metavar=('float'),
                              nargs=3, type=float,
                              help="A list of 3 floats specifying voxel sizes for each dimension of output image")
        list_atlases = ['AAL2', 'LPBA40', 'Neuromorphometrics', 'AICHA', 'Hammers']
        advanced.add_argument("-atlases", "--atlases",
                              nargs='+', type=str,
                              default=list_atlases, choices=list_atlases,
                              help='A list of atlases to use to calculate the mean GM concentration at each region')

    def run_command(self, args):
        """
        """
        from t1_volume_new_template_pipeline import T1VolumeNewTemplate

        pipeline = T1VolumeNewTemplate (bids_directory=self.absolute_path(args.bids_directory),
                                 caps_directory=self.absolute_path(args.caps_directory),
                                 tsv_file=self.absolute_path(args.subjects_sessions_tsv),
                                 group_id=args.group_id
                                 )

        pipeline.parameters.update({
            'tissue_classes': args.tissue_classes,
            'dartel_tissues': args.dartel_tissues,
            'tpm': args.tissue_probability_maps,
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
        elif args.slurm:
            pipeline.run(plugin='SLURM')
        else:
            pipeline.run()
