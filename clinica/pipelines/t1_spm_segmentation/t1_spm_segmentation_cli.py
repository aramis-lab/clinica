# coding: utf8

"""T1 SPM Segmentation - Clinica Command Line Interface.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details: https://gitlab.icm-institute.org/aramis/clinica/wikis/docs/InteractingWithClinica.
"""


import clinica.engine as ce

__author__ = "Jorge Samper Gonzalez"
__copyright__ = "Copyright 2016-2018, The Aramis Lab Team"
__credits__ = [ "Jorge Samper Gonzalez"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


class T1SPMSegmentationCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipelines.
        """
        self._name = 't1-volume-segmentation'

    def define_description(self):
        """Define a description of this pipeline.
        """
        self._description = 'Tissue segmentation of T1w images with SPM:\nhttp://clinica.run/doc/Pipelines/T1_SPM_Full/'

    def define_options(self):
        """Define the sub-command arguments
        """
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')
        self._args.add_argument("caps_directory",
                                help='Path to the CAPS directory.')
        self._args.add_argument("-tsv", "--subjects_sessions_tsv",
                                help='TSV file containing the subjects with their sessions.')
        self._args.add_argument("-ti", "--tissue_classes", nargs='+', type=int, default=[1, 2, 3], choices=range(1, 7),
                                help="Tissue classes (gray matter, GM; white matter, WM; cerebro-spinal fluid, CSF...) to save. Up to 6 tissue classes can be saved. Ex: 1 2 3 is GM, WM and CSF")
        self._args.add_argument("-dt", "--dartel_tissues", nargs='+', type=int, default=[1, 2, 3], choices=range(1, 7),
                                help='Tissues to use for DARTEL template calculation. Ex: 1 is only GM')
        self._args.add_argument("-tpm", "--tissue_probability_maps",
                                help='Tissue probability maps to use for segmentation.')
        self._args.add_argument("-swu", "--save_warped_unmodulated", action='store_true', default=True,
                                help="Save warped unmodulated images for tissues specified in --tissue_classes")
        self._args.add_argument("-swm", "--save_warped_modulated", action='store_true',
                                help="Save warped modulated images for tissues specified in --tissue_classes")
        # self._args.add_argument("-wdf", "--write_deformation_fields", nargs=2, type=bool,
        #                         help="Option to save the deformation fields from Unified Segmentation. Both inverse and forward fields can be saved. Format: a list of 2 booleans. [Inverse, Forward]")
        self._args.add_argument("-wd", "--working_directory",
                                help='Temporary directory to store pipelines intermediate results')
        self._args.add_argument("-np", "--n_procs", type=int,
                                help='Number of cores used to run in parallel')
        self._args.add_argument("-sl", "--slurm", action='store_true',
                                help='Run the pipelines using SLURM')

    def run_command(self, args):
        """
        """
        from t1_spm_segmentation_pipeline import T1SPMSegmentation

        pipeline = T1SPMSegmentation(bids_directory=self.absolute_path(args.bids_directory),
                                     caps_directory=self.absolute_path(args.caps_directory),
                                     tsv_file=self.absolute_path(args.subjects_sessions_tsv))

        pipeline.parameters.update({'tissue_classes': args.tissue_classes,
                                    'dartel_tissues': args.dartel_tissues,
                                    'tpm': args.tissue_probability_maps,
                                    'save_warped_unmodulated': args.save_warped_unmodulated,
                                    'save_warped_modulated': args.save_warped_modulated,
                                    'write_deformation_fields': [True, True],  # args.write_deformation_fields
                                    'save_t1_mni': True
                                    })

        pipeline.base_dir = self.absolute_path(args.working_directory)

        if args.n_procs:
            pipeline.run(plugin='MultiProc', plugin_args={'n_procs': args.n_procs})
        elif args.slurm:
            pipeline.run(plugin='SLURM')
        else:
            pipeline.run()
