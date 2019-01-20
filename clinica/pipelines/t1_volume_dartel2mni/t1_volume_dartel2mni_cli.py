# coding: utf8

import clinica.engine as ce

__author__ = "Jorge Samper Gonzalez"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Jorge Samper Gonzalez"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


class T1VolumeDartel2MNICLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipelines.
        """
        self._name = 't1-volume-dartel2mni'

    def define_description(self):
        """Define a description of this pipeline.
        """
        self._description = 'Register DARTEL template to MNI space with SPM:\n'\
                            + 'http://clinica.run/doc/Pipelines/T1_Volume/'

    def define_options(self):
        """Define the sub-command arguments
        """
        self._args.add_argument("bids_directory",
                                help='Path to the BIDS directory.')
        self._args.add_argument("caps_directory",
                                help='Path to the CAPS directory.')
        self._args.add_argument("group_id",
                                help='User-defined identifier for the provided group of subjects.')
        self._args.add_argument("-tsv", "--subjects_sessions_tsv",
                                help='TSV file containing a list of subjects with their sessions.')
        self._args.add_argument("-t", "--tissues", nargs='+', type=int, default=[1, 2, 3], choices=range(1, 7),
                                help='Tissues to register into MNI space. Ex: 1 is only GM')
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
                                help='Temporary directory to store pipelines intermediate results')
        self._args.add_argument("-np", "--n_procs", type=int,
                                help='Number of cores used to run in parallel')

    def run_command(self, args):
        """
        """

        from clinica.pipelines.t1_volume_dartel2mni.t1_volume_dartel2mni_pipeline import T1VolumeDartel2MNI

        pipeline = T1VolumeDartel2MNI(
                bids_directory=self.absolute_path(args.bids_directory),
                caps_directory=self.absolute_path(args.caps_directory),
                tsv_file=self.absolute_path(args.subjects_sessions_tsv),
                group_id=args.group_id
                )

        pipeline.parameters.update({
            'tissues': args.tissues,
            # 'bounding_box': None,
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
