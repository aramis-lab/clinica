# coding: utf8

__author__ = "Simona Bottani"
__copyright__ = "Copyright 2016-2018, The Aramis Lab Team"
__credits__ = ["Simona Bottani"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Simona Bottani"
__email__ = "simona.bottani@icm-institute.org"
__status__ = "Development"
import clinica.engine as ce


class T1VolumeParcellationCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """
        self._name = 't1-volume-parcellation'

    def define_description(self):
        """Define a description of this pipeline.
        """
        self._description = 'Computation of mean GM concentration for a set of regions:\nhttp://clinica.run/doc/Pipelines/T1_Volume/'

    def define_options(self):
        """Define the sub-command arguments
        """

        #self._args.add_argument("bids_directory",
        #                        help='Path to the BIDS directory.')
        self._args.add_argument("caps_directory",
                                help='Path to the CAPS directory.')
        self._args.add_argument("group_id",
                                help='User-defined identifier for the provided group of subjects.')
        self._args.add_argument("-tsv", "--subjects_sessions_tsv",
                                help='TSV file containing a list of subjects with their sessions.')
        #self._args.add_argument("-im_type", "--image_type", type = str, default = 'T1',
        #                        choices =['T1', 'pet'],
        #                        help = 'image type. Possible values are T1 and pet')
        self._args.add_argument("-pet", "--pet_type", type = str, default = 'fdg',
                                choices = ['fdg', 'av45'],
                                help = 'PET image type. Possible values are fdg and av45.')
        self._args.add_argument("-m", "--modulate", default='on',
                                choices=['on', 'off'],
                                help='A boolean. Modulate output images - no modulation preserves concentrations')
        self._args.add_argument("-atlases", "--atlases", nargs='+', type=str,
                                default=['AAL2', 'LPBA40', 'Neuromorphometrics', 'AICHA', 'Hammers'],
                                choices=['AAL2', 'LPBA40', 'Neuromorphometrics', 'AICHA', 'Hammers'],
                                help='A list of atlases to use to calculate the mean GM concentration at each region')
        self._args.add_argument("-wd", "--working_directory",
                                help='Temporary directory to store pipeline intermediate results')
        self._args.add_argument("-np", "--n_procs", type=int,
                                help='Number of cores used to run in parallel')

    def run_command(self, args):
        """
        """

        from tempfile import mkdtemp
        from t1_volume_parcellation_pipeline import T1VolumeParcellation

        # Most of the time, you will want to instantiate your pipeline with a
        # BIDS and CAPS directory as inputs:
        pipeline = T1VolumeParcellation(
             bids_directory='./4',
             caps_directory=self.absolute_path(args.caps_directory),
             tsv_file=self.absolute_path(args.subjects_sessions_tsv))
        #pipeline = spm_parcellation()
        pipeline.parameters = {
            # Add your own pipeline parameters here to use them inside your
            # pipeline. See the file `spm_parcellation_pipeline.py` to
            # see an example of use.
            'group_id': args.group_id,
            'pet_type': args.pet_type,
            'modulate': args.modulate,
            'atlases': args.atlases,
            'wd': self.absolute_path(args.working_directory),
            'n_procs': args.n_procs
        }

        if args.working_directory is None:
            args.working_directory = mkdtemp()
        pipeline.base_dir = self.absolute_path(args.working_directory)
        if args.n_procs:
            pipeline.run(plugin='MultiProc', plugin_args={'n_procs': args.n_procs})
        else:
            pipeline.run()
