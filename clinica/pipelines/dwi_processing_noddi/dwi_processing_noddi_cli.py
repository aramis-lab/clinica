# coding: utf8

__author__ = "Junhao Wen"
__copyright__ = "Copyright 2016-2018, The Aramis Lab Team"
__credits__ = ["Junhao Wen"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Junhao Wen"
__email__ = "Junhao.Wen@inria.fr"
__status__ = "Development"

import clinica.engine as ce


class DwiProcessingNoddiCLI(ce.CmdParser):

    def __init__(self):
        super(DwiProcessingNoddiCLI, self).__init__()

    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """

        self._name = 'dwi-processing-noddi'

    def define_description(self):
        """Define a description of this pipeline.
        """
        self._description = ('NODDI-based processing of DWI datasets:\n'
                            'http://clinica.run/doc/DWIProcessing')

    def define_options(self):
        """Define the sub-command arguments
        """

        from clinica.engine.cmdparser import PIPELINE_CATEGORIES

        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])

        clinica_comp.add_argument("caps_directory",
                                help='Path to the CAPS directory.')
        clinica_comp.add_argument("list_bvalues", type=str,
                                help='String listing all the shells (i.e. the b-values) in the corrected DWI datasets comma separated (e.g, 0,300,700,2200)')
        # Optional arguments
        clinica_opt = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_OPTIONAL'])

        clinica_opt.add_argument("-wd", "--working_directory",
                                help='Temporary directory to store pipeline intermediate results')
        clinica_opt.add_argument("-np", "--n_procs", type=int, default=4,
                                help='Number of cores used to run in parallel')
        clinica_opt.add_argument("-tsv", "--subjects_sessions_tsv",
                                 help='TSV file containing a list of subjects with their sessions.')

    def run_command(self, args):
        """
        """

        from tempfile import mkdtemp
        from clinica.pipelines.dwi_processing_noddi.dwi_processing_noddi_pipeline import DwiProcessingNoddi

        pipeline = DwiProcessingNoddi(
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv))

        from clinica.utils.check_dependency import check_noddi_matlab_toolbox, check_nifti_matlib_toolbox
        noddi_matlab_toolbox = check_noddi_matlab_toolbox()
        nifti_matlib_toolbox = check_nifti_matlib_toolbox()

        pipeline.parameters = {
            'bvalue_str': dict(
                [('bvalue_str', args.list_bvalues)]),
            'n_procs': dict(
                [('n_procs', args.n_procs or 4)]),
            'noddi_toolbox_dir': dict(
                [('noddi_toolbox_dir', noddi_matlab_toolbox)]),
            'nifti_matlib_dir': dict(
                [('nifti_matlib_dir', nifti_matlib_toolbox)]),
        }

        if args.working_directory is None:
            args.working_directory = mkdtemp()
        pipeline.base_dir = self.absolute_path(args.working_directory)
        if args.n_procs:
            pipeline.run(plugin='MultiProc',
                         plugin_args={'n_procs': args.n_procs})
        else:
            pipeline.run()
