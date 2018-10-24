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


class DwiPreprocessingNoddiCLI(ce.CmdParser):

    def __init__(self):
        super(DwiPreprocessingNoddiCLI, self).__init__()

    def define_name(self):
        """Define the sub-command name to run this pipeline.
        """

        self._name = 'dwi-preprocessing-multi-shell'

    def define_description(self):
        """Define a description of this pipeline.
        """
        self._description = ('Preprocessing of raw DWI datasets with '
                              'multi-shell acquisitions and opposite phase '
                              'encoding directions:\n'
                              'http://clinica.run/doc/Pipelines/DWI_Preprocessing/')

    def define_options(self):
        """Define the sub-command arguments
        """

        from clinica.engine.cmdparser import PIPELINE_CATEGORIES

        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("bids_directory",
                                help='Path to the BIDS directory.')
        clinica_comp.add_argument("caps_directory",
                                help='Path to the CAPS directory.')
        clinica_comp.add_argument("echo_spacing", type=float,
                                help='The echo spacing such that EffectiveEchoSpacing=echo_spacing/acceleration_factor (see BIDS specifications and https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/topup/Faq for details).')
        clinica_comp.add_argument("acceleration_factor", type=int,
                                help='Acceleration factor such that EffectiveEchoSpacing=echo_spacing/acceleration_factor (see BIDS specifications and https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/topup/Faq for details).')
        clinica_comp.add_argument("phase_encoding_direction", type=str,
                                help='Phase encoding direction using FSL convention (e.g. y- or y). Be careful, this is currently not the BIDS convention (i.e. j- or j).')
        clinica_comp.add_argument("phase_encoding_direction_alternative", type=str,
                                help='The opposite phase encoding direction (e.g. if phase_encoding_direction is y- then phase_encoding_direction_alternative will be y).')
        clinica_comp.add_argument("epi_factor", type=int,
                                help='EPI factor (e.g. 128) used for the computation of the TotalReadoutTime (see https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy/Faq for details)')

        # Optional arguments
        clinica_opt = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_OPTIONAL'])
        clinica_opt.add_argument("-wd", "--working_directory",
                                help='Temporary directory to store pipeline intermediate results')
        clinica_opt.add_argument("-np", "--n_procs", type=int,
                                help='Number of cores used to run in parallel')
        clinica_opt.add_argument("-tsv", "--subjects_sessions_tsv",
                                help='TSV file containing a list of subjects with their sessions.')

    def run_command(self, args):
        """
        """

        from tempfile import mkdtemp
        from clinica.pipelines.dwi_preprocessing_noddi.dwi_preprocessing_noddi_pipeline import DwiPreprocessingNoddi
        import os, errno
        from clinica.iotools.utils.data_handling import create_subs_sess_list
        from nipype import config
        cfg = dict(execution={'parameterize_dirs': False})
        config.update_config(cfg)

        if args.subjects_sessions_tsv is None:
            try:
                temp_dir = mkdtemp()
            except OSError as exception:
                if exception.errno != errno.EEXIST:
                    raise
            create_subs_sess_list(args.bids_directory, temp_dir)
            args.subjects_sessions_tsv = os.path.join(temp_dir, 'subjects_sessions_list.tsv')

        pipeline = DwiPreprocessingNoddi(
            bids_directory=self.absolute_path(args.bids_directory),
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv))

        pipeline.parameters = {
            # pass these args by using self.parameters in a dictionary
            'epi_param': dict([('echospacing', args.echo_spacing), ('acc_factor', args.acceleration_factor), ('enc_dir', args.phase_encoding_direction), ('epi_factor', args.epi_factor)]),
            'alt_epi_params': dict([('echospacing', args.echo_spacing), ('acc_factor', args.acceleration_factor), ('enc_dir_alt', args.phase_encoding_direction_alternative), ('epi_factor', args.epi_factor)])
        }

        if args.working_directory is None:
            args.working_directory = mkdtemp()
        pipeline.base_dir = self.absolute_path(args.working_directory)

        if args.n_procs:
            pipeline.run(plugin='MultiProc',
                         plugin_args={'n_procs': args.n_procs})
        else:
            pipeline.run()
