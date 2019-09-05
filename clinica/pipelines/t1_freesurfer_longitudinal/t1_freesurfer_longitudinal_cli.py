# coding: utf8

import clinica.engine as ce

__author__ = "Alexis Guyot"
__copyright__ = "Copyright 2016-2019, The Aramis Lab Team"
__credits__ = ["Alexis Guyot"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Alexis Guyot"
__email__ = "alexis.guyot@icm-institute.org"
__status__ = "Development"


class T1FreeSurferLongitudinalCLI(ce.CmdParser):

    def define_name(self):
        """Define the sub-command name to run this pipelines.
        """
        self._name = 't1-freesurfer-longitudinal'

    def define_description(self):
        """Define a description of this pipeline.
        """
        self._description = 'Longitudinal pre-processing of T1w images with FreeSurfer:\nhttp://clinica.run/doc/Pipelines/T1_FreeSurfer_Longitudinal/'

    def define_options(self):
        """Define the sub-command arguments
        """
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES
        # clinica compulsory arguments (e.g. bids, caps, group_id)
        clinica_comp = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_COMPULSORY'])
        clinica_comp.add_argument("caps_directory",
                                  help='path to the caps directory.')
        # clinica standard arguments (e.g. --n_procs)
        clinica_opt = self._args.add_argument_group(PIPELINE_CATEGORIES['CLINICA_OPTIONAL'])
        clinica_opt.add_argument("-tsv", "--subjects_sessions_tsv",
                                 help='tsv file containing a list of subjects with their sessions.')
        clinica_opt.add_argument("-wd", "--working_directory",
                                 help='temporary directory to store pipelines intermediate results')
        clinica_opt.add_argument("-np", "--n_procs",
                                 metavar=('n'), type=int,
                                 help='number of cores used to run in parallel')
        # clinica pipeline-specific arguments
        clinica_opt.add_argument("-oc", "--overwrite_caps",
                                 type=str, default='false',
                                 help='Overwrite existing data in CAPS dir (Possible values: true, True, TRUE, false, False, FALSE. Default: -oc false)')

    def run_command(self, args):
        """Run the pipelines with defined args
        """
        import nipype.pipeline.engine as npe
        from nipype.pipeline.engine import Workflow
        from .t1_freesurfer_template_pipeline import T1FreeSurferTemplate
        from .t1_freesurfer_longitudinal_correction_pipeline import T1FreeSurferLongitudinalCorrection
        from tempfile import mkdtemp

        template_pipeline = T1FreeSurferTemplate(
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv))
        longcorr_pipeline = T1FreeSurferLongitudinalCorrection(
            caps_directory=self.absolute_path(args.caps_directory),
            tsv_file=self.absolute_path(args.subjects_sessions_tsv))

        # add working directory and n_procs to pipeline parameters so we
        # know, while running the pipeline, which arguments the user
        # passed to the CLI
        template_pipeline.parameters = {
            'working_directory': args.working_directory,
            'n_procs': args.n_procs,
            'overwrite_caps': args.overwrite_caps
            }
        longcorr_pipeline.parameters = {
            'working_directory': args.working_directory,
            'n_procs': args.n_procs,
            'overwrite_caps': args.overwrite_caps
            }

        # make sure if working_directory is not defined, using a temp
        # folder to the working directory.
        if args.working_directory is None:
            args.working_directory = mkdtemp()
        template_pipeline.base_dir = self.absolute_path(args.working_directory)
        longcorr_pipeline.base_dir = self.absolute_path(args.working_directory)

        # separately build the two template and longitudinal-correction pipelines
        template_pipeline.build()
        longcorr_pipeline.build()

        # create overall workflow
        longitudinal_workflow = npe.Workflow(name='T1FreeSurferLongitudinal')
        longitudinal_workflow.base_dir = self.absolute_path(args.working_directory)
        # connect the template pipeline to the longitudinal-correction pipeline
        longitudinal_workflow.connect(
            template_pipeline, '5_sendto_longcorr.out_unpcssd_sublist',
            longcorr_pipeline, '0_receivefrom_template.unpcssd_sublist')
        longitudinal_workflow.connect(
            template_pipeline, '5_sendto_longcorr.out_pcssd_capstargetlist',
            longcorr_pipeline, '0_receivefrom_template.pcssd_capstargetlist')
        longitudinal_workflow.connect(
            template_pipeline, '5_sendto_longcorr.out_overwrite_tsv',
            longcorr_pipeline, '0_receivefrom_template.overwrite_tsv')
        # run overall workflow
        if args.n_procs:
            Workflow.run(
                longitudinal_workflow,
                plugin_args={'n_procs': args.n_procs},
                plugin='MultiProc')
        else:
            Workflow.run(longitudinal_workflow)
