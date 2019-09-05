# coding: utf8

__author__ = "Alexis Guyot"
__copyright__ = "Copyright 2016-2019, The Aramis Lab Team"
__credits__ = ["Alexis Guyot"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Alexis Guyot"
__email__ = "alexis.guyot@icm-institute.org"
__status__ = "Development"

import clinica.pipelines.engine as cpe


class T1FreeSurferLongitudinalCorrection(cpe.Pipeline):
    """FreeSurfer Longitudinal correction class

    Creates a pipeline that runs the Freesurfer longitudinal
    (correction) processing module for each subjects in a .tsv-defined
    list of subjects/sessions. This requires a prior run of
    t1-freesurfer-cross-sectional on the .tsv, followed by a run of
    t1-freesurfer-template on the same .tsv. For each subject, all the
    timepoints (sessions) are re-processed based on a template computed
    with t1-freesurfer-template for that specific subject.

    Warnings:
        - A WARNING.

    Todos: N/A

    Args:
        caps_directory: Path to the CAPS directory
        tsv: TSV file containing the subjects with their sessions
        wd: Temporary directory to store pipelines intermediate results
        np: Number of cores used to run in parallel

    Returns:
        A clinica pipeline object containing the T1 FreeSurfer
            pipeline

    Raises:


    Example:
        >>> from t1_freesurfer_longitudinal_correction import T1FreeSurferLongitudinalCorrection
        >>> pipelines = T1FreeSurferLongitudinalCorrection('~/MYDATASET_BIDS', '~/MYDATASET_CAPS', 'TSV')
        >>> pipelines.base_dir = '/tmp/'
        >>> pipelines.run()
    """

    def check_custom_dependencies(self):
        """Check dependencies that cannot be listed in `info.json`.
        """
        pass

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Returns:
            A list of (string) input fields name.
        """

        return ['subject_list',
                'session_list',
                'caps_target_list',
                'caps_dir',
                'overwrite_warning',
                'overwrite_caps']

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """

        return ['subject',
                'session',
                'subject_dir',
                'caps_target',
                'caps_dir',
                'overwrite_warning',
                'overwrite_caps',
                'longitudinal_result',
                'stats_path']

    def build_input_node(self):
        """Build and connect an input node to the pipeline.

        This node is supposedly used to load BIDS inputs when this
        pipeline is not already connected to the output of a previous
        Clinica pipeline. For the purpose of the example, we simply
        read input arguments given by the command line interface and
        transmitted here through the `self.parameters` dictionary and
        pass it to the `self.input_node` to further by used as input of
        the core nodes.
        """
        import clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_longitudinal_correction_utils as utils
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        # Step 0: receive data from the template pipeline
        #         (subjects that have been detected as
        #         processed/non-processed and corresponding sessions and
        #         CAPS locations)
        # ======
        receivefrom_template_node_name = '0_receivefrom_template'
        receivefrom_template_node = npe.Node(
            name=receivefrom_template_node_name,
            interface=nutil.IdentityInterface(
                fields=['unpcssd_sublist',
                        'pcssd_capstargetlist',
                        'overwrite_tsv'])
            )

        # check if cross-sectional pipeline run on all subjects
        checkinput_node_name = '1_check_input'
        checkinput_node = npe.Node(name=checkinput_node_name,
                                   interface=nutil.Function(
                                       input_names=['in_caps_dir',
                                                    'in_subject_list',
                                                    'in_session_list',
                                                    'in_unpcssd_sublist',
                                                    'in_pcssd_capstargetlist',
                                                    'in_overwrite_caps',
                                                    'in_working_directory',
                                                    'in_n_procs',
                                                    'in_overwrite_tsv'],
                                       output_names=['out_subject_list',
                                                     'out_session_list',
                                                     'out_caps_target_list',
                                                     'out_overwrite_warning'],
                                       function=utils.process_input_node))
        checkinput_node.inputs.in_caps_dir = self.caps_directory
        checkinput_node.inputs.in_subject_list = self.subjects
        checkinput_node.inputs.in_session_list = self.sessions
        checkinput_node.inputs.in_working_directory = self.parameters[
            'working_directory']
        checkinput_node.inputs.in_overwrite_caps = self.parameters[
            'overwrite_caps']
        checkinput_node.inputs.in_n_procs = self.parameters['n_procs']

        self.connect([
            (
                receivefrom_template_node, checkinput_node,
                [('unpcssd_sublist', 'in_unpcssd_sublist'),
                 ('pcssd_capstargetlist', 'in_pcssd_capstargetlist'),
                 ('overwrite_tsv', 'in_overwrite_tsv')]),
            (
                checkinput_node, self.input_node,
                [('out_subject_list', 'subject_list'),
                 ('out_session_list', 'session_list'),
                 ('out_caps_target_list', 'caps_target_list'),
                 ('out_overwrite_warning', 'overwrite_warning')])
        ])

    def build_output_node(self):
        """copy the longitudinal-correction folder to CAPS DIR

        Will take all the data stored into the node in working directory
        that contains the built templates and copy those to the CAPS
        directory
        """
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_longitudinal_correction_utils as utils

        # define the node to copy data to CAPS folder
        copy2caps_node_name = '7_copy2caps'
        copy2caps_node = npe.MapNode(
            name=copy2caps_node_name,
            interface=nutil.Function(
                input_names=['in_subject',
                             'in_session',
                             'in_subject_dir',
                             'in_caps_target',
                             'in_caps_dir',
                             'in_overwrite_warning',
                             'in_overwrite_caps',
                             'in_stats_path'],
                output_names=['out_copy2_caps'],
                function=utils.copy_to_caps
                ),
            iterfield=['in_subject',
                       'in_session',
                       'in_caps_target',
                       'in_stats_path'])
        copy2caps_node.inputs.in_caps_dir = self.caps_directory
        copy2caps_node.inputs.in_overwrite_caps = self.parameters[
            'overwrite_caps']

        self.connect([
            (
                self.output_node, copy2caps_node,
                [('subject', 'in_subject'),
                 ('session', 'in_session'),
                 ('subject_dir', 'in_subject_dir'),
                 ('caps_target', 'in_caps_target'),
                 ('overwrite_warning', 'in_overwrite_warning'),
                 ('stats_path', 'in_stats_path')])
        ])

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipelines.
        """
        import clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_longitudinal_correction_utils as utils
        import clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_template_utils as template_utils
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from nipype.interfaces.freesurfer.preprocess import ReconAll

        # Per-subject creation of the template
        # Step 0: receivefrom_template_node
        # Step 1: check_input_node
        # Step 2: get longitudinal subfolder names
        # ======
        long_subdirname_node_name = '2_long_subdirname'
        long_subdirname_node = npe.Node(
            name=long_subdirname_node_name,
            interface=nutil.Function(
                input_names=['in_subject_list',
                             'in_session_list'],
                output_names=['out_longsubdirnames_dict'],
                function=template_utils.get_longsubdir_dict)
            )

        # Step 3: create symbolic links for each subject
        # ======
        symlink_node_name = '3_symlink'
        symlink_node = npe.Node(name=symlink_node_name,
                                interface=nutil.Function(
                                    input_names=['in_caps_dir',
                                                 'in_subject_list',
                                                 'in_session_list',
                                                 'in_longsubdirnames_dict'],
                                    output_names=['out_symlink_path'],
                                    function=utils.create_symlinks))
        symlink_node.inputs.in_caps_dir = self.caps_directory

        # Step 4: prepare longitudinal recon-all flags (prior to
        #         running recon-all -long)
        # ======
        prepare_longitudinal_node_name = '4_prepare_longitudinal'
        prepare_longitudinal_node = npe.MapNode(
            name=prepare_longitudinal_node_name,
            interface=nutil.Function(
                input_names=['in_subject',
                             'in_session'],
                output_names=['out_reconalllong_flags'],
                function=utils.get_reconalllong_flags
                ),
            iterfield=['in_subject',
                       'in_session'])

        # Step 5: run freesurfer longitudinal for all {subject,session}
        #         (i.e. run recon-all -long)
        # ======
        run_longitudinal_node_name = '5_run_longitudinal'
        run_longitudinal_node = npe.MapNode(
            name=run_longitudinal_node_name,
            interface=nutil.Function(
                input_names=['in_subject',
                             'in_session',
                             'in_reconalllong_flags',
                             'in_symlink_path'],
                output_names=['out_longitudinal_result'],
                function=utils.run_reconalllong
                ),
            iterfield=['in_subject',
                       'in_session',
                       'in_reconalllong_flags'])

        # Step 6: Create statistics .tsv files for all {subject,session}
        # ======
        write_stats_node_name = '6_write_stats'
        write_stats_node = npe.MapNode(
            name=write_stats_node_name,
            interface=nutil.Function(
                input_names=['in_subject',
                             'in_session',
                             'in_symlink_path',
                             'in_longitudinal_result'],
                output_names=['out_stats_path'],
                function=utils.write_stats_files,
                ),
            iterfield=['in_subject',
                       'in_session'])

        # Connection
        # ==========
        self.connect([
            # long_subdirname_node
            (
                self.input_node, long_subdirname_node,
                [('subject_list', 'in_subject_list'),
                 ('session_list', 'in_session_list')]),
            # symlink_node
            (
                self.input_node, symlink_node,
                [('subject_list', 'in_subject_list'),
                 ('session_list', 'in_session_list')]),
            (
                long_subdirname_node, symlink_node,
                [('out_longsubdirnames_dict', 'in_longsubdirnames_dict')]),
            # prepare_longitudinal_node
            (
                self.input_node, prepare_longitudinal_node,
                [('subject_list', 'in_subject'),
                 ('session_list', 'in_session')]),
            # run_longitudinal_node
            (
                self.input_node, run_longitudinal_node,
                [('subject_list', 'in_subject'),
                 ('session_list', 'in_session')]),
            (
                prepare_longitudinal_node, run_longitudinal_node,
                [('out_reconalllong_flags', 'in_reconalllong_flags')]),
            (
                symlink_node, run_longitudinal_node,
                [('out_symlink_path', 'in_symlink_path')]),
            # write_stats_node
            (
                self.input_node, write_stats_node,
                [('subject_list', 'in_subject'),
                 ('session_list', 'in_session')]),
            (
                symlink_node, write_stats_node,
                [('out_symlink_path', 'in_symlink_path')]),
            (
                run_longitudinal_node, write_stats_node,
                [('out_longitudinal_result', 'in_longitudinal_result')]),
            # self.output_node (copy2caps_node)
            (
                self.input_node, self.output_node,
                [('subject_list', 'subject'),
                 ('session_list', 'session')]),
            (
                symlink_node, self.output_node,
                [('out_symlink_path', 'subject_dir')]),
            (
                self.input_node, self.output_node,
                [('caps_target_list', 'caps_target'),
                 ('overwrite_warning', 'overwrite_warning')]),
            (
                write_stats_node, self.output_node,
                [('out_stats_path', 'stats_path')])
        ])

        return self
