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


class T1FreeSurferTemplate(cpe.Pipeline):
    """FreeSurfer Longitudinal template class

    Creates a pipelines that runs the Freesurfer longitudinal template
    initialisation module for each subjects in a .tsv-defined list of
    subjects/sessions. This requires a prior-run of
    t1-freesurfer-cross-sectional on the .tsv.


    Warnings:
        - A WARNING.

    Todos: N/A

    Args:
        caps_directory: Path to the CAPS directory.
        tsv: TSV file containing the subjects with their sessions
        wd: Temporary directory to store pipelines intermediate results
        np: Number of cores used to run in parallel

    Returns:
        A clinica pipeline object containing the T1 FreeSurfer pipeline.

    Raises:


    Example:
        >>> from t1_freesurfer_template import T1FreeSurferTemplate
        >>> pipelines = T1FreeSurferTemplate('~/MYDATASET_BIDS', '~/MYDATASET_CAPS', 'TSV')
        >>> pipelines.base_dir = '/tmp/'
        >>> pipelines.run()
    """

    def check_custom_dependencies(self):
        """Check dependencies that cannot be listed in `info.json`
        """
        pass

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Returns:
            A list of (string) input fields name.
        """

        return ['subject_list',
                'session_list2',
                'caps_target_list',
                'unpcssd_sublist',
                'pcssd_capstargetlist',
                'overwrite_tsv']

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """

        return []

    def build_input_node(self):
        """Build and connect an input node to the pipeline.

        This node is supposedly used to load BIDS inputs when this
        pipeline is not already connected to the output of a previous
        Clinica pipeline.  For the purpose of the example, we simply
        read input arguments given by the command line interface and
        transmitted here through the `self.parameters` dictionary and
        pass it to the `self.input_node` to further by used as input of
        the core nodes.
        """
        import os
        import clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_template_utils as utils
#        import t1_freesurfer_template_utils as utils
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        # check if cross-sectional pipeline run on all subjects
        checkinput_node_name = '0_check_input'
        checkinput_node = npe.Node(name=checkinput_node_name,
                                   interface=nutil.Function(
                                       input_names=['in_caps_dir',
                                                    'in_subject_list',
                                                    'in_session_list',
                                                    'in_overwrite_caps',
                                                    'in_working_directory',
                                                    'in_n_procs',
                                                    'in_rundir'],
                                       output_names=['out_subject_list',
                                                     'out_session_list2',
                                                     'out_caps_target_list',
                                                     'out_unpcssd_sublist',
                                                     'out_pcssd_capstargetlist',
                                                     'out_overwrite_tsv'],
                                       function=utils.process_input_node))
        checkinput_node.inputs.in_caps_dir = self.caps_directory
        checkinput_node.inputs.in_subject_list = self.subjects
        checkinput_node.inputs.in_session_list = self.sessions
        checkinput_node.inputs.in_working_directory = self.parameters[
            'working_directory']
        checkinput_node.inputs.in_overwrite_caps = self.parameters[
            'overwrite_caps']
        checkinput_node.inputs.in_n_procs = self.parameters['n_procs']
        checkinput_node.inputs.in_rundir = os.getcwd()

        self.connect([
            (
                checkinput_node, self.input_node,
                [('out_subject_list', 'subject_list'),
                 ('out_session_list2', 'session_list2'),
                 ('out_caps_target_list', 'caps_target_list'),
                 ('out_unpcssd_sublist', 'unpcssd_sublist'),
                 ('out_pcssd_capstargetlist', 'pcssd_capstargetlist'),
                 ('out_overwrite_tsv', 'overwrite_tsv')])
        ])

    def build_output_node(self):
        """copy the template folder to CAPS DIR

        Will take all the data stored into the node in working directory
        that contains the built templates and copy those to the CAPS
        directory
        """
        pass

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline.
        """

        import clinica.pipelines.t1_freesurfer_longitudinal.t1_freesurfer_template_utils as utils
#        import t1_freesurfer_template_utils as utils
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        # Per-subject creation of the template

        # Step 1: create folder to store recon-all -base results
        # ======
        store_reconallbase_node_name = '1_store_reconallbase'
        store_reconallbase_node = npe.Node(
            name=store_reconallbase_node_name,
            interface=nutil.Function(
                input_names=['in_subject_list'],
                output_names=['out_workdirstore_path'],
                function=utils.store_reconallbase_results))

        # Step 2: prepare unbiased template flags for all subjects
        #         (i.e. create flags for recon-all -base)
        # ======
        prepare_template_node_name = '2_prepare_template'
        prepare_template_node = npe.MapNode(
            name=prepare_template_node_name,
            interface=nutil.Function(
                input_names=['in_subject',
                             'in_session_list'],
                output_names=['out_reconallbase_flags'],
                function=utils.get_reconallbase_flags
                ),
            iterfield=['in_subject',
                       'in_session_list'])

        # Step 3: create unbiased template for all subjects (i.e. run
        #         recon-all -base)
        # ======
        create_template_node_name = '3_create_template'
        create_template_node = npe.MapNode(
            name=create_template_node_name,
            interface=nutil.Function(
                input_names=['in_caps_dir',
                             'in_subject',
                             'in_session_list',
                             'in_reconallbase_flags',
                             'in_workdirstore_path'],
                output_names=['out_template_created'],
                function=utils.run_reconallbase
                ),
            iterfield=['in_subject',
                       'in_session_list',
                       'in_reconallbase_flags'])
        create_template_node.inputs.in_caps_dir = self.caps_directory

        # define the node to copy data to CAPS folder
        copy2caps_node_name = '4_copy2caps'
        copy2caps_node = npe.MapNode(
            name=copy2caps_node_name,
            interface=nutil.Function(
                input_names=['in_subject',
                             'in_subject_dir',
                             'in_caps_target',
                             'in_caps_dir',
                             'in_overwrite_caps',
                             'in_template_created'],
                output_names=['out_copy2_caps'],
                function=utils.copy_to_caps
                ),
            iterfield=['in_subject',
                       'in_caps_target'])
        copy2caps_node.inputs.in_caps_dir = self.caps_directory
        copy2caps_node.inputs.in_overwrite_caps = self.parameters[
            'overwrite_caps']

        # send data to the longitudinal correction pipeline (subjects
        # that have been detected as processed/non-processed and
        # corresponding sessions and CAPS locations)
        sendto_longcorr_node_name = '5_sendto_longcorr'
        sendto_longcorr_node = npe.Node(
            name=sendto_longcorr_node_name,
            interface=nutil.Function(
                input_names=['in_copy2_caps',
                             'in_unpcssd_sublist',
                             'in_pcssd_capstargetlist',
                             'in_overwrite_tsv'],
                output_names=['out_unpcssd_sublist',
                              'out_pcssd_capstargetlist',
                              'out_overwrite_tsv'],
                function=utils.sendto_longcorr))

        # Connection
        # ==========
        self.connect([
            # store_reconallbase_node
            (
                self.input_node, store_reconallbase_node,
                [('subject_list', 'in_subject_list')]),
            # prepare_template_node
            (
                self.input_node, prepare_template_node,
                [('subject_list', 'in_subject'),
                 ('session_list2', 'in_session_list')]),
            # create_template_node
            (
                self.input_node, create_template_node,
                [('subject_list', 'in_subject'),
                 ('session_list2', 'in_session_list')]),
            (
                prepare_template_node, create_template_node,
                [('out_reconallbase_flags', 'in_reconallbase_flags')]),
            (
                store_reconallbase_node, create_template_node,
                [('out_workdirstore_path', 'in_workdirstore_path')]),
            # copy2caps_node
            (
                self.input_node, copy2caps_node,
                [('subject_list', 'in_subject'),
                 ('caps_target_list', 'in_caps_target')]),
            (
                store_reconallbase_node, copy2caps_node,
                [('out_workdirstore_path', 'in_subject_dir')]),
            (
                create_template_node, copy2caps_node,
                [('out_template_created', 'in_template_created')]),
            # self.output_node (copy2caps_node)
            # data for sendto_longcorr_node
            (
                self.input_node, sendto_longcorr_node,
                [('unpcssd_sublist', 'in_unpcssd_sublist'),
                 ('pcssd_capstargetlist', 'in_pcssd_capstargetlist'),
                 ('overwrite_tsv', 'in_overwrite_tsv')]),
            (
                copy2caps_node, sendto_longcorr_node,
                [('out_copy2_caps', 'in_copy2_caps')])
        ])

        return self
