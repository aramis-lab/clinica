# coding: utf8

__author__ = "Junhao Wen"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Junhao Wen", "Alexandre Routier"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Junhao Wen"
__email__ = "Junhao.Wen@inria.fr"
__status__ = "Development"

import clinica.pipelines.engine as cpe


class T1FreeSurferCrossSectional(cpe.Pipeline):
    """Performs cortical surface extraction, segmentation of subcortical
    structures, cortical thickness estimation, spatial normalization onto the
    FreeSurfer surface template (FsAverage), and parcellation of cortical
    regions based on different atlases using the FreeSurfer recon-all command.
    It also generate TSV files containing a summary of the regional statistics
    (e.g. regional volume, mean cortical thickness) to ease subsequent
    statistical analysis.

    Args:
        bids_directory: Path to the BIDS directory.
        caps_directory: Path to the CAPS directory.
        tsv: TSV file containing the subjects with their sessions
        ras: additional flags for recon-all command line, default will be -qcache
        wd: Temporary directory to store pipelines intermediate results
        np: Number of cores used to run in parallel

    Returns:
        A clinica pipeline object containing the T1FreeSurferCrossSectional pipeline.
    """

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file.
        """
        pass

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipelines.

        Returns:
            A list of (string) input fields name.
        """
        return ['recon_all_args', 'subject_list', 'session_list', 'subject_id', 'subject_dir', 'anat_t1']

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipelines.

        Returns:
            A list of (string) output fields name.
        """
        return ['subject_id']

    def build_input_node(self):
        """Build and connect an input node to the pipelines.
        """
        import clinica.pipelines.t1_freesurfer_cross_sectional.t1_freesurfer_cross_sectional_utils as utils
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from clinica.utils.stream import cprint

        # This node is supposedly used to load BIDS inputs when this pipelines is
        # not already connected to the output of a previous Clinica pipelines.
        # For the purpose of the example, we simply read input arguments given
        # by the command line interface and transmitted here through the
        # `self.parameters` dictionary and pass it to the `self.input_node` to
        # further by used as input of the core nodes.

        # An IdentityInterface to get the input fields for the pipelines
        read_parameters_node = npe.Node(name="LoadingCLIArguments",
                                        interface=nutil.IdentityInterface(
                                            fields=self.get_input_fields(),
                                            mandatory_inputs=True))
        read_parameters_node.inputs.recon_all_args = self.parameters['recon_all_args']

        # check if a specific subject has already run recon-all, and return the related lists in the tsv file
        dataprepr = npe.Node(name='dataprepr',
                             interface=nutil.Function(
                                 input_names=['output_dir', 'subject_list', 'session_list'],
                                 output_names=['nouse0', 'nouse1', 'subject_dir', 'subject_id',
                                               'subject_list', 'session_list'],
                                 function=utils.get_dirs_check_reconalled))
        dataprepr.inputs.output_dir = self.caps_directory
        dataprepr.inputs.subject_list = self.subjects
        dataprepr.inputs.session_list = self.sessions

        # use pybids to grab the anat_t1 images based on the tsv file
        datagrabbernode = npe.Node(name='datagrabbernode',
                                   interface=nutil.Function(
                                       function=utils.bids_datagrabber,
                                       input_names=['input_dir', 'subject_list', 'session_list'],
                                       output_names=['anat_t1']))
        datagrabbernode.inputs.input_dir = self.bids_directory

        self.connect([
            (read_parameters_node, self.input_node, [('recon_all_args', 'recon_all_args')]),
            (dataprepr, datagrabbernode, [('subject_list', 'subject_list')]),
            (dataprepr, datagrabbernode, [('session_list', 'session_list')]),
            (dataprepr, self.input_node, [('subject_list', 'subject_list')]),
            (dataprepr, self.input_node, [('session_list', 'session_list')]),
            (dataprepr, self.input_node, [('subject_id', 'subject_id')]),
            (dataprepr, self.input_node, [('subject_dir', 'subject_dir')]),
            (datagrabbernode, self.input_node, [('anat_t1', 'anat_t1')]),
        ])

        cprint('The pipeline will last approximately 10 hours per image.')

    def build_output_node(self):
        """Build and connect an output node to the pipelines.
        """
        pass

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipelines.
        """
        import clinica.pipelines.t1_freesurfer_cross_sectional.t1_freesurfer_cross_sectional_utils as utils
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from nipype.interfaces.freesurfer.preprocess import ReconAll

        # MapNode to check out if we need -cw256 for every subject, and -qcache is default for every subject.
        flagnode = npe.MapNode(name='flagnode',
                               iterfield=['t1_list'],
                               interface=nutil.Function(
                                   input_names=['t1_list', 'recon_all_args'],
                                   output_names=['output_flags'],
                                   function=utils.check_fov))

        # MapNode to transfer every subject's flag to string.
        create_flags = npe.MapNode(interface=nutil.Function(
            input_names=['input_flags'],
            output_names=['output_str'],
            function=utils.create_flags_str),
            name='create_flags_string',
            iterfield=['input_flags'])

        # MapNode to implement recon-all.
        recon_all = npe.MapNode(interface=ReconAll(),
                                name='recon_all',
                                iterfield=['subject_id', 'T1_files', 'subjects_dir', 'flags'])
        recon_all.inputs.directive = 'all'

        # MapNode to create the statistical tsv file for each subject
        tsvmapnode = npe.MapNode(name='tsvmapnode',
                                 iterfield=['subject_id'],
                                 interface=nutil.Function(
                                     input_names=['subject_id', 'output_dir'],
                                     output_names=[],
                                     function=utils.write_tsv_files,
                                     imports=['import os', 'import errno']))
        tsvmapnode.inputs.output_dir = self.caps_directory

        print_begin_message = npe.MapNode(
            interface=nutil.Function(
                input_names=['subject_id'],
                function=utils.print_begin_pipeline),
            iterfield='subject_id',
            name='Write-Begin_Message')

        print_end_message = npe.MapNode(
            interface=nutil.Function(
                input_names=['subject_id', 'final_file'],
                function=utils.print_end_pipeline),
            iterfield='subject_id',
            name='Write-End_Message')

        # Connection
        # ==========
        self.connect([
            (self.input_node, print_begin_message, [('subject_id', 'subject_id')]),
            (self.input_node, flagnode, [('recon_all_args', 'recon_all_args')]),
            (self.input_node, flagnode, [('anat_t1',        't1_list')]),
            (flagnode, create_flags, [('output_flags', 'input_flags')]),
            (create_flags,    recon_all, [('output_str',  'flags')]),
            (self.input_node, recon_all, [('subject_dir', 'subjects_dir')]),
            (self.input_node, recon_all, [('subject_id',  'subject_id')]),
            (self.input_node, recon_all, [('anat_t1',     'T1_files')]),
            #            (self.input_node, lognode, [('subject_list', 'subject_list')]),
            #            (self.input_node, lognode, [('session_list', 'session_list')]),
            (recon_all, tsvmapnode, [('subject_id', 'subject_id')]),
            #            (recon_all, lognode, [('subject_id', 'subject_id')]),
            (recon_all, self.output_node, [('subject_id', 'subject_id')]),
            (self.input_node, print_end_message, [('subject_id', 'subject_id')]),
            (recon_all,       print_end_message, [('subject_id', 'final_file')]),
        ])

        return self
