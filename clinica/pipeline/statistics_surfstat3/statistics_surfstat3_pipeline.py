"""Statistics Surfstat3 - Clinica Pipeline.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details: https://gitlab.icm-institute.org/aramis/clinica/wikis/docs/InteractingWithClinica.
"""

# WARNING: Don't put any import statement here except if it's absolutly
# necessary. Put it *inside* the different methods.
# Otherwise it will slow down the dynamic loading of the pipelines list by the
# command line tool.
import clinica.pipeline.engine as cpe


class StatisticsSurfstat3(cpe.Pipeline):
    """Statistics Surfstat3 SHORT DESCRIPTION.

    Warnings:
        - A WARNING.

    Todos:
        - [x] A FILLED TODO ITEM.
        - [ ] AN ON-GOING TODO ITEM.

    Args:
        input_dir: A BIDS directory.
        output_dir: An empty output directory where CAPS structured data will be written.
        subjects_sessions_list: The Subjects-Sessions list file (in .tsv format).

    Returns:
        A clinica pipeline object containing the Statistics Surfstat3 pipeline.

    Raises:


    Example:
        >>> from statistics_surfstat3 import StatisticsSurfstat3
        >>> pipeline = FMRIPreprocessing('~/MYDATASET_BIDS', '~/MYDATASET_CAPS')
        >>> pipeline.parameters = {
        >>>     # ...
        >>> }
        >>> pipeline.base_dir = '/tmp/'
        >>> pipeline.run()
    """


    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file.
        """
        pass


    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Returns:
            A list of (string) input fields name.
        """

        return ['design_matrix', 'contrast', 'str_format', 'group_label', 'glm_type', 'full_width_at_half_maximum', 'threshold_uncorrected_pvalue', 'threshold_corrected_pvalue', 'cluster_threshold'] # Fill here the list


    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """

        return [] # Fill here the list


    def build_input_node(self):
        """Build and connect an input node to the pipeline.
        """

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        # This node is supposedly used to load BIDS inputs when this pipeline is
        # not already connected to the output of a previous Clinica pipeline.
        # For the purpose of the example, we simply read input arguments given
        # by the command line interface and transmitted here through the
        # `self.parameters` dictionary and pass it to the `self.input_node` to
        # further by used as input of the core nodes.
        read_parameters_node = npe.Node(name="LoadingCLIArguments",
                                        interface=nutil.IdentityInterface(
                                            fields=self.get_input_fields(),
                                            mandatory_inputs=True))
        read_parameters_node.inputs.design_matrix = self.parameters['design_matrix']
        read_parameters_node.inputs.contrast = self.parameters['contrast']
        read_parameters_node.inputs.str_format = self.parameters['str_format']
        read_parameters_node.inputs.group_label = self.parameters['group_label']
        read_parameters_node.inputs.glm_type = self.parameters['glm_type']
        read_parameters_node.inputs.full_width_at_half_maximum = self.parameters['full_width_at_half_maximum']
        read_parameters_node.inputs.threshold_uncorrected_pvalue = self.parameters['threshold_uncorrected_pvalue']
        read_parameters_node.inputs.threshold_corrected_pvalue = self.parameters['threshold_corrected_pvalue']
        read_parameters_node.inputs.cluster_threshold = self.parameters['cluster_threshold']

        self.connect([
            (read_parameters_node,      self.input_node,    [('design_matrix',    'design_matrix')]),
            (read_parameters_node,      self.input_node,    [('contrast',    'contrast')]),
            (read_parameters_node,      self.input_node,    [('str_format',    'str_format')]),
            (read_parameters_node,      self.input_node,    [('group_label',    'group_label')]),
            (read_parameters_node,      self.input_node,    [('glm_type',    'glm_type')]),
            (read_parameters_node,      self.input_node,    [('full_width_at_half_maximum',    'full_width_at_half_maximum')]),
            (read_parameters_node,      self.input_node,    [('threshold_uncorrected_pvalue',    'threshold_uncorrected_pvalue')]),
            (read_parameters_node,      self.input_node,    [('threshold_corrected_pvalue',    'threshold_corrected_pvalue')]),
            (read_parameters_node,      self.input_node,    [('cluster_threshold',    'cluster_threshold')]),
        ])


    def build_output_node(self):
        """Build and connect an output node to the pipeline.
        """

        # In the same idea as the input node, this output node is supposedly
        # used to write the output fields in a CAPS. It should be executed only
        # if this pipeline output is not already connected to a next Clinica
        # pipeline.

        pass


    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline.
        """

        import statistics_surfstat3_utils as utils
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from nipype.interfaces.io import JSONFileSink

        # Node to fetch the input vars.
        data_prep = npe.Node(name='inputnode',
                             interface=nutil.Function(
                                 input_names=['input_directory', 'subjects_visits_tsv', 'group_label', 'glm_type'],
                                 output_names=['path_to_matscript', 'surfstat_input_dir', 'output_directory',
                                               'freesurfer_home', 'out_json'],
                                 function=utils.data_prep))
        data_prep.inputs.input_directory = self.caps_directory
        data_prep.inputs.subjects_visits_tsv = self.tsv_file

        # Node to wrap the surfstat matlab script.
        surfstat = npe.Node(name='surfstat',
                            interface=nutil.Function(
                                input_names=['input_directory', 'output_directory', 'subjects_visits_tsv',
                                             'design_matrix',
                                             'contrast', 'str_format', 'glm_type', 'group_label',
                                             'freesurfer_home', 'path_to_matscript',
                                             'full_width_at_half_maximum', 'threshold_uncorrected_pvalue',
                                             'threshold_corrected_pvalue', 'cluster_threshold'],
                                output_names=['out_images'],
                                function=utils.runmatlab))
        surfstat.inputs.subjects_visits_tsv = self.tsv_file

        ### Node to create the dic for JSONFileSink
        jsondict = npe.Node(name='Jsondict',
                            interface=nutil.Function(
                                input_names=['glm_type', 'design_matrix', 'str_format', 'contrast', 'group_label', 'full_width_at_half_maximum',
                                             'threshold_uncorrected_pvalue', 'threshold_corrected_pvalue', 'cluster_threshold'],
                                output_names=['json_dict'],
                                function=utils.json_dict_create))

        # Node to write the GLM infor into a json file
        jsonsink = npe.Node(JSONFileSink(input_names=['out_file']), name='jsonsinker')

        # Connection
        # ==========
        self.connect([
            (self.input_node, data_prep, [('group_label', 'group_label')]),
            (self.input_node, data_prep, [('glm_type', 'glm_type')]),
            (self.input_node, surfstat, [('design_matrix', 'design_matrix')]),
            (self.input_node, surfstat, [('contrast', 'contrast')]),
            (self.input_node, surfstat, [('str_format', 'str_format')]),
            (self.input_node, surfstat, [('glm_type', 'glm_type')]),
            (self.input_node, surfstat, [('group_label', 'group_label')]),
            (self.input_node, surfstat, [('full_width_at_half_maximum', 'full_width_at_half_maximum')]),
            (self.input_node, surfstat, [('threshold_uncorrected_pvalue', 'threshold_uncorrected_pvalue')]),
            (self.input_node, surfstat, [('threshold_corrected_pvalue', 'threshold_corrected_pvalue')]),
            (self.input_node, surfstat, [('cluster_threshold', 'cluster_threshold')]),
            (data_prep, surfstat, [('surfstat_input_dir', 'input_directory')]),
            (data_prep, surfstat, [('path_to_matscript', 'path_to_matscript')]),
            (data_prep, surfstat, [('output_directory', 'output_directory')]),
            (data_prep, surfstat, [('freesurfer_home', 'freesurfer_home')]),
            (self.input_node, jsondict, [('glm_type', 'glm_type')]),
            (self.input_node, jsondict, [('design_matrix', 'design_matrix')]),
            (self.input_node, jsondict, [('str_format', 'str_format')]),
            (self.input_node, jsondict, [('contrast', 'contrast')]),
            (self.input_node, jsondict, [('group_label', 'group_label')]),
            (self.input_node, jsondict, [('full_width_at_half_maximum', 'full_width_at_half_maximum')]),
            (self.input_node, jsondict, [('threshold_uncorrected_pvalue', 'threshold_uncorrected_pvalue')]),
            (self.input_node, jsondict, [('threshold_corrected_pvalue', 'threshold_corrected_pvalue')]),
            (self.input_node, jsondict, [('cluster_threshold', 'cluster_threshold')]),
            (data_prep, jsonsink, [('out_json', 'out_file')]),
            (jsondict, jsonsink, [('json_dict', 'in_dict')]),
        ])
