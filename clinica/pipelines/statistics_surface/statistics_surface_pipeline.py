# coding: utf8

import clinica.pipelines.engine as cpe


class StatisticsSurface(cpe.Pipeline):
    """
    Based on the Matlab toolbox [SurfStat](http://www.math.mcgill.ca/keith/surfstat/), which performs statistical
    analyses of univariate and multivariate surface and volumetric data using the generalized linear model (GLM),
    this pipeline performs analyses including group comparison and correlation with surface-based features e.g.
    cortical thickness from t1-freesurfer or map of activity from PET data from pet-surface pipeline.

    TODO: Refactor StatisticsSurface
        [ ] build_input_node
            [ ] Remove current read_parameters_node
            [ ] With the help of statistics_surface_utils.py::check_inputs, use new clinica_file_reader function
                to check and extract surface-based features
            [ ] Delete statistics_surface_utils.py::check_inputs function
            [ ] Move statistics_surface_cli.py checks of input data in this method
            [ ] Handle overwrite case
            [ ] Display participants, covariates and info regarding the GLM.
        [ ] build_core_nodes
            [ ] Use surfstat.inputs.full_width_at_half_maximum = self.parameter['full_width_at_half_maximum']
                instead of connecting read_parameters_node.inputs.full_width_at_half_maximum
                to surfstat.inputs.full_width_at_half_maximum
            [ ] Repeat for other keys
            [ ] Use working directory
        [ ] build_output_node
            [ ] Remove path_to_matscript and freesurfer_home: it should be set in runmatlab function
            [ ] Copy results from <WD> to <CAPS>
        [ ] Clean/adapt statistics_surface_utils.py

    Args:
        caps_directory: str, the output folder of recon-all which will contain the result files: ?h.thickness.fwhm**.mgh.
        tsv_file: str, Path to the tsv containing the information for GLM.
        design_matrix: str, the linear model that fits into the GLM, for example '1+group'.
        contrast: string, the contrast matrix for GLM, if the factor you choose is categorized variable, clinica_surfstat will create two contrasts,
                  for example, contrast = 'Label', this will create contrastpos = Label.AD - Label.CN, contrastneg = Label.CN - Label.AD; if the fac-
                  tory that you choose is a continuous factor, clinica_surfstat will just create one contrast, for example, contrast = 'Age', but note,
                  the string name that you choose should be exactly the same with the columns names in your subjects_visits_tsv.
        str_format: string, the str_format which uses to read your tsv file, the typy of the string should corresponds exactly with the columns in the tsv file.
            Defaut parameters, we set these parameters to be some default values, but you can also set it by yourself:
        group_label: current group name for this analysis
        glm_type: based on the hypothesis, you should define one of the glm types, "group_comparison", "correlation"
        full_width_at_half_maximum: fwhm for the surface smoothing, default is 20, integer.
        threshold_uncorrected_pvalue: threshold to display the uncorrected Pvalue, float, default is 0.001.
        threshold_corrected_pvalue: the threshold to display the corrected cluster, default is 0.05, float.
        cluster_threshold: threshold to define a cluster in the process of cluster-wise correction, default is 0.001, float.
        working_directory: define where to put the infomartion of the nipype workflow.
        n_procs: define how many cores to run this workflow.

    Returns:
        A clinica pipeline object containing the StatisticsSurface pipeline.
    """
    def check_pipeline_parameters(self):
        """Check pipeline parameters."""
        from .statistics_surface_utils import get_t1_freesurfer_custom_file
        from clinica.utils.exceptions import ClinicaException
        from clinica.utils.group import check_group_label

        if 'custom_file' not in self.parameters.keys():
            self.parameters['custom_file'] = get_t1_freesurfer_custom_file()
        if 'feature_label' not in self.parameters.keys():
            self.parameters['feature_label'] = 'ct',
        if 'full_width_at_half_maximum' not in self.parameters.keys():
            self.parameters['full_width_at_half_maximum'] = 20
        if 'threshold_uncorrected_pvalue' not in self.parameters.keys():
            self.parameters['threshold_uncorrected_pvalue'] = 0.001
        if 'threshold_corrected_pvalue' not in self.parameters.keys():
            self.parameters['threshold_corrected_pvalue'] = 0.05,
        if 'cluster_threshold' not in self.parameters.keys():
            self.parameters['cluster_threshold'] = 0.001,

        check_group_label(self.parameters['group_label'])
        if self.parameters['glm_type'] not in ['group_comparison', 'correlation']:
            raise ClinicaException("The glm_type you specified is wrong: it should be group_comparison or "
                                   "correlation (given value: %s)." % self.parameters['glm_type'])
        if self.parameters['full_width_at_half_maximum'] not in [0, 5, 10, 15, 20]:
            raise ClinicaException(
                "FWHM for the surface smoothing you specified is wrong: it should be 0, 5, 10, 15 or 20 "
                "(given value: %s)." % self.parameters['full_width_at_half_maximum'])
        if self.parameters['threshold_uncorrected_pvalue'] < 0 or self.parameters['threshold_uncorrected_pvalue'] > 1:
            raise ClinicaException("Uncorrected p-value threshold should be a lower than 1 "
                                   "(given value: %s)." % self.parameters['threshold_uncorrected_pvalue'])
        if self.parameters['threshold_corrected_pvalue'] < 0 or self.parameters['threshold_corrected_pvalue'] > 1:
            raise ClinicaException("Corrected p-value threshold should be between 0 and 1 "
                                   "(given value: %s)." % self.parameters['threshold_corrected_pvalue'])
        if self.parameters['cluster_threshold'] < 0 or self.parameters['cluster_threshold'] > 1:
            raise ClinicaException("Cluster threshold should be between 0 and 1 "
                                   "(given value: %s)." % self.parameters['cluster_threshold'])

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file.
        """
        pass

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipelines.

        Returns:
            A list of (string) input fields name.
        """

        return ['design_matrix',
                'contrast',
                'str_format',
                'group_label',
                'glm_type',
                'surface_file',
                'full_width_at_half_maximum',
                'threshold_uncorrected_pvalue',
                'threshold_corrected_pvalue',
                'cluster_threshold',
                'feature_label']

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipelines.

        Returns:
            A list of (string) output fields name.
        """

        return []

    def build_input_node(self):
        """Build and connect an input node to the pipelines.
        """

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        read_parameters_node = npe.Node(name="LoadingCLIArguments",
                                        interface=nutil.IdentityInterface(
                                            fields=self.get_input_fields(),
                                            mandatory_inputs=True))
        read_parameters_node.inputs.design_matrix = self.parameters['design_matrix']
        read_parameters_node.inputs.contrast = self.parameters['contrast']
        read_parameters_node.inputs.str_format = self.parameters['str_format']
        read_parameters_node.inputs.group_label = self.parameters['group_label']
        read_parameters_node.inputs.glm_type = self.parameters['glm_type']
        read_parameters_node.inputs.surface_file = self.parameters['custom_file']
        read_parameters_node.inputs.full_width_at_half_maximum = self.parameters['full_width_at_half_maximum']
        read_parameters_node.inputs.threshold_uncorrected_pvalue = self.parameters['threshold_uncorrected_pvalue']
        read_parameters_node.inputs.threshold_corrected_pvalue = self.parameters['threshold_corrected_pvalue']
        read_parameters_node.inputs.cluster_threshold = self.parameters['cluster_threshold']
        read_parameters_node.inputs.feature_label = self.parameters['feature_label']

        self.connect([
            (read_parameters_node, self.input_node, [('design_matrix',                'design_matrix')]),  # noqa
            (read_parameters_node, self.input_node, [('surface_file',                 'surface_file')]),  # noqa
            (read_parameters_node, self.input_node, [('contrast',                     'contrast')]),  # noqa
            (read_parameters_node, self.input_node, [('str_format',                   'str_format')]),  # noqa
            (read_parameters_node, self.input_node, [('group_label',                  'group_label')]),  # noqa
            (read_parameters_node, self.input_node, [('glm_type',                     'glm_type')]),  # noqa
            (read_parameters_node, self.input_node, [('full_width_at_half_maximum',   'full_width_at_half_maximum')]),  # noqa
            (read_parameters_node, self.input_node, [('threshold_uncorrected_pvalue', 'threshold_uncorrected_pvalue')]),  # noqa
            (read_parameters_node, self.input_node, [('threshold_corrected_pvalue',   'threshold_corrected_pvalue')]),  # noqa
            (read_parameters_node, self.input_node, [('cluster_threshold',            'cluster_threshold')]),  # noqa
            (read_parameters_node, self.input_node, [('feature_label',                'feature_label')])
        ])

    def build_output_node(self):
        """Build and connect an output node to the pipelines.
        """

        pass

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipelines.
        """

        import clinica.pipelines.statistics_surface.statistics_surface_utils as utils
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from nipype.interfaces.io import JSONFileSink

        # Node to fetch the input variables
        data_prep = npe.Node(name='inputnode',
                             interface=nutil.Function(
                                 input_names=['input_directory', 'subjects_visits_tsv', 'group_label', 'glm_type'],
                                 output_names=['path_to_matscript', 'surfstat_input_dir', 'output_directory',
                                               'freesurfer_home', 'out_json'],
                                 function=utils.prepare_data))
        data_prep.inputs.input_directory = self.caps_directory
        data_prep.inputs.subjects_visits_tsv = self.tsv_file

        # Node to wrap the SurfStat matlab script
        surfstat = npe.Node(name='surfstat',
                            interface=nutil.Function(
                                input_names=['input_directory',
                                             'output_directory',
                                             'subjects_visits_tsv',
                                             'design_matrix',
                                             'contrast',
                                             'str_format',
                                             'glm_type',
                                             'group_label',
                                             'freesurfer_home',
                                             'surface_file',
                                             'path_to_matscript',
                                             'full_width_at_half_maximum',
                                             'threshold_uncorrected_pvalue',
                                             'threshold_corrected_pvalue',
                                             'cluster_threshold',
                                             'feature_label'],
                                output_names=['out_images'],
                                function=utils.runmatlab))
        surfstat.inputs.subjects_visits_tsv = self.tsv_file

        # Node to create the dictionary for JSONFileSink
        jsondict = npe.Node(name='Jsondict',
                            interface=nutil.Function(
                                input_names=['glm_type', 'design_matrix', 'str_format', 'contrast', 'group_label', 'full_width_at_half_maximum',
                                             'threshold_uncorrected_pvalue', 'threshold_corrected_pvalue', 'cluster_threshold'],
                                output_names=['json_dict'],
                                function=utils.json_dict_create))

        # Node to write the GLM information into a JSON file
        jsonsink = npe.Node(JSONFileSink(input_names=['out_file']), name='jsonsinker')

        # Connection
        # ==========
        self.connect([
            (self.input_node, data_prep, [('group_label', 'group_label')]),  # noqa
            (self.input_node, data_prep, [('glm_type', 'glm_type')]),  # noqa
            (self.input_node, surfstat, [('design_matrix', 'design_matrix')]),  # noqa
            (self.input_node, surfstat, [('contrast', 'contrast')]),  # noqa
            (self.input_node, surfstat, [('str_format', 'str_format')]),  # noqa
            (self.input_node, surfstat, [('glm_type', 'glm_type')]),  # noqa
            (self.input_node, surfstat, [('group_label', 'group_label')]),  # noqa
            (self.input_node, surfstat, [('full_width_at_half_maximum', 'full_width_at_half_maximum')]),  # noqa
            (self.input_node, surfstat, [('threshold_uncorrected_pvalue', 'threshold_uncorrected_pvalue')]),  # noqa
            (self.input_node, surfstat, [('threshold_corrected_pvalue', 'threshold_corrected_pvalue')]),  # noqa
            (self.input_node, surfstat, [('cluster_threshold', 'cluster_threshold')]),  # noqa
            (self.input_node, surfstat, [('surface_file', 'surface_file')]),  # noqa
            (self.input_node, surfstat, [('feature_label', 'feature_label')]),
            (data_prep, surfstat, [('surfstat_input_dir', 'input_directory')]),  # noqa
            (data_prep, surfstat, [('path_to_matscript', 'path_to_matscript')]),  # noqa
            (data_prep, surfstat, [('output_directory', 'output_directory')]),  # noqa
            (data_prep, surfstat, [('freesurfer_home', 'freesurfer_home')]),  # noqa
            (self.input_node, jsondict, [('glm_type', 'glm_type')]),  # noqa
            (self.input_node, jsondict, [('design_matrix', 'design_matrix')]),  # noqa
            (self.input_node, jsondict, [('str_format', 'str_format')]),  # noqa
            (self.input_node, jsondict, [('contrast', 'contrast')]),  # noqa
            (self.input_node, jsondict, [('group_label', 'group_label')]),  # noqa
            (self.input_node, jsondict, [('full_width_at_half_maximum', 'full_width_at_half_maximum')]),  # noqa
            (self.input_node, jsondict, [('threshold_uncorrected_pvalue', 'threshold_uncorrected_pvalue')]),  # noqa
            (self.input_node, jsondict, [('threshold_corrected_pvalue', 'threshold_corrected_pvalue')]),  # noqa
            (self.input_node, jsondict, [('cluster_threshold', 'cluster_threshold')]),  # noqa
            (data_prep, jsonsink, [('out_json', 'out_file')]),  # noqa
            (jsondict, jsonsink, [('json_dict', 'in_dict')]),  # noqa
        ])
