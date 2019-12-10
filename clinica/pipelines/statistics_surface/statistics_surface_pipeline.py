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
            [X] Remove current read_parameters_node
            [X] With the help of statistics_surface_utils.py::check_inputs, use new clinica_file_reader function
                to check and extract surface-based features
            [X] Delete statistics_surface_utils.py::check_inputs function
            [X] Move statistics_surface_cli.py checks of input data in this method
            [ ] Handle overwrite case
            [ ] Display participants, covariates and info regarding the GLM.
        [ ] build_core_nodes
            [X] Use surfstat.inputs.full_width_at_half_maximum = self.parameter['full_width_at_half_maximum']
                instead of connecting read_parameters_node.inputs.full_width_at_half_maximum
                to surfstat.inputs.full_width_at_half_maximum
            [X] Repeat for other keys
            [ ] Use working directory
        [ ] build_output_node
            [X] Remove path_to_matscript and freesurfer_home: it should be set in runmatlab function
            [ ] Copy results from <WD> to <CAPS>
        [ ] Clean/adapt statistics_surface_utils.py

    Note:
        The `tsv_file` attribute is overloaded for this pipeline. It must contain a list of subjects
        with their sessions and all the covariates and factors needed for the GLM.

        Pipeline parameters are explained in StatisticsSurfaceCLI.define_options()

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

        return []

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipelines.

        Returns:
            A list of (string) output fields name.
        """

        return []

    def build_input_node(self):
        """Build and connect an input node to the pipelines.
        """
        import os
        from clinica.utils.inputs import clinica_file_reader
        from clinica.utils.exceptions import ClinicaException
        from clinica.utils.stream import cprint

        # Check if the group label has been existed, if yes, give an error to the users
        # Note(AR): if the user wants to compare Cortical Thickness measure with PET measure
        # using the group_id, Clinica won't allow it.
        # TODO: Modify this behaviour
        if os.path.exists(os.path.join(self.caps_directory, 'groups', 'group-' + self.parameters['group_label'])):
            error_message = ('Group ID %s already exists, please choose another one or delete the existing folder and '
                             'also the working directory and rerun the pipeline') % self.parameters['group_label']
            raise ClinicaException(error_message)

        # Check input files
        all_errors = []
        # clinica_files_reader expects regexp to start at subjects/ so sub-*/ses-*/ is removed here
        pattern_hemisphere = self.parameters['custom_file'].replace(
            '@subject', 'sub-*').replace(
            '@session', 'ses-*').replace(
            '@fwhm', str(self.parameters['full_width_at_half_maximum'])).replace(
            'sub-*/ses-*/', ''
        )
        # Files on left hemisphere
        lh_surface_based_info = {
            'pattern': pattern_hemisphere.replace('@hemi', 'lh'),
            'description': 'surface-based features on left hemisphere at FWHM = %s' %
                           self.parameters['full_width_at_half_maximum'],
        }
        try:
            clinica_file_reader(self.subjects, self.sessions, self.caps_directory, lh_surface_based_info)
        except ClinicaException as e:
            all_errors.append(e)
        rh_surface_based_info = {
            'pattern': pattern_hemisphere.replace('@hemi', 'rh'),
            'description': 'surface-based features on right hemisphere at FWHM = %s' %
                           self.parameters['full_width_at_half_maximum'],
        }
        try:
            clinica_file_reader(self.subjects, self.sessions, self.caps_directory, rh_surface_based_info)
        except ClinicaException as e:
            all_errors.append(e)
        # Raise all errors if something happened
        if len(all_errors) > 0:
            error_message = 'Clinica faced errors while trying to read files in your CAPS directory.\n'
            for msg in all_errors:
                error_message += str(msg)
            raise RuntimeError(error_message)

        # Print GLM information
        cprint("Parameters used for this pipeline:")
        cprint(self.parameters)

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
                                 output_names=['surfstat_input_dir', 'output_directory', 'out_json'],
                                 function=utils.prepare_data))
        data_prep.inputs.input_directory = self.caps_directory
        data_prep.inputs.subjects_visits_tsv = self.tsv_file
        data_prep.inputs.group_label = self.parameters['group_label']
        data_prep.inputs.glm_type = self.parameters['glm_type']

        # Node to wrap the SurfStat matlab script
        surfstat = npe.Node(name='surfstat',
                            interface=nutil.Function(
                                input_names=['input_directory',
                                             'output_directory',
                                             'subjects_visits_tsv',
                                             'pipeline_parameters',
                                             ],
                                output_names=['out_images'],
                                function=utils.run_matlab))
        surfstat.inputs.subjects_visits_tsv = self.tsv_file
        surfstat.inputs.pipeline_parameters = self.parameters

        # Node to create the dictionary for JSONFileSink
        json_dict = npe.Node(name='Jsondict',
                             interface=nutil.Function(
                                 input_names=['pipeline_parameters'],
                                 output_names=['json_dict'],
                                 function=utils.create_glm_info_dictionary))
        json_dict.inputs.pipeline_parameters = self.parameters

        # Node to write the GLM information into a JSON file
        json_datasink = npe.Node(JSONFileSink(
            input_names=['out_file']),
            name='json_datasink')

        # Connection
        # ==========
        self.connect([
            (data_prep, surfstat, [('surfstat_input_dir', 'input_directory')]),
            (data_prep, surfstat, [('output_directory', 'output_directory')]),
            (data_prep, json_datasink, [('out_json', 'out_file')]),
            (json_dict, json_datasink, [('json_dict', 'in_dict')]),
        ])
