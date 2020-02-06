"""Statistics_Volume_Correction - Clinica Pipeline.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details: https://gitlab.icm-institute.org/aramislab/clinica/wikis/docs/InteractingWithClinica.
"""

import clinica.pipelines.engine as cpe


class StatisticsVolumeCorrection(cpe.Pipeline):
    """Statistics_Volume_Correction


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

        return ['t_map']


    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """

        return []


    def build_input_node(self):
        """Build and connect an input node to the pipeline.
        """

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from clinica.utils.inputs import clinica_group_reader
        from clinica.utils.exceptions import ClinicaException

        all_errors = []
        try:
            t_map = clinica_group_reader(self.caps_directory, {'pattern': self.parameters['t_map'] + '*',
                                                               'description': 'statistics t map',
                                                               'needed_pipeline': 'statistics-volume'})
        except ClinicaException as e:
            all_errors.append(e)

        read_parameters_node = npe.Node(name="LoadingCLIArguments",
                                        interface=nutil.IdentityInterface(
                                            fields=self.get_input_fields(),
                                            mandatory_inputs=True))
        read_parameters_node.inputs.t_map = t_map

        self.connect([
            (read_parameters_node,      self.input_node,    [('t_map', 't_map')])
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

        import statistics_volume_correction_utils as utils
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        peak_correction_FWE = npe.Node(name='peak_correction_FWE',
                                       interface=nutil.Function(
                                           input_names=['t_map', 't_threshold'],
                                           output_names=['output'],
                                           function=utils.peak_correction))
        peak_correction_FDR = peak_correction_FWE.clone(name='peak_correction_FDR')

        cluster_correction_FWE = npe.Node(name='cluster_correction_FWE',
                                      interface=nutil.Function(
                                          input_names=['t_map', 't_thresh', 'c_thresh'],
                                          output_names=['output'],
                                          function=utils.cluster_correction))
        cluster_correction_FDR = cluster_correction_FWE.clone(name='cluster_correction_FDR')


        # Connection
        # ==========
        self.connect([
            # STEP 1
            (self.input_node,      node1,    [('hello_word',    'in_hello_word')]),
            # STEP 2
            (self.input_node,      node2,    [('hello_word',    'in_hello_word')])
        ])