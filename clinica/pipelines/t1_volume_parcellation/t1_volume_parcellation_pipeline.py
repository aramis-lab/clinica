# coding: utf8

import clinica.pipelines.engine as cpe

__author__ = "Simona Bottani"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Simona Bottani"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Simona Bottani"
__email__ = "simona.bottani@icm-institute.org"
__status__ = "Development"


class T1VolumeParcellation(cpe.Pipeline):
    """T1VolumeParcellation - Computation of mean GM concentration for a set of regions

    Args:
        input_dir: A BIDS directory.
        output_dir: An empty output directory where CAPS structured data will be written.
        subjects_sessions_list: The Subjects-Sessions list file (in .tsv format).

    Returns:
        A clinica pipeline object containing the T1VolumeParcellation pipeline.
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

        return ['file_list', 'atlas_list']

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """

        pass

    def build_input_node(self):
        """Build and connect an input node to the pipeline.
        """

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from clinica.utils.inputs import clinica_file_reader
        from clinica.utils.exceptions import ClinicaCAPSError, ClinicaException

        # Get gray matter map from t1w preprocessing (from t1-volume)
        try:
            gm_mni = clinica_file_reader(self.subjects,
                                         self.sessions,
                                         self.caps_directory,
                                         {'pattern': 't1/spm/dartel/group-' + self.parameters['group_id']
                                                     + '/*_T1w_segm-graymatter_space-Ixi549Space_modulated-'
                                                     + 'on_probability.nii*',
                                          'description': ' grey matter map in MNI space (Ixi549) with modulation',
                                          'needed_pipeline': 't1-volume'})
        except ClinicaException as e:
            final_error_str = 'Clinica faced error(s) while trying to read files in your CAPS directory.\n'
            final_error_str += str(e)
            raise ClinicaCAPSError(final_error_str)

        read_parameters_node = npe.Node(name="LoadingCLIArguments",
                                        interface=nutil.IdentityInterface(
                                            fields=self.get_input_fields(),
                                            mandatory_inputs=True))

        read_parameters_node.inputs.file_list = gm_mni
        read_parameters_node.inputs.atlas_list = self.parameters['atlases']

        self.connect([(read_parameters_node, self.input_node, [('file_list', 'file_list')]),
                      (read_parameters_node, self.input_node, [('atlas_list', 'atlas_list')])
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
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio
        from ..t1_volume_parcellation import t1_volume_parcellation_utils as parcellation_utils

        atlas_stats_node = npe.MapNode(nutil.Function(input_names=['in_image',
                                                                   'atlas_list'],
                                                      output_names=['atlas_statistics'],
                                                      function=parcellation_utils.atlas_statistics),
                                       name='atlas_stats_node',
                                       iterfield=['in_image'])
        outputnode = npe.Node(nutil.IdentityInterface(fields=['atlas_statistics']),
                              name='outputnode',
                              mandatory_inputs=True)

        datasink = npe.Node(nio.DataSink(),
                            name='datasink')

        datasink.inputs.base_directory = self.caps_directory
        datasink.inputs.parameterization = True
        datasink.inputs.regexp_substitutions = [
            (r'(.*)(atlas_statistics)/.*/(sub-(.*)_ses-(.*)_T1.*)$',
             r'\1/subjects/sub-\4/ses-\5/t1/spm/dartel/group-' + self.parameters['group_id'] + r'/\2/\3')]

        # Connection
        # ==========
        self.connect([
            (self.input_node,      atlas_stats_node,    [('file_list',    'in_image')]),
            (self.input_node,      atlas_stats_node,    [('atlas_list',    'atlas_list')]),
            (atlas_stats_node,     outputnode,          [('atlas_statistics',  'atlas_statistics')]),
            (outputnode,           datasink,            [('atlas_statistics', 'atlas_statistics')])
        ])
