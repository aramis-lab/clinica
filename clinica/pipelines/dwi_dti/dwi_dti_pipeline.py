# coding: utf8

import clinica.pipelines.engine as cpe
from nipype import config

# Use hash instead of parameters for iterables folder names
# Otherwise path will be too long and generate OSError
cfg = dict(execution={'parameterize_dirs': False})
config.update_config(cfg)


class DwiDti(cpe.Pipeline):
    """DTI-based processing of DWI datasets.

    Returns:
        A clinica pipeline object containing the DwiDti pipeline.

    Raises:

    """

    def check_pipeline_parameters(self):
        """Check pipeline parameters."""
        pass

    def check_custom_dependencies(self): pass

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipelines.

        Returns:
            A list of (string) input fields name.
        """
        return ['preproc_dwi', 'preproc_bvec', 'preproc_bval', 'b0_mask']

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipelines.

        Returns:
            A list of (string) output fields name.
        """
        output_list_dti = [
            'dti', 'fa', 'md', 'ad', 'rd', 'decfa',
            'registered_fa', 'registered_md', 'registered_ad', 'registered_rd',
            'statistics_fa', 'statistics_md', 'statistics_ad', 'statistics_rd',
            'b_spline_transform', 'affine_matrix'
        ]

        return output_list_dti

    def build_input_node(self):
        """Build and connect an input node to the pipelines.
        """

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from clinica.utils.inputs import clinica_file_reader
        from clinica.utils.exceptions import ClinicaCAPSError, ClinicaException
        import clinica.utils.input_files as input_files
        from clinica.utils.stream import cprint

        all_errors = []

        # b0 Mask
        try:
            b0_mask = clinica_file_reader(self.subjects,
                                          self.sessions,
                                          self.caps_directory,
                                          input_files.DWI_PREPROC_BRAINMASK)
        except ClinicaException as e:
            all_errors.append(e)

        # DWI preprocessing NIfTI
        try:
            dwi_caps = clinica_file_reader(self.subjects,
                                           self.sessions,
                                           self.caps_directory,
                                           input_files.DWI_PREPROC_NII)
        except ClinicaException as e:
            all_errors.append(e)

        # bval files
        try:
            bval_files = clinica_file_reader(self.subjects,
                                             self.sessions,
                                             self.caps_directory,
                                             input_files.DWI_PREPROC_BVAL)
        except ClinicaException as e:
            all_errors.append(e)

        # bvec files
        try:
            bvec_files = clinica_file_reader(self.subjects,
                                             self.sessions,
                                             self.caps_directory,
                                             input_files.DWI_PREPROC_BVEC)
        except ClinicaException as e:
            all_errors.append(e)

        if len(all_errors) > 0:
            error_message = 'Clinica faced error(s) while trying to read files in your CAPS directory.\n'
            for msg in all_errors:
                error_message += str(msg)
            raise ClinicaCAPSError(error_message)

        read_input_node = npe.Node(name="LoadingCLIArguments",
                                   interface=nutil.IdentityInterface(
                                       fields=self.get_input_fields(),
                                       mandatory_inputs=True),
                                   iterables=[('preproc_dwi', dwi_caps),
                                              ('preproc_bvec', bvec_files),
                                              ('preproc_bval', bval_files),
                                              ('b0_mask', b0_mask)],
                                   synchronize=True)

        self.connect([
            (read_input_node, self.input_node, [('b0_mask', 'b0_mask')]),
            (read_input_node, self.input_node, [('preproc_dwi', 'preproc_dwi')]),
            (read_input_node, self.input_node, [('preproc_bval', 'preproc_bval')]),
            (read_input_node, self.input_node, [('preproc_bvec', 'preproc_bvec')])
        ])
        cprint('The pipeline will last approximately 20 minutes per image.')

    def build_output_node(self):
        """Build and connect an output node to the pipelines.
        """
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio
        from clinica.utils.nipype import fix_join

        import clinica.pipelines.dwi_dti.dwi_dti_utils as utils

        # Find container path from filename
        container_path = npe.Node(nutil.Function(
            input_names=['dwi_filename'],
            output_names=['container'],
            function=utils.dwi_container_from_filename),
            name='container_path')

        rename_into_caps = npe.Node(nutil.Function(
            input_names=['in_caps_dwi', 'in_norm_fa', 'in_norm_md',
                         'in_norm_ad', 'in_norm_rd', 'in_b_spline_transform',
                         'in_affine_matrix'],
            output_names=['out_caps_fa', 'out_caps_md', 'out_caps_ad',
                          'out_caps_rd', 'out_caps_b_spline_transform',
                          'out_caps_affine_matrix'],
            function=utils.rename_into_caps),
            name='rename_into_caps')

        # Writing results into CAPS
        write_results = npe.Node(name='write_results',
                                 interface=nio.DataSink())
        write_results.inputs.base_directory = self.caps_directory
        write_results.inputs.parameterization = False

        self.connect([
           (self.input_node, container_path, [('preproc_dwi', 'dwi_filename')]),  # noqa

           (container_path,   write_results, [(('container', fix_join, 'dwi', 'dti_based_processing'), 'container')]),  # noqa
           (self.output_node, write_results, [('dti',   'native_space.@dti')]),  # noqa
           (self.output_node, write_results, [('fa',    'native_space.@fa'),  # noqa
                                              ('md',    'native_space.@md'),  # noqa
                                              ('ad',    'native_space.@ad'),  # noqa
                                              ('rd',    'native_space.@rd'),  # noqa
                                              ('decfa', 'native_space.@decfa')]),  # noqa

           (self.input_node,  rename_into_caps, [('preproc_dwi',        'in_caps_dwi')]),  # noqa
           (self.output_node, rename_into_caps, [('registered_fa',      'in_norm_fa'),  # noqa
                                                 ('registered_md',      'in_norm_md'),  # noqa
                                                 ('registered_ad',      'in_norm_ad'),  # noqa
                                                 ('registered_rd',      'in_norm_rd'),  # noqa
                                                 ('affine_matrix',      'in_affine_matrix'),  # noqa
                                                 ('b_spline_transform', 'in_b_spline_transform')]),  # noqa

           (rename_into_caps, write_results, [('out_caps_fa',                 'normalized_space.@registered_fa'),  # noqa
                                              ('out_caps_md',                 'normalized_space.@registered_md'),  # noqa
                                              ('out_caps_ad',                 'normalized_space.@registered_ad'),  # noqa
                                              ('out_caps_rd',                 'normalized_space.@registered_rd'),  # noqa
                                              ('out_caps_affine_matrix',      'normalized_space.@affine_matrix'),  # noqa
                                              ('out_caps_b_spline_transform', 'normalized_space.@b_spline_transform')]),  # noqa

           (self.output_node, write_results, [('statistics_fa', 'atlas_statistics.@statistics_fa'),  # noqa
                                              ('statistics_md', 'atlas_statistics.@statistics_md'),  # noqa
                                              ('statistics_ad', 'atlas_statistics.@statistics_ad'),  # noqa
                                              ('statistics_rd', 'atlas_statistics.@statistics_rd')])   # noqa
        ])

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipelines.
        """
        import clinica.pipelines.dwi_dti.dwi_dti_workflows as workflows
        import clinica.pipelines.dwi_dti.dwi_dti_utils as utils

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.fsl as fsl
        import nipype.interfaces.mrtrix as mrtrix
        from clinica.lib.nipype.interfaces.mrtrix3.utils import TensorMetrics
        from clinica.lib.nipype.interfaces.mrtrix.preprocess import DWI2Tensor

        # Nodes creation
        # ==============
        get_bids_identifier = npe.Node(interface=nutil.Function(
            input_names=['caps_dwi_filename'],
            output_names=['bids_identifier'],
            function=utils.extract_bids_identifier_from_caps_filename),
            name='0-Get_BIDS_Identifier')

        get_caps_filenames = npe.Node(interface=nutil.Function(
            input_names=['caps_dwi_filename'],
            output_names=['bids_source', 'out_dti',
                          'out_fa', 'out_md', 'out_ad', 'out_rd'],
            function=utils.get_caps_filenames),
            name='0-CAPS_Filenames')

        convert_gradients = npe.Node(
            interface=mrtrix.FSL2MRTrix(),
            name='0-Convert_FSL_Gradient')

        dwi_to_dti = npe.Node(
            # interface=mrtrix.DWI2Tensor(),
            interface=DWI2Tensor(),
            name='1-Compute_DTI')

        dti_to_metrics = npe.Node(
            #            interface=mrtrix3.TensorMetrics(),
            interface=TensorMetrics(),
            name='2-DTI-based_Metrics')

        register_on_jhu_atlas = workflows.register_dti_maps_on_atlas(
            working_directory=self.base_dir,
            name="3-Register_DTI_Maps_On_JHU")

        scalar_analysis = npe.Node(
            interface=nutil.Function(
                input_names=['in_registered_map', 'name_map', 'prefix_file'],
                output_names=['atlas_statistics_list'],
                function=utils.statistics_on_atlases),
            name='4-Scalar_Analysis')
        scalar_analysis_fa = scalar_analysis.clone('4-Scalar_Analysis_FA')
        scalar_analysis_fa.inputs.name_map = 'FA'
        scalar_analysis_md = scalar_analysis.clone('4-Scalar_Analysis_MD')
        scalar_analysis_md.inputs.name_map = 'MD'
        scalar_analysis_ad = scalar_analysis.clone('4-Scalar_Analysis_AD')
        scalar_analysis_ad.inputs.name_map = 'AD'
        scalar_analysis_rd = scalar_analysis.clone('4-Scalar_Analysis_RD')
        scalar_analysis_rd.inputs.name_map = 'RD'

        thres_map = npe.Node(fsl.Threshold(thresh=0.0),
                             iterfield=['in_file'],
                             name='5-Remove_Negative')
        thres_fa = thres_map.clone('5-Remove_Negative_FA')
        thres_md = thres_map.clone('5-Remove_Negative_MD')
        thres_ad = thres_map.clone('5-Remove_Negative_AD')
        thres_rd = thres_map.clone('5-Remove_Negative_RD')

        print_begin_message = npe.Node(
            interface=nutil.Function(
                input_names=['in_bids_or_caps_file'],
                function=utils.print_begin_pipeline),
            name='Write-Begin_Message')

        print_end_message = npe.Node(
            interface=nutil.Function(
                input_names=['in_bids_or_caps_file', 'final_file_1', 'final_file_2'],
                function=utils.print_end_pipeline),
            name='Write-End_Message')

        # Connection
        # ==========
        self.connect([
            (self.input_node, get_caps_filenames, [('preproc_dwi', 'caps_dwi_filename')]),  # noqa
            # Print begin message
            (self.input_node, print_begin_message, [('preproc_dwi', 'in_bids_or_caps_file')]),  # noqa
            # Get BIDS/CAPS identifier from filename
            (self.input_node, get_bids_identifier, [('preproc_dwi', 'caps_dwi_filename')]),  # noqa
            # Convert FSL gradient files (bval/bvec) to MRtrix format
            (self.input_node, convert_gradients, [('preproc_bval', 'bval_file'),  # noqa
                                                  ('preproc_bvec', 'bvec_file')]),  # noqa
            # Computation of the DTI model
            (self.input_node,    dwi_to_dti, [('b0_mask',       'in_mask'),  # noqa
                                              ('preproc_dwi',   'in_file')]),  # noqa
            (convert_gradients,  dwi_to_dti, [('encoding_file', 'encoding_file')]),  # noqa
            (get_caps_filenames, dwi_to_dti, [('out_dti',       'out_filename')]),  # noqa
            # Computation of the different metrics from the DTI
            (get_caps_filenames, dti_to_metrics, [('out_fa',  'out_fa')]),  # noqa
            (get_caps_filenames, dti_to_metrics, [('out_md',  'out_adc')]),  # noqa
            (get_caps_filenames, dti_to_metrics, [('out_ad',  'out_ad')]),  # noqa
            (get_caps_filenames, dti_to_metrics, [('out_rd',  'out_rd')]),  # noqa
            (self.input_node,    dti_to_metrics, [('b0_mask', 'in_mask')]),  # noqa
            (dwi_to_dti,         dti_to_metrics, [('tensor',  'in_file')]),  # noqa
            # Register DTI maps on JHU atlas
            (dti_to_metrics,    register_on_jhu_atlas, [('out_fa',  'inputnode.in_fa'),  # noqa
                                                        ('out_adc', 'inputnode.in_md'),  # noqa
                                                        ('out_ad',  'inputnode.in_ad'),  # noqa
                                                        ('out_rd',  'inputnode.in_rd')]),  # noqa
            # Generate regional TSV files
            (get_bids_identifier,   scalar_analysis_fa, [('bids_identifier',        'prefix_file')]),  # noqa
            (register_on_jhu_atlas, scalar_analysis_fa, [('outputnode.out_norm_fa', 'in_registered_map')]),  # noqa
            (get_bids_identifier,   scalar_analysis_md, [('bids_identifier',        'prefix_file')]),  # noqa
            (register_on_jhu_atlas, scalar_analysis_md, [('outputnode.out_norm_md', 'in_registered_map')]),  # noqa
            (get_bids_identifier,   scalar_analysis_ad, [('bids_identifier',        'prefix_file')]),  # noqa
            (register_on_jhu_atlas, scalar_analysis_ad, [('outputnode.out_norm_ad', 'in_registered_map')]),  # noqa
            (get_bids_identifier,   scalar_analysis_rd, [('bids_identifier',        'prefix_file')]),  # noqa
            (register_on_jhu_atlas, scalar_analysis_rd, [('outputnode.out_norm_rd', 'in_registered_map')]),  # noqa
            # Remove negative values from the DTI maps:
            (get_caps_filenames, thres_fa, [('out_fa',  'out_file')]),  # noqa
            (dti_to_metrics,     thres_fa, [('out_fa',  'in_file')]),  # noqa

            (get_caps_filenames, thres_md, [('out_md',  'out_file')]),  # noqa
            (dti_to_metrics,     thres_md, [('out_adc', 'in_file')]),  # noqa

            (get_caps_filenames, thres_ad, [('out_ad',  'out_file')]),  # noqa
            (dti_to_metrics,     thres_ad, [('out_ad',  'in_file')]),  # noqa

            (get_caps_filenames, thres_rd, [('out_rd',  'out_file')]),  # noqa
            (dti_to_metrics,     thres_rd, [('out_rd',  'in_file')]),  # noqa
            # Outputnode
            (dwi_to_dti,            self.output_node, [('tensor',   'dti')]),  # noqa
            (thres_fa,              self.output_node, [('out_file', 'fa')]),  # noqa
            (thres_md,              self.output_node, [('out_file', 'md')]),  # noqa
            (thres_ad,              self.output_node, [('out_file', 'ad')]),  # noqa
            (thres_rd,              self.output_node, [('out_file', 'rd')]),  # noqa
            (dti_to_metrics,        self.output_node, [('out_evec', 'decfa')]),  # noqa
            (register_on_jhu_atlas, self.output_node, [('outputnode.out_norm_fa',            'registered_fa'),  # noqa
                                                       ('outputnode.out_norm_md',            'registered_md'),  # noqa
                                                       ('outputnode.out_norm_ad',            'registered_ad'),  # noqa
                                                       ('outputnode.out_norm_rd',            'registered_rd'),  # noqa
                                                       ('outputnode.out_affine_matrix',      'affine_matrix'),  # noqa
                                                       ('outputnode.out_b_spline_transform', 'b_spline_transform')]),  # noqa
            (scalar_analysis_fa,    self.output_node, [('atlas_statistics_list', 'statistics_fa')]),  # noqa
            (scalar_analysis_md,    self.output_node, [('atlas_statistics_list', 'statistics_md')]),  # noqa
            (scalar_analysis_ad,    self.output_node, [('atlas_statistics_list', 'statistics_ad')]),  # noqa
            (scalar_analysis_rd,    self.output_node, [('atlas_statistics_list', 'statistics_rd')]),   # noqa
            # Print end message
            (self.input_node,    print_end_message, [('preproc_dwi',           'in_bids_or_caps_file')]),  # noqa
            (thres_rd,           print_end_message, [('out_file',              'final_file_1')]),  # noqa
            (scalar_analysis_rd, print_end_message, [('atlas_statistics_list', 'final_file_2')]),  # noqa
        ])
