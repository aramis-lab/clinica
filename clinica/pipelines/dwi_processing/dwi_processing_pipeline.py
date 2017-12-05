# coding: utf8

import clinica.pipelines.engine as cpe


class DWIProcessing(cpe.Pipeline):
    """DWI Processing SHORT DESCRIPTION.

    Todos:
        - [ ] CAPS
        - [ ] CAPS Layout for data reading
        - [ ] Refactoring Statistics code
        - [ ] Tractogram-based processing with automatic optional registration
        - [ ] Option to the CLI to choose the DTI or Tracto processing

    Args:
        input_dir: A BIDS directory.
        output_dir: An empty output directory where CAPS structured data will
            be written.
        subjects_sessions_list: The Subjects-Sessions list file (in .tsv
            format).

    Returns:
        A clinica pipelines object containing the DWI Processing pipelines.

    Raises:


    Example:
        >>> from dwi_processing import DWIProcessing
        >>> pipelines = DWIProcessing('~/MYDATASET_BIDS', '~/MYDATASET_CAPS')
        >>> pipelines.parameters = {
        >>>     # ...
        >>> }
        >>> pipelines.base_dir = '/tmp/'
        >>> pipelines.run()
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
        input_list = ['preproc_dwi', 'preproc_bvec', 'preproc_bval', 'b0_mask']
        return input_list

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipelines.

        Returns:
            A list of (string) output fields name.
        """
        output_list_dti = [
            'dti', 'fa', 'md', 'ad', 'rd',
            'registered_fa', 'registered_md', 'registered_ad', 'registered_rd',
            'statistics_fa', 'statistics_md', 'statistics_ad', 'statistics_rd'
        ]

        return output_list_dti

    def build_input_node(self):
        """Build and connect an input node to the pipelines.
        """

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio

        from clinica.utils.stream import cprint
        from clinica.lib.pycaps.caps_layout import CAPSLayout

        caps_layout = CAPSLayout(self.caps_directory)

        cprint('Reading CAPS dataset for %s image(s)' % len(self.subjects))
        for i in range(len(self.subjects)):
            # cprint('------- SUBJECT %s SESSION %s -------'
            #        % (self.subjects[i], self.sessions[i]))

            # Check b-val file:
            bval_file = caps_layout.get(
                return_type='file',
                dwi_preprocessing_file='preproc',
                extensions=['bval'],
                regex_search=True,
                session=self.sessions[i].replace('ses-', ''),
                subject=self.subjects[i].replace('sub-', '')
            )
            if len(bval_file) != 1:
                raise IOError('Expected to find 1 bval file file for subject '
                              + self.subjects[i]
                              + ' and session '
                              + self.sessions[i]
                              + ' but found '
                              + str(len(bval_file))
                              + ' bval instead.')

            # Check b-vec file:
            bvec_file = caps_layout.get(
                return_type='file',
                dwi_preprocessing_file='preproc',
                extensions=['bvec'],
                regex_search=True,
                session=self.sessions[i].replace('ses-', ''),
                subject=self.subjects[i].replace('sub-', '')
            )
            if len(bvec_file) != 1:
                raise IOError('Expected to find 1 bvec file file for subject '
                              + self.subjects[i]
                              + ' and session '
                              + self.sessions[i]
                              + ' but found '
                              + str(len(bvec_file))
                              + ' bvec instead.')

            # Check DWI file:
            dwi_file = caps_layout.get(
                return_type='file',
                dwi_preprocessing_file='preproc',
                extensions=['.nii|.nii.gz'],
                regex_search=True,
                session=self.sessions[i].replace('ses-', ''),
                subject=self.subjects[i].replace('sub-', '')
            )
            if len(dwi_file) != 1:
                raise IOError('Expected to find 1 dwi file for subject '
                              + self.subjects[i]
                              + ' and session '
                              + self.sessions[i]
                              + ' but found '
                              + str(len(dwi_file))
                              + ' dwi instead.')

            # Check B0 file:
            b0_mask_file = caps_layout.get(
                return_type='file',
                dwi_preprocessing_file='brainmask',
                extensions=['.nii|.nii.gz'],
                regex_search=True,
                session=self.sessions[i].replace('ses-', ''),
                subject=self.subjects[i].replace('sub-', '')
            )
            if len(b0_mask_file) != 1:
                raise IOError('Expected to find 1 B0 brainmask for subject '
                              + self.subjects[i]
                              + ' and session '
                              + self.sessions[i]
                              + ' but found '
                              + str(len(b0_mask_file))
                              + ' dwi instead.')

        # Iterables:
        iterables_node = npe.Node(name="LoadingCLIArguments",
                                  interface=nutil.IdentityInterface(
                                      fields=['subject_id', 'session_id'],
                                      mandatory_inputs=True)
                                  )
        iterables_node.iterables = [('subject_id', self.subjects),
                                    ('session_id', self.sessions)]
        iterables_node.synchronize = True

        # B0 mask
        b0_mask_caps_reader = npe.Node(
            nio.DataGrabber(infields=['subject_id', 'session',
                                      'subject_repeat', 'session_repeat'],
                            outfields=['out_files']),
            name='b0_mask_caps_reader')
        b0_mask_caps_reader.inputs.base_directory = self.caps_directory
        b0_mask_caps_reader.inputs.template = 'subjects/%s/%s/dwi/preprocessing/%s_%s_*dwi_space-*_brainmask.nii.gz'  # noqa
        b0_mask_caps_reader.inputs.sort_filelist = False

        # DWI DataGrabber
        dwi_caps_reader = npe.Node(
            nio.DataGrabber(infields=['subject_id', 'session',
                                      'subject_repeat', 'session_repeat'],
                            outfields=['out_files']), name='dwi_caps_reader')
        dwi_caps_reader.inputs.base_directory = self.caps_directory
        dwi_caps_reader.inputs.template = 'subjects/%s/%s/dwi/preprocessing/%s_%s_*dwi_space-*_preproc.nii.gz'  # noqa
        dwi_caps_reader.inputs.sort_filelist = False

        # Bval DataGrabber
        bval_caps_reader = npe.Node(
            nio.DataGrabber(infields=['subject_id', 'session',
                                      'subject_repeat', 'session_repeat'],
                            outfields=['out_files']), name='bval_caps_reader')
        bval_caps_reader.inputs.base_directory = self.caps_directory
        bval_caps_reader.inputs.template = 'subjects/%s/%s/dwi/preprocessing/%s_%s_*dwi_space-*_preproc.bval'  # noqa
        bval_caps_reader.inputs.sort_filelist = False

        # Bvec dataGrabber
        bvec_caps_reader = npe.Node(
            nio.DataGrabber(infields=['subject_id', 'session',
                                      'subject_repeat', 'session_repeat'],
                            outfields=['out_files']), name='bvec_caps_reader')
        bvec_caps_reader.inputs.base_directory = self.caps_directory
        bvec_caps_reader.inputs.template = 'subjects/%s/%s/dwi/preprocessing/%s_%s_*dwi_space-*_preproc.bvec'  # noqa
        bvec_caps_reader.inputs.sort_filelist = False

        self.connect([
            # Iterables:
            (iterables_node,  b0_mask_caps_reader, [('subject_id',       'subject_id'),  # noqa
                                                    ('session_id',          'session'),  # noqa
                                                    ('subject_id',   'subject_repeat'),  # noqa
                                                    ('session_id', 'session_repeat')]),  # noqa
            (iterables_node,     dwi_caps_reader,  [('subject_id',       'subject_id'),  # noqa
                                                    ('session_id',          'session'),  # noqa
                                                    ('subject_id',   'subject_repeat'),  # noqa
                                                    ('session_id', 'session_repeat')]),  # noqa
            (iterables_node,     bval_caps_reader, [('subject_id',       'subject_id'),  # noqa
                                                    ('session_id',          'session'),  # noqa
                                                    ('subject_id',   'subject_repeat'),  # noqa
                                                    ('session_id', 'session_repeat')]),  # noqa
            (iterables_node,     bvec_caps_reader, [('subject_id',       'subject_id'),  # noqa
                                                    ('session_id',          'session'),  # noqa
                                                    ('subject_id',   'subject_repeat'),  # noqa
                                                    ('session_id', 'session_repeat')]),  # noqa
            # Inputnode:
            (b0_mask_caps_reader, self.input_node, [('out_files',      'b0_mask')]),  # noqa
            (dwi_caps_reader,     self.input_node, [('out_files',  'preproc_dwi')]),  # noqa
            (bval_caps_reader,    self.input_node, [('out_files', 'preproc_bval')]),  # noqa
            (bvec_caps_reader,    self.input_node, [('out_files', 'preproc_bvec')])  # noqa
        ])

    def build_output_node(self):
        """Build and connect an output node to the pipelines.
        """
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio
        from os.path import join

        import dwi_processing_utils as utils

        # Find container path from filename
        # =================================
        container_path = npe.Node(nutil.Function(
            input_names=['dwi_filename'],
            output_names=['container'],
            function=utils.dwi_container_from_filename),
            name='container_path')

        # Writing results into CAPS
        # =========================
        write_results = npe.Node(name='write_results',
                                 interface=nio.DataSink())
        write_results.inputs.base_directory = self.caps_directory
        write_results.inputs.parameterization = False
        write_results.inputs.regexp_substitutions = [
#            (r'(.*/)native_space/r(sub-.*)(\.nii(\.gz)?)$', r'\1\2_space-T1w_pet\3'),  # noqa
#            (r'(.*/)pet_pvc/pvc-rbv_r(sub-.*)(\.nii(\.gz)?)$', r'\1\2_space-T1w_pvc-rbv_pet\3'),
#            (r'(.*/)pet_mni/wr(sub-.*)(\.nii(\.gz)?)$', r'\1\2_space-Ixi549Space_pet\3'),
#            (r'(.*/)pet_pvc_mni/wpvc-rbv_r(sub-.*)(\.nii(\.gz)?)$', r'\1\2_space-Ixi549Space_pvc-rbv_pet\3'),
#            (r'(.*/)pet_suvr/suvr_wr(sub-.*)(\.nii(\.gz)?)$', r'\1\2_space-Ixi549Space_suvr-' + re.escape(self._suvr_region) + r'_pet\3'),
#            (r'(.*/)pet_pvc_suvr/suvr_wpvc-rbv_r(sub-.*)(\.nii(\.gz)?)$', r'\1\2_space-Ixi549Space_pvc-rbv_suvr-' + re.escape(self._suvr_region) + r'_pet\3'),
#            (r'(.*/)pet_suvr_masked/masked_suvr_wr(sub-.*)(\.nii(\.gz)?)$', r'\1\2_space-Ixi549Space_suvr-' + re.escape(self._suvr_region) + r'_mask-brain_pet\3'),
#            (r'(.*/)pet_pvc_suvr_masked/masked_suvr_wpvc-rbv_r(sub-.*)(\.nii(\.gz)?)$', r'\1\2_space-Ixi549Space_pvc-rbv_suvr-' + re.escape(self._suvr_region) + r'_mask-brain_pet\3'),
#
#
#            (r'(.*/)pet_suvr_masked_smoothed/(fwhm-[0-9]+mm)_masked_suvr_wr(sub-.*)(\.nii(\.gz)?)$', r'\1\3_space-Ixi549Space_suvr-' + re.escape(self._suvr_region) + r'_mask-brain_\2_pet\4'),
#            (r'(.*/)pet_pvc_suvr_masked_smoothed/(fwhm-[0-9]+mm)_masked_suvr_wpvc-rbv_r(sub-.*)(\.nii(\.gz)?)$', r'\1\3_space-Ixi549Space_pvc-rbv_suvr-' + re.escape(self._suvr_region) + r'_mask-brain_\2_pet\4'),
#
#            (r'(.*/)binary_mask/(sub-.*_T1w_).*(space-[a-zA-Z0-9]+).*(_brainmask\.nii(\.gz)?)$', r'\1\2\3\4')
        ]

        self.connect([
           (self.input_node, container_path, [('preproc_dwi', 'dwi_filename')]),  # noqa

           (container_path,   write_results, [(('container', join, 'dwi', 'dti_based_processing'), 'container')]),  # noqa
           (self.output_node, write_results, [('dti', 'native_space.@dti')]),
           (self.output_node, write_results, [('fa', 'native_space.@fa'),
                                              ('md', 'native_space.@md'),
                                              ('ad', 'native_space.@ad'),
                                              ('rd', 'native_space.@rd')]),

           (self.output_node, write_results, [('registered_fa', 'normalized_space.@registered_fa')]),  # noqa
           (self.output_node, write_results, [('registered_md', 'normalized_space.@registered_md')]),  # noqa
           (self.output_node, write_results, [('registered_ad', 'normalized_space.@registered_ad')]),  # noqa
           (self.output_node, write_results, [('registered_rd', 'normalized_space.@registered_rd')]),  # noqa

           (self.output_node, write_results, [('statistics_fa', 'atlas_statistics.@statistics_fa')]),  # noqa
           (self.output_node, write_results, [('statistics_md', 'atlas_statistics.@statistics_md')]),  # noqa
           (self.output_node, write_results, [('statistics_ad', 'atlas_statistics.@statistics_ad')]),  # noqa
           (self.output_node, write_results, [('statistics_rd', 'atlas_statistics.@statistics_rd')])   # noqa
        ])

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipelines.
        """
        import dwi_processing_workflows as workflows
        import dwi_processing_utils as utils

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from clinica.utils.atlas import JHUTracts01mm

        # Nodes creation
        # ==============
        nthreads = 1

        convert_to_mrtrix_format = npe.Node(interface=nutil.Function(
            input_names=['in_dwi_nii', 'in_bvals', 'in_bvecs', 'nthreads'],
            output_names=['out_dwi_mif'],
            function=utils.convert_nifti_to_mrtrix_format),
            name='convert_to_mrtrix_format')
        convert_to_mrtrix_format.inputs.nthreads = nthreads

        dwi_to_tensor = npe.Node(interface=nutil.Function(
            input_names=['in_dwi_mif', 'in_b0_mask', 'nthreads'],
            output_names=['out_dti'],
            function=utils.dwi_to_tensor), name='dwi_to_tensor')
        dwi_to_tensor.inputs.nthreads = nthreads

        get_bids_identifier = npe.Node(interface=nutil.Function(
            input_names=['caps_dwi_filename'],
            output_names=['bids_identifier'],
            function=utils.extract_bids_identifier_from_caps_filename),
            name='get_bids_identifier')
        tensor_to_metrics = npe.Node(interface=nutil.Function(
            input_names=['in_dti', 'in_b0_mask', 'nthreads', 'prefix_file'],
            output_names=['out_fa', 'out_md', 'out_ad', 'out_rd', 'out_ev'],
            function=utils.tensor_to_metrics), name='tensor_to_metrics')
        # Register DTI-maps on JHU atlas
        register_on_jhu_atlas = workflows.register_dti_maps_on_atlas(
            working_directory=self.base_dir, name="RegisterDTIMapsOnJHU")

        scalar_analysis_fa = npe.Node(interface=nutil.Function(
            input_names=['in_registered_map', 'name_map', 'prefix_file'],
            output_names=['atlas_statistics_list'],
            function=utils.statistics_on_atlases), name='ScalarAnalysisForFA')
        scalar_analysis_fa.inputs.name_map = 'fa'

        scalar_analysis_md = npe.Node(interface=nutil.Function(
            input_names=['in_registered_map', 'name_map', 'prefix_file'],
            output_names=['atlas_statistics_list'],
            function=utils.statistics_on_atlases), name='ScalarAnalysisForMD')
        scalar_analysis_md.inputs.name_map = 'md'

        scalar_analysis_ad = npe.Node(interface=nutil.Function(
            input_names=['in_registered_map', 'name_map', 'prefix_file'],
            output_names=['atlas_statistics_list'],
            function=utils.statistics_on_atlases), name='ScalarAnalysisForAD')
        scalar_analysis_ad.inputs.name_map = 'ad'

        scalar_analysis_rd = npe.Node(interface=nutil.Function(
            input_names=['in_registered_map', 'name_map', 'prefix_file'],
            output_names=['atlas_statistics_list'],
            function=utils.statistics_on_atlases), name='ScalarAnalysisForRD')
        scalar_analysis_rd.inputs.name_map = 'rd'

        # Connection
        # ==========
        self.connect([
            # Conversion to MRtrix format
            (self.input_node, convert_to_mrtrix_format, [('preproc_dwi',  'in_dwi_nii'),  # noqa
                                                         ('preproc_bval',   'in_bvals'),  # noqa
                                                         ('preproc_bvec', 'in_bvecs')]),  # noqa
            # Computation of the DTI
            (self.input_node,            dwi_to_tensor, [('b0_mask',     'in_b0_mask')]),  # noqa
            (convert_to_mrtrix_format,   dwi_to_tensor, [('out_dwi_mif', 'in_dwi_mif')]),  # noqa
            # Get BIDS identifier from filename
            (self.input_node,      get_bids_identifier, [('preproc_dwi', 'caps_dwi_filename')]),  # noqa
            # Computation of the different metrics from the DTI
            (get_bids_identifier,  tensor_to_metrics, [('bids_identifier', 'prefix_file')]),  # noqa
            (self.input_node,      tensor_to_metrics, [('b0_mask', 'in_b0_mask')]),  # noqa
            (dwi_to_tensor,        tensor_to_metrics, [('out_dti', 'in_dti')]),  # noqa
            # Register DTI maps on JHU atlas
            (tensor_to_metrics,    register_on_jhu_atlas, [('out_fa', 'inputnode.in_fa'),  # noqa
                                                           ('out_md', 'inputnode.in_md'),  # noqa
                                                           ('out_ad', 'inputnode.in_ad'),  # noqa
                                                           ('out_rd', 'inputnode.in_rd')]),  # noqa
            (get_bids_identifier,   scalar_analysis_fa, [('bids_identifier', 'prefix_file')]),  # noqa
            (register_on_jhu_atlas, scalar_analysis_fa, [('outputnode.out_registered_fa', 'in_registered_map')]),  # noqa
            (get_bids_identifier,   scalar_analysis_md, [('bids_identifier', 'prefix_file')]),  # noqa
            (register_on_jhu_atlas, scalar_analysis_md, [('outputnode.out_registered_md', 'in_registered_map')]),  # noqa
            (get_bids_identifier,   scalar_analysis_ad, [('bids_identifier', 'prefix_file')]),  # noqa
            (register_on_jhu_atlas, scalar_analysis_ad, [('outputnode.out_registered_ad', 'in_registered_map')]),  # noqa
            (get_bids_identifier,   scalar_analysis_rd, [('bids_identifier', 'prefix_file')]),  # noqa
            (register_on_jhu_atlas, scalar_analysis_rd, [('outputnode.out_registered_rd', 'in_registered_map')]),  # noqa
            # Outputnode
            (dwi_to_tensor,           self.output_node, [('out_dti', 'dti')]),
            (tensor_to_metrics,       self.output_node, [('out_fa',   'fa'),
                                                         ('out_md',   'md'),
                                                         ('out_ad',   'ad'),
                                                         ('out_rd',   'rd')]),
            (register_on_jhu_atlas,   self.output_node, [('outputnode.out_registered_fa', 'registered_fa')]),  # noqa
            (register_on_jhu_atlas,   self.output_node, [('outputnode.out_registered_md', 'registered_md')]),  # noqa
            (register_on_jhu_atlas,   self.output_node, [('outputnode.out_registered_ad', 'registered_ad')]),  # noqa
            (register_on_jhu_atlas,   self.output_node, [('outputnode.out_registered_rd', 'registered_rd')]),  # noqa
            (scalar_analysis_fa,      self.output_node, [('atlas_statistics_list',        'statistics_fa')]),  # noqa
            (scalar_analysis_md,      self.output_node, [('atlas_statistics_list',        'statistics_md')]),  # noqa
            (scalar_analysis_ad,      self.output_node, [('atlas_statistics_list',        'statistics_ad')]),  # noqa
            (scalar_analysis_rd,      self.output_node, [('atlas_statistics_list',        'statistics_rd')])   # noqa

        ])
