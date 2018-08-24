# coding: utf8

import clinica.pipelines.engine as cpe

__author__ = ["Alexandre Routier", "Junhao Wen", "Thomas Jacquemont"]
__copyright__ = "Copyright 2016-2018 The Aramis Lab Team"
__credits__ = ["Nipype"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__status__ = "Development"


class DWIPreprocessingUsingT1(cpe.Pipeline):
    """DWI Preprocessing using T1 image for susceptibility distortion step.

    Warnings:
        - Do not use this pipelines if you have fieldmap data in your dataset.

    Args:
        input_dir(str): Input directory in a BIDS hierarchy.
        output_dir(str): Output directory in a CAPS hierarchy.
        subjects_sessions_list(str): The Subjects-Sessions list file (in .tsv
            format).

    Returns:
        A clinica pipelines object containing the DWI Preprocessing using T1
        pipelines.

    Raises:


    Example:
        >>> pipeline = DWIPreprocessingUsingT1('~/MYDATASET_BIDS',
        >>>                                    '~/MYDATASET_CAPS')
        >>> pipeline.parameters = {
        >>>     'low_bval': 10
        >>> }
        >>> pipeline.base_dir = '/tmp/'
        >>> pipeline.run()
    """
    def __init__(self, bids_directory=None, caps_directory=None, tsv_file=None,
                 name=None, low_bval=5):
        """

        Args:
            bids_directory(str): Input directory in a BIDS hierarchy.
            caps_directory(str): Output directory in a CAPS hierarchy.
            tsv_file(str): TSV file containing the list of participants
                (participant_id) with their sessions (session_id).
            name(optional[str]): Name of the pipeline
            low_bval (int): Define the b0 volumes as all volume
                bval <= lowbval. (Default = 5)
        """
        import warnings

        super(DWIPreprocessingUsingT1, self).__init__(
            bids_directory=bids_directory,
            caps_directory=caps_directory,
            tsv_file=tsv_file,
            name=name)

        self._low_bval = low_bval

        if self._low_bval < 0:
            raise ValueError('The low_bval is equals to '
                             + str(self._low_bval)
                             + ': it should be zero or close to zero.')

        if self._low_bval > 100:
            warnings.warn('Warning: The low_bval parameter is é ('
                          + str(self._low_bval)
                          + '), it should be close to zero', UserWarning)

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file.
        """
        pass

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipelines.

        Returns:
            A list of (string) input fields name.
        """
        input_list = ['T1w', 'dwi', 'bvec', 'bval', 'epi_param']
        return input_list

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipelines.

        Returns:
            A list of (string) output fields name.
        """
        output_list = ['preproc_dwi', 'preproc_bvec', 'preproc_bval',
                       'b0_mask']

        return output_list

    def build_input_node(self):
        """Build and connect an input node to the pipelines.
        """
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio

        from clinica.utils.stream import cprint
        from clinica.utils.dwi import check_dwi_volume

        cprint('Found %s image(s) in BIDS dataset' % len(self.subjects))
        for i in range(len(self.subjects)):
            # cprint('------- SUBJECT %s SESSION %s -------'
            #        % (self.subjects[i], self.sessions[i]))

            # Check b-val file and compute the nb of b0 from file:
            bval_file = self.bids_layout.get(
                return_type='file',
                type='dwi',
                extensions=['bval'],
                session=self.sessions[i].replace('ses-', ''),
                subject=self.subjects[i].replace('sub-', '')
            )
            if len(bval_file) != 1:
                raise IOError('Expected to find 1 bval file for subject '
                              + self.subjects[i]
                              + ' and session '
                              + self.sessions[i]
                              + ' but found '
                              + str(len(bval_file))
                              + ' bval instead.')

            # Check b-vec file:
            bvec_file = self.bids_layout.get(
                return_type='file',
                type='dwi',
                extensions=['bvec'],
                session=self.sessions[i].replace('ses-', ''),
                subject=self.subjects[i].replace('sub-', '')
            )
            if len(bvec_file) != 1:
                raise IOError('Expected to find 1 bvec file for subject '
                              + self.subjects[i]
                              + ' and session '
                              + self.sessions[i]
                              + ' but found '
                              + str(len(bvec_file))
                              + ' bvec instead.')

            # Check DWI file:
            dwi_file = self.bids_layout.get(
                return_type='file',
                type='dwi',
                extensions=['.nii|.nii.gz'],
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

            # Check that the number of DWI, b-vecs & b-val are the same:
            check_dwi_volume(
                in_dwi=dwi_file[0], in_bvec=bvec_file[0], in_bval=bval_file[0])

            # Check T1w file:
            t1_file = self.bids_layout.get(
                return_type='file',
                type='T1w',
                extensions=['.nii|.nii.gz'],
                session=self.sessions[i].replace('ses-', ''),
                subject=self.subjects[i].replace('sub-', '')
            )
            if len(t1_file) != 1:
                raise IOError('Expected to find 1 T1w file for subject '
                              + self.subjects[i]
                              + ' and session '
                              + self.sessions[i]
                              + ' but found '
                              + str(len(t1_file))
                              + ' T1w instead.')

        # Iterables:
        iterables_node = npe.Node(name="LoadingCLIArguments",
                                  interface=nutil.IdentityInterface(
                                      fields=['subject_id', 'session_id'],
                                      mandatory_inputs=True)
                                  )
        iterables_node.iterables = [('subject_id', self.subjects),
                                    ('session_id', self.sessions)]
        iterables_node.synchronize = True

        # T1 DataGrabber
        t1_bids_reader = npe.Node(
            nio.DataGrabber(infields=['subject_id', 'session',
                                      'subject_repeat', 'session_repeat'],
                            outfields=['out_files']), name='t1_bids_reader')
        t1_bids_reader.inputs.base_directory = self.bids_directory
        t1_bids_reader.inputs.template = '%s/%s/anat/%s_%s_*T1w.nii*'
        t1_bids_reader.inputs.sort_filelist = False

        # DWI DataGrabber
        dwi_bids_reader = npe.Node(
            nio.DataGrabber(infields=['subject_id', 'session',
                                      'subject_repeat', 'session_repeat'],
                            outfields=['out_files']), name='dwi_bids_reader')
        dwi_bids_reader.inputs.base_directory = self.bids_directory
        dwi_bids_reader.inputs.template = '%s/%s/dwi/%s_%s_*dwi.nii*'
        dwi_bids_reader.inputs.sort_filelist = False

        # Bval DataGrabber
        bval_bids_reader = npe.Node(
            nio.DataGrabber(infields=['subject_id', 'session',
                                      'subject_repeat', 'session_repeat'],
                            outfields=['out_files']), name='bval_bids_reader')
        bval_bids_reader.inputs.base_directory = self.bids_directory
        bval_bids_reader.inputs.template = '%s/%s/dwi/%s_%s_*dwi.bval'
        bval_bids_reader.inputs.sort_filelist = False

        # Bvec dataGrabber
        bvec_bids_reader = npe.Node(
            nio.DataGrabber(infields=['subject_id', 'session',
                                      'subject_repeat', 'session_repeat'],
                            outfields=['out_files']), name='bvec_bids_reader')
        bvec_bids_reader.inputs.base_directory = self.bids_directory
        bvec_bids_reader.inputs.template = '%s/%s/dwi/%s_%s_*dwi.bvec'
        bvec_bids_reader.inputs.sort_filelist = False

        self.connect([
            # Iterables:
            (iterables_node,      t1_bids_reader,  [('subject_id',       'subject_id'),  # noqa
                                                    ('session_id',          'session'),  # noqa
                                                    ('subject_id',   'subject_repeat'),  # noqa
                                                    ('session_id', 'session_repeat')]),  # noqa
            (iterables_node,     dwi_bids_reader,  [('subject_id',       'subject_id'),  # noqa
                                                    ('session_id',          'session'),  # noqa
                                                    ('subject_id',   'subject_repeat'),  # noqa
                                                    ('session_id', 'session_repeat')]),  # noqa
            (iterables_node,     bval_bids_reader, [('subject_id',       'subject_id'),  # noqa
                                                    ('session_id',          'session'),  # noqa
                                                    ('subject_id',   'subject_repeat'),  # noqa
                                                    ('session_id', 'session_repeat')]),  # noqa
            (iterables_node,     bvec_bids_reader, [('subject_id',       'subject_id'),  # noqa
                                                    ('session_id',          'session'),  # noqa
                                                    ('subject_id',   'subject_repeat'),  # noqa
                                                    ('session_id', 'session_repeat')]),  # noqa
            # Inputnode:
            (t1_bids_reader,     self.input_node,  [('out_files',             'T1w')]),  # noqa
            (dwi_bids_reader,    self.input_node,  [('out_files',             'dwi')]),  # noqa
            (bval_bids_reader,   self.input_node,  [('out_files',            'bval')]),  # noqa
            (bvec_bids_reader,   self.input_node,  [('out_files',            'bvec')])   # noqa
        ])

    def build_output_node(self):
        """Build and connect an output node to the pipelines.
        """
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio
        from os.path import join
        import dwi_preprocessing_using_t1_utils as utils

        # Find container path from DWI filename
        # =====================================
        container_path = npe.Node(nutil.Function(
            input_names=['bids_dwi_filename'],
            output_names=['container'],
            function=utils.dwi_container_from_filename),
            name='container_path')

        rename_into_caps = npe.Node(nutil.Function(
            input_names=['in_bids_dwi', 'fname_dwi', 'fname_bval',
                         'fname_bvec', 'fname_brainmask'],
            output_names=['out_caps_dwi', 'out_caps_bval', 'out_caps_bvec',
                          'out_caps_brainmask'],
            function=utils.rename_into_caps),
            name='rename_into_caps')

        # Writing results into CAPS
        # =========================
        write_results = npe.Node(name='write_results',
                                 interface=nio.DataSink())
        write_results.inputs.base_directory = self.caps_directory
        write_results.inputs.parameterization = False

        self.connect([
            (self.input_node, container_path,    [('dwi',                    'bids_dwi_filename')]),  # noqa
            (self.input_node,  rename_into_caps, [('dwi',                          'in_bids_dwi')]),  # noqa
            (self.output_node, rename_into_caps, [('preproc_dwi',                      'fname_dwi'),  # noqa
                                                  ('preproc_bval',                    'fname_bval'),  # noqa
                                                  ('preproc_bvec',                    'fname_bvec'),  # noqa
                                                  ('b0_mask',                  'fname_brainmask')]),  # noqa
            (container_path, write_results,      [(('container', join, 'dwi'),       'container')]),  # noqa
            (rename_into_caps, write_results,    [('out_caps_dwi',    'preprocessing.@preproc_dwi'),  # noqa
                                                  ('out_caps_bval',  'preprocessing.@preproc_bval'),  # noqa
                                                  ('out_caps_bvec',  'preprocessing.@preproc_bvec'),  # noqa
                                                  ('out_caps_brainmask', 'preprocessing.@b0_mask')])  # noqa
        ])

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipelines.
        """

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.fsl as fsl

        from clinica.utils.dwi import prepare_reference_b0
        from clinica.workflows.dwi_preprocessing import eddy_fsl_pipeline
        from clinica.workflows.dwi_preprocessing import epi_pipeline
        from clinica.workflows.dwi_preprocessing import apply_all_corrections
        from clinica.workflows.dwi_preprocessing import remove_bias

        # Nodes creation
        # ==============

        # Prepare b0 image for further corrections
        prepare_b0 = npe.Node(name="PrepareB0", interface=nutil.Function(
            input_names=['in_dwi', 'in_bval', 'in_bvec', 'low_bval'],
            output_names=['out_reference_b0', 'out_b0_dwi_merge',
                          'out_updated_bval', 'out_updated_bvec'],
            function=prepare_reference_b0))
        prepare_b0.inputs.low_bval = self._low_bval
        # Mask b0 for computations purposes
        mask_b0_pre = npe.Node(fsl.BET(frac=0.3, mask=True, robust=True),
                               name='PreMaskB0')
        # Head-motion correction + Eddy-currents correction
        eddy_fsl = eddy_fsl_pipeline(self.parameters['epi_param'], name='HeadMtionandEddyCurrentCorrection')
        # Susceptibility distortion correction using T1w image
        sdc = epi_pipeline(name='SusceptibilityDistortionCorrection')
        # Remove bias correction
        bias = remove_bias(name='RemoveBias')
        # Apply all corrections
        aac = apply_all_corrections(name='ApplyAllCorrections')

        # Connection
        # ==========

        self.connect([
            # Preliminary step (possible computation of a mean b0):
            (self.input_node, prepare_b0, [('dwi',  'in_dwi'),  # noqa
                                           ('bval', 'in_bval'),  # noqa
                                           ('bvec', 'in_bvec')]),  # noqa
            # Mask b0 before corrections
            (prepare_b0, mask_b0_pre, [('out_reference_b0', 'in_file')]),  # noqa
            # Head-motion correction + eddy current correction
            (prepare_b0, eddy_fsl, [('out_b0_dwi_merge', 'inputnode.in_file'),
                                    ('out_updated_bval', 'inputnode.in_bval'),
                                    ('out_updated_bvec', 'inputnode.in_bvec'),
                                    ('out_reference_b0', 'inputnode.ref_b0')]),
            (mask_b0_pre, eddy_fsl, [('mask_file', 'inputnode.in_mask')]),
            # Magnetic susceptibility correction
            (eddy_fsl, sdc, [('outputnode.out_corrected', 'inputnode.DWI')]),
            (self.input_node, sdc, [('T1w', 'inputnode.T1')]),
            (eddy_fsl, sdc, [('outputnode.out_rotated_bvecs', 'inputnode.bvec')]),
            # Apply all corrections
            (prepare_b0,      aac, [('out_b0_dwi_merge',    'inputnode.in_dwi')]),  # noqa
            (eddy_fsl, aac, [('outputnode.out_parameter', 'inputnode.eddy_parameters')]),  # noqa
            (sdc, aac, [('outputnode.warp_epi', 'inputnode.in_epi')]),  # noqa
            (self.input_node, aac, [('T1w', 'inputnode.T1')]),
            (prepare_b0, aac, [('out_reference_b0', 'inputnode.in_hmc_ref')]),
            # Bias correction
            (aac, bias, [('outputnode.out_file', 'inputnode.in_file')]),
            # Outputnode:
            (bias,       self.output_node, [('outputnode.out_file', 'preproc_dwi')]),  # noqa
            (sdc,        self.output_node, [('outputnode.out_bvec', 'preproc_bvec')]),  # noqa
            (prepare_b0, self.output_node, [('out_updated_bval',    'preproc_bval')]),  # noqa
            (bias,       self.output_node, [('outputnode.b0_mask',  'b0_mask')])   # noqa
        ])
