# coding: utf8

import clinica.pipelines.engine as cpe

# Use hash instead of parameters for iterables folder names
# Otherwise path will be too long and generate OSError
from nipype import config
cfg = dict(execution={'parameterize_dirs': False})
config.update_config(cfg)


class DwiPreprocessingUsingPhaseDiffFieldmap(cpe.Pipeline):
    """DWI Preprocessing using phase difference fieldmap.

    Todo:
        [ ] Refactor input_node (cf AM's comments on https://github.com/aramis-lab/clinica/pull/8)
        [ ] Detect & check FSL version
        [ ] Add CLI flag for eddy_cuda8.0 / eddy_cuda9.1
        [ ] Decide bias correction (FSL FAST vs N4)
        [ ] Core nodes
              [X] - Brain masking b0 (dwi2mask)
              [/] - Calibration FMap
              [ ] - Registration FMap <-> b0
              [X] - Run eddy
              [/] - Bias correction
              [ ] - Compute mean b=0?
        [ ] CAPS
              [X] - Check than CAPS does not change
              [X] - Add calibrated FMap
              [ ] - Add mean b=0?
        [ ] CI
              [ ] - Chose parameters so that FSL eddy is reproducible
                    [ ] - Set random init (--init)
                    [ ] - Other
              [ ] - Data CI
        [ ] Wiki page

    Returns:
        A clinica pipeline object containing the DwiPreprocessingUsingPhaseDiffFieldmap pipeline.
    """
    def __init__(self, bids_directory=None, caps_directory=None, tsv_file=None,
                 name=None, low_bval=5):
        """
        Init the pipeline

        Args:
            bids_directory(str): Input directory in a BIDS hierarchy.
            caps_directory(str): Output directory in a CAPS hierarchy.
            tsv_file(str): TSV file containing the list of participants
                (participant_id) with their sessions (session_id).
            name(optional[str]): Name of the pipeline
            low_bval (int): Define the b0 volumes as all volume
                bval <= low_bval. (Default = 5)
        """
        import warnings

        super(DwiPreprocessingUsingPhaseDiffFieldmap, self).__init__(
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
            warnings.warn('Warning: The low_bval parameter is huge ('
                          + str(self._low_bval)
                          + '), it should be close to zero', UserWarning)

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file.
        """
        pass

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Returns:
            A list of (string) input fields name.
        """
        input_list = ['dwi', 'bvec', 'bval',
                      'total_readout_time', 'phase_encoding_direction',
                      'fmap_magnitude', 'fmap_phasediff',
                      'delta_echo_time']

        return input_list

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """
        output_list = ['preproc_dwi', 'preproc_bvec', 'preproc_bval',
                       'b0_mask', 'calibrated_fmap', 'dilate_b0_mask']

        return output_list

    def build_input_node(self):
        """Build and connect an input node to the pipeline.
        """
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from colorama import Fore
        from clinica.utils.stream import cprint
        from clinica.utils.io import check_input_bids_file, extract_metadata_from_json
        from clinica.utils.dwi import check_dwi_volume
        from clinica.utils.epi import bids_dir_to_fsl_dir

        list_dwi_files = []
        list_bvec_files = []
        list_bval_files = []
        list_total_readout_times = []
        list_phase_encoding_directions = []
        list_fmap_magnitude_files = []
        list_fmap_phasediff_files = []
        list_delta_echo_times = []

        cprint('Reading input files...')
        for i in range(len(self.subjects)):
            cprint('\t...subject \'' + str(
                    self.subjects[i][4:]) + '\', session \'' + str(
                    self.sessions[i][4:]) + '\'')

            # Inputs from dwi/ folder
            # =======================
            # Bval file:
            bval_file = self.bids_layout.get(type='dwi', return_type='file', extensions=['bval'],
                                             subject=self.subjects[i][4:], session=self.sessions[i][4:])
            check_input_bids_file(bval_file, "DWI_BVAL",
                                  self.bids_directory, self.subjects[i], self.sessions[i])
            list_bval_files.append(bval_file[0])

            # Bvec file:
            bvec_file = self.bids_layout.get(type='dwi', return_type='file', extensions=['bvec'],
                                             subject=self.subjects[i][4:], session=self.sessions[i][4:])
            check_input_bids_file(bvec_file, "DWI_BVEC",
                                  self.bids_directory, self.subjects[i], self.sessions[i])
            list_bvec_files.append(bvec_file[0])

            # DWI file:
            dwi_file = self.bids_layout.get(type='dwi', return_type='file', extensions=['.nii|.nii.gz'],
                                            subject=self.subjects[i][4:], session=self.sessions[i][4:])
            check_input_bids_file(bvec_file, "DWI_NII",
                                  self.bids_directory, self.subjects[i], self.sessions[i])
            list_dwi_files.append(dwi_file[0])

            # Check that the number of DWI, bvec & bval are the same:
            check_dwi_volume(in_dwi=dwi_file[0], in_bvec=bvec_file[0], in_bval=bval_file[0])

            # DWI JSON file:
            dwi_json = self.bids_layout.get(type='dwi', return_type='file', extensions=['.json'],
                                            subject=self.subjects[i][4:], session=self.sessions[i][4:])
            check_input_bids_file(bvec_file, "DWI_JSON",
                                  self.bids_directory, self.subjects[i], self.sessions[i])
            [total_readout_time, enc_direction] = extract_metadata_from_json(dwi_json[0], ['TotalReadoutTime', 'PhaseEncodingDirection'])
            list_total_readout_times.append(total_readout_time)
            list_phase_encoding_directions.append(bids_dir_to_fsl_dir(enc_direction))

            # Inputs from fmap/ folder
            # ========================
            # Magnitude1 file:
            fmap_mag_file = self.bids_layout.get(type='magnitude1', return_type='file', extensions=['.nii|.nii.gz'],
                                                 subject=self.subjects[i][4:], session=self.sessions[i][4:])
            check_input_bids_file(fmap_mag_file, "FMAP_MAGNITUDE1_NII",
                                  self.bids_directory, self.subjects[i], self.sessions[i])
            list_fmap_magnitude_files.append(fmap_mag_file[0])

            # PhaseDiff file:
            fmap_phasediff_file = self.bids_layout.get(type='phasediff', return_type='file', extensions=['.nii|.nii.gz'],
                                                       subject=self.subjects[i][4:], session=self.sessions[i][4:])
            check_input_bids_file(fmap_phasediff_file, "FMAP_PHASEDIFF_NII",
                                  self.bids_directory, self.subjects[i], self.sessions[i])
            list_fmap_phasediff_files.append(fmap_phasediff_file[0])

            # PhaseDiff JSON file:
            fmap_phasediff_json = self.bids_layout.get(type='phasediff', return_type='file', extensions=['.json'],
                                                       subject=self.subjects[i][4:], session=self.sessions[i][4:])
            check_input_bids_file(fmap_phasediff_json, "FMAP_PHASEDIFF_JSON",
                                  self.bids_directory, self.subjects[i], self.sessions[i])
            [echo_time_1, echo_time_2] = extract_metadata_from_json(fmap_phasediff_json[0], ['EchoTime1', 'EchoTime2'])
            list_delta_echo_times.append(abs(echo_time_2 - echo_time_1))

            cprint('From JSON files: TotalReadoutTime = %s, PhaseEncodingDirection = %s, '
                   'EchoTime1 = %s, EchoTime2 = %s (DeltaEchoTime = %s)' %
                   (total_readout_time, enc_direction, echo_time_1, echo_time_2, abs(echo_time_2 - echo_time_1)))

        if len(list_dwi_files) == 0:
            import sys
            cprint('%s\nEither all the images were already run by the pipeline or no image was found to run the pipeline. '
                   'The program will now exit.%s' % (Fore.BLUE, Fore.RESET))
            sys.exit(0)
        else:
            cprint('Found %s image(s) in BIDS dataset' % len(self.subjects))

        read_node = npe.Node(name="ReadingFiles",
                             iterables=[
                                 ('dwi', list_dwi_files),
                                 ('bvec', list_bvec_files),
                                 ('bval', list_bval_files),
                                 ('total_readout_time', list_total_readout_times),
                                 ('phase_encoding_direction', list_phase_encoding_directions),
                                 ('fmap_magnitude', list_fmap_magnitude_files),
                                 ('fmap_phasediff', list_fmap_phasediff_files),
                                 ('delta_echo_time', list_delta_echo_times),
                             ],
                             synchronize=True,
                             interface=nutil.IdentityInterface(
                                 fields=self.get_input_fields()))
        self.connect([
            (read_node, self.input_node, [('dwi', 'dwi')]),
            (read_node, self.input_node, [('bvec', 'bvec')]),
            (read_node, self.input_node, [('bval', 'bval')]),
            (read_node, self.input_node, [('total_readout_time', 'total_readout_time')]),
            (read_node, self.input_node, [('phase_encoding_direction', 'phase_encoding_direction')]),
            (read_node, self.input_node, [('fmap_magnitude', 'fmap_magnitude')]),
            (read_node, self.input_node, [('fmap_phasediff', 'fmap_phasediff')]),
            (read_node, self.input_node, [('delta_echo_time', 'delta_echo_time')]),
        ])

    def build_output_node(self):
        """Build and connect an output node to the pipeline.
        """
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio
        from clinica.utils.io import fix_join
        import clinica.pipelines.dwi_preprocessing_using_fmap.dwi_preprocessing_using_fmap_utils as utils

        # Find container path from DWI filename
        # =====================================
        container_path = npe.Node(nutil.Function(
            input_names=['dwi_filename'],
            output_names=['container'],
            function=utils.dwi_container_from_filename),
            name='container_path')

        rename_into_caps = npe.Node(nutil.Function(
            input_names=['in_bids_dwi', 'fname_dwi', 'fname_bval',
                         'fname_bvec', 'fname_brainmask',  'fname_calibrated_fmap'],
            output_names=['out_caps_dwi', 'out_caps_bval', 'out_caps_bvec',
                          'out_caps_brainmask', 'out_caps_calibrated_fmap'],
            function=utils.rename_into_caps),
            name='rename_into_caps')

        # Writing results into CAPS
        # =========================
        write_results = npe.Node(name='write_results',
                                 interface=nio.DataSink())
        write_results.inputs.base_directory = self.caps_directory
        write_results.inputs.parameterization = False

        self.connect([
            (self.input_node, container_path,    [('dwi', 'dwi_filename')]),  # noqa
            (self.input_node,  rename_into_caps, [('dwi',             'in_bids_dwi')]),  # noqa
            (self.output_node, rename_into_caps, [('preproc_dwi',     'fname_dwi'),  # noqa
                                                  ('preproc_bval',    'fname_bval'),  # noqa
                                                  ('preproc_bvec',    'fname_bvec'),  # noqa
                                                  ('b0_mask',         'fname_brainmask'),  # noqa
                                                  ('calibrated_fmap', 'fname_calibrated_fmap')]),  # noqa
            (container_path, write_results,      [(('container', fix_join, 'dwi'), 'container')]),  # noqa
            (rename_into_caps, write_results,    [('out_caps_dwi',             'preprocessing.@preproc_dwi'),  # noqa
                                                  ('out_caps_bval',            'preprocessing.@preproc_bval'),  # noqa
                                                  ('out_caps_bvec',            'preprocessing.@preproc_bvec'),  # noqa
                                                  ('out_caps_brainmask',       'preprocessing.@b0_mask'),  # noqa
                                                  ('out_caps_calibrated_fmap', 'preprocessing.@calibrated_fmap')]),  # noqa
            (self.output_node, write_results, [('dilate_b0_mask', 'preprocessing.@dilate_b0_mask')]),
        ])

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline.
        """
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.fsl as fsl
        import nipype.interfaces.ants as ants
        import nipype.interfaces.mrtrix3 as mrtrix3

        from clinica.utils.dwi import generate_acq_file, generate_index_file
        from clinica.utils.fmap import remove_filename_extension
        from clinica.workflows.dwi_preprocessing import remove_bias

        from .dwi_preprocessing_using_fmap_workflows import prepare_phasediff_fmap
        import clinica.pipelines.dwi_preprocessing_using_fmap.dwi_preprocessing_using_fmap_utils as utils

        # Nodes creation
        # ==============
        # Get <subject_id> (e.g. sub-CLNC01_ses-M00) from input_node
        # and print begin message
        init_node = npe.Node(interface=nutil.Function(
            input_names=self.get_input_fields(),
            output_names=['subject_id'] + self.get_input_fields(),
            function=utils.init_input_node),
            name='0-InitNode')

        # Compute brain mask from uncorrected dataset
        pre_mask_b0 = npe.Node(mrtrix3.BrainMask(),
                               name='0-PreMaskB0')
        pre_mask_b0.inputs.out_file = 'brainmask.nii.gz'  # On default, .mif file is generated

        # Generate <subject_id>_acq.txt for eddy
        get_grad_fsl = npe.Node(nutil.Function(
            input_names=['bval', 'bvec'],
            output_names=['grad_fsl'],
            function=utils.get_grad_fsl),
            name='0-GetFslGrad')

        # Step 0 - BET & Bias field correction of the magnitude image
        n4 = npe.Node(ants.N4BiasFieldCorrection(dimension=3), name='FMap-0-N4MagnitudeFmap')

        bet_mag_fmap = npe.Node(fsl.BET(frac=0.4, mask=True), name='FMap-0-BetN4MagnitudeFmap')

        # Calibrate FMap
        calibrate_fmap = prepare_phasediff_fmap(name='FMap-1-CalibrateFMap')

        # Node when we register the fmap onto the b0
        fmm2b0 = npe.Node(ants.Registration(
            output_warped_image=True), name="FmapMagnitudeToB0")
        fmm2b0.inputs.transforms = ['Rigid'] * 2
        fmm2b0.inputs.transform_parameters = [(1.0,)] * 2
        fmm2b0.inputs.number_of_iterations = [[50], [20]]
        fmm2b0.inputs.dimension = 3
        fmm2b0.inputs.metric = ['Mattes', 'Mattes']
        fmm2b0.inputs.metric_weight = [1.0] * 2
        fmm2b0.inputs.radius_or_number_of_bins = [64, 64]
        fmm2b0.inputs.sampling_strategy = ['Regular', 'Random']
        fmm2b0.inputs.sampling_percentage = [None, 0.2]
        fmm2b0.inputs.convergence_threshold = [1.e-5, 1.e-8]
        fmm2b0.inputs.convergence_window_size = [20, 10]
        fmm2b0.inputs.smoothing_sigmas = [[6.0], [2.0]]
        fmm2b0.inputs.sigma_units = ['vox'] * 2
        fmm2b0.inputs.shrink_factors = [[6], [1]]  # ,[1] ]
        fmm2b0.inputs.use_estimate_learning_rate_once = [True] * 2
        fmm2b0.inputs.use_histogram_matching = [True] * 2
        fmm2b0.inputs.initial_moving_transform_com = 0
        fmm2b0.inputs.collapse_output_transforms = True
        fmm2b0.inputs.winsorize_upper_quantile = 0.995

        from clinica.utils.fmap import resample_fmap_to_b0
        res_fmap = npe.Node(nutil.Function(
            input_names=['in_fmap', 'in_b0', 'out_file'],
            output_names=['out_resampled_fmap'],
            function=resample_fmap_to_b0), name='ResampleFmap')

        smoothing = npe.Node(name='Smoothing',
                             interface=fsl.maths.IsotropicSmooth())
        smoothing.inputs.fwhm = 4.0

        # Generate <subject_id>_acq.txt for eddy
        gen_acq_txt = npe.Node(nutil.Function(
            input_names=['in_dwi', 'fsl_phase_encoding_direction', 'total_readout_time', 'subject_id'],
            output_names=['out_acq'],
            function=generate_acq_file),
            name='0-GenerateAcqFile')

        # Generate <subject_id>_index.txt for eddy
        gen_index_txt = npe.Node(nutil.Function(
            input_names=['in_bval', 'low_bval', 'subject_id'],
            output_names=['out_index'],
            function=generate_index_file),
            name='0-GenerateIndexFile')
        gen_index_txt.inputs.low_bval = self._low_bval

        # Remove ".nii.gz" from fieldmap filename for eddy
        rm_extension = npe.Node(interface=nutil.Function(
            input_names=['in_file'],
            output_names=['file_without_extension'],
            function=remove_filename_extension),
            name="0-RemoveFNameExtension")

        dilate = npe.Node(fsl.maths.MathsCommand(nan2zeros=True,
                                                 args='-kernel sphere 5 -dilM'),
                          name='dilate')

        # Run FSL eddy
        eddy = npe.Node(name='2b-Eddy',
                        interface=fsl.Eddy())
        eddy.inputs.repol = True
        eddy.inputs.use_cuda = True

        # Bias correction
        # TODO: Be careful with b0 extraction
        bias = remove_bias(name='RemoveBias')

        # Print end message
        print_end_message = npe.Node(
            interface=nutil.Function(
                input_names=['subject_id', 'final_file'],
                function=utils.print_end_pipeline),
            name='99-WriteEndMessage')

        # Connection
        # ==========
        self.connect([
            (self.input_node, init_node, [('dwi', 'dwi'),
                                          ('bvec', 'bvec'),
                                          ('bval', 'bval'),
                                          ('total_readout_time', 'total_readout_time'),
                                          ('phase_encoding_direction', 'phase_encoding_direction'),
                                          ('fmap_magnitude', 'fmap_magnitude'),
                                          ('fmap_phasediff', 'fmap_phasediff'),
                                          ('delta_echo_time', 'delta_echo_time')]),  # noqa

            (init_node, get_grad_fsl, [('bval', 'bval'),  # noqa
                                       ('bvec', 'bvec')]),  # noqa
#            (self.input_node, pre_mask_b0, [('dwi',  'in_file'),  # noqa
#                                            ('bval', 'in_bval'),  # noqa
#                                            ('bvec', 'in_bvec')]),  # noqa
            (get_grad_fsl,    pre_mask_b0, [('grad_fsl',  'grad_fsl')]),  # noqa
            (init_node, pre_mask_b0, [('dwi',  'in_file')]),  # noqa
            # Calibration of the FMap
            (init_node, n4, [('fmap_magnitude', 'input_image')]),  # noqa
            (n4, bet_mag_fmap, [('output_image', 'in_file')]),  # noqa
            (bet_mag_fmap,    calibrate_fmap, [('mask_file', 'input_node.fmap_mask')]),  # noqa
            (init_node, calibrate_fmap, [('fmap_phasediff', 'input_node.fmap_phasediff')]),  # noqa
            (init_node, calibrate_fmap, [('fmap_magnitude', 'input_node.fmap_magnitude')]),  # noqa
            (init_node, calibrate_fmap, [('delta_echo_time', 'input_node.delta_echo_time')]),  # noqa

            (pre_mask_b0,    res_fmap, [('out_file', 'in_b0')]),  # noqa
#            (prepare_b0, res_fmap, [('out_reference_b0', 'in_b0')]),  # noqa
            (calibrate_fmap, res_fmap, [('output_node.calibrated_fmap', 'in_fmap')]),  # noqa
            (res_fmap, smoothing, [('out_resampled_fmap', 'in_file')]),  # noqa

            # Remove ".nii.gz" from fieldmap filename for eddy
            (smoothing, rm_extension, [('out_file', 'in_file')]),  # noqa
            # Generate <subject_id>_acq.txt for eddy
            (init_node, gen_acq_txt, [('dwi', 'in_dwi'),  # noqa
                                      ('total_readout_time', 'total_readout_time'),  # noqa
                                      ('phase_encoding_direction', 'fsl_phase_encoding_direction'),  # noqa
                                      ('subject_id', 'subject_id')]),  # noqa
            # Generate <subject_id>_index.txt for eddy
            (init_node, gen_index_txt, [('bval', 'in_bval'),  # noqa
                                        ('subject_id', 'subject_id')]),  # noqa
            # Run FSL eddy
            (init_node,     eddy, [('dwi', 'in_file'),  # noqa
                                   ('bval', 'in_bval'),  # noqa
                                   ('bvec', 'in_bvec'),  # noqa
                                   ('subject_id', 'out_base')]),  # noqa
            (gen_acq_txt,   eddy, [('out_acq', 'in_acqp')]),  # noqa
            (gen_index_txt, eddy, [('out_index', 'in_index')]),  # noqa
            (rm_extension,  eddy, [('file_without_extension', 'field')]),  # noqa
            (pre_mask_b0,   eddy, [('out_file', 'in_mask')]),  # noqa

            (pre_mask_b0, dilate, [('out_file', 'in_file')]),  # noqa
#            (dilate, eddy, [('out_file', 'in_mask')]),  # noqa
            # Bias correction
            (eddy, bias, [('out_corrected', 'inputnode.in_file')]),
            # Print end message
            (init_node, print_end_message, [('subject_id', 'subject_id')]),  # noqa
            (bias,      print_end_message, [('outputnode.out_file', 'final_file')]),  # noqa
            # Output node
            (init_node, self.output_node, [('bval', 'preproc_bval')]),  # noqa
            (eddy,      self.output_node, [('out_rotated_bvecs',   'preproc_bvec')]),  # noqa
            (bias,      self.output_node, [('outputnode.out_file', 'preproc_dwi')]),  # noqa
            (bias,      self.output_node, [('outputnode.b0_mask',  'b0_mask')]),  # noqa
            # (bias,     self.output_node, [('outputnode.out_file', 'preproc_dwi')]),  # noqa
            (res_fmap,  self.output_node, [('out_resampled_fmap',  'calibrated_fmap')]),  # noqa
            (dilate,    self.output_node, [('out_file',            'dilate_b0_mask')]),  # noqa
        ])
