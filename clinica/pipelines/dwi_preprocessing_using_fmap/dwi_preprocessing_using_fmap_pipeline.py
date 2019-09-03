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
        [X] Refactor input_node (cf AM's comments on https://github.com/aramis-lab/clinica/pull/8)
        [ ] Detect & check FSL version
        [X] Add CLI flag for eddy_cuda8.0 / eddy_cuda9.1
        [ ] Decide bias correction (FSL FAST vs N4)
        [/] Core nodes
              [X] - Brain masking b0 (dwi2mask)
              [/] - Calibration FMap
              [/] - Registration FMap <-> b0
              [X] - Run eddy
              [/] - Bias correction
        [/] CAPS
              [X] - Check than CAPS does not change
              [X] - Add calibrated FMap
              [ ] - Add mean b=0?
        [ ] Update space_required_by_pipeline.csv info
        [ ] CI
        [ ] Wiki page

    Note:
        Some reading regarding the reproducibility of FSL eddy command:
        https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;1ccf038f.1608

    Returns:
        A clinica pipeline object containing the DwiPreprocessingUsingPhaseDiffFieldmap pipeline.
    """
    def __init__(self, bids_directory=None, caps_directory=None, tsv_file=None, name=None,
                 low_bval=5, use_cuda_8_0=False, use_cuda_9_1=False, seed_fsl_eddy=None):
        """Init the pipeline.

        Args:
            bids_directory(str): Input directory in a BIDS hierarchy.
            caps_directory(str): Output directory in a CAPS hierarchy.
            tsv_file(str): TSV file containing the list of participants (participant_id)
                with their sessions (session_id).
            name(optional[str]): Name of the pipeline
            low_bval(optional[int]): Define the b0 volumes as all the volumes with bval <= low_bval (default: 5)
            use_cuda_8_0(optional[bool]): Use the CUDA 8.0 implementation of FSL eddy (default: False).
            use_cuda_9_1(optional[bool]): Use the CUDA 9.1 implementation of FSL eddy (default: False).
            seed_fsl_eddy(optional[int]): Set the seed of the random number generator used
                when estimating hyperparameters in FSL eddy (default: None).

        Raise:
            ClinicaException: If low_bval is not ranging from 0 to 100.
            ClinicaException: If CUDA 8.0 and CUDA 9.1 are chosen.
       """
        from colorama import Fore
        from clinica.utils.exceptions import ClinicaException

        super(DwiPreprocessingUsingPhaseDiffFieldmap, self).__init__(
            bids_directory=bids_directory,
            caps_directory=caps_directory,
            tsv_file=tsv_file,
            name=name)

        self._low_bval = low_bval
        if self._low_bval < 0 or self._low_bval > 100:
            raise ClinicaException('\n%s[Error] The low_bval parameter should be zero or close to zero '
                                   '(found value: %s).%s' % (Fore.RED, self._low_bval, Fore.RESET))

        self._use_cuda_8_0 = use_cuda_8_0
        self._use_cuda_9_1 = use_cuda_9_1
        if self._use_cuda_8_0 and self._use_cuda_9_1:
            raise ClinicaException('\n%s[Error] You must choose between CUDA 8.0 or CUDA 9.1, not both.%s' %
                                   (Fore.RED, Fore.RESET))

        self._seed_fsl_eddy = seed_fsl_eddy

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file.
        """
        pass

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Note:
            The list of inputs of the DwiPreprocessingUsingPhaseDiffFieldmap pipeline is:
                * dwi (str): Path of the diffusion weighted image in BIDS format
                * bvec (str): Path of the bvec file in BIDS format
                * bval (str): Path of the bval file in BIDS format
                * dwi_json (str): Path to the DWI JSON file in BIDS format and containing
                    TotalReadoutTime and PhaseEncodingDirection metadata (see BIDS specifications)
                * fmap_magnitude (str): Path of the magnitude (1st) image in BIDS format
                * fmap_phasediff (str): Path of the phase difference image in BIDS format
                * fmap_phasediff_json (str): Path of the phase difference JSON file in BIDS format
                    and containing EchoTime1 & EchoTime2 metadata (see BIDS specifications)

        Returns:
            A list of (string) input fields name.
        """
        input_list = ['dwi', 'bvec', 'bval', 'dwi_json',
                      'fmap_magnitude', 'fmap_phasediff', 'fmap_phasediff_json']

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
        """Build and connect an input node to the pipeline."""
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        from colorama import Fore
        from clinica.utils.exceptions import ClinicaBIDSError
        from clinica.utils.stream import cprint
        from clinica.utils.io import check_input_bids_files

        # Remove 'sub-' prefix from participant IDs
        participant_labels = '|'.join(sub[4:] for sub in self.subjects)
        # Remove 'ses-' prefix from session IDs
        session_labels = '|'.join(ses[4:] for ses in self.sessions)

        error_message = ""
        # Inputs from anat/ folder
        # ========================
        # T1w file:
        t1w_files = self.bids_layout.get(type='T1w', return_type='file', extensions=['.nii|.nii.gz'],
                                         subject=participant_labels, session=session_labels)
        error_message += check_input_bids_files(t1w_files, "T1W_NII",
                                                self.bids_directory, self.subjects, self.sessions)

        # Inputs from dwi/ folder
        # =======================
        # Bval file:
        bval_files = self.bids_layout.get(type='dwi', return_type='file', extensions='bval',
                                          subject=participant_labels, session=session_labels)
        error_message += check_input_bids_files(bval_files, "DWI_BVAL",
                                                self.bids_directory, self.subjects, self.sessions)

        # Bvec file:
        bvec_files = self.bids_layout.get(type='dwi', return_type='file', extensions='bvec',
                                          subject=participant_labels, session=session_labels)
        error_message += check_input_bids_files(bvec_files, "DWI_BVEC",
                                                self.bids_directory, self.subjects, self.sessions)

        # DWI file:
        dwi_files = self.bids_layout.get(type='dwi', return_type='file', extensions=['.nii|.nii.gz'],
                                         subject=participant_labels, session=session_labels)
        error_message += check_input_bids_files(dwi_files, "DWI_NII",
                                                self.bids_directory, self.subjects, self.sessions)

        # DWI JSON file:
        dwi_json_files = self.bids_layout.get(type='dwi', return_type='file', extensions=['.json'],
                                              subject=participant_labels, session=session_labels)
        error_message += check_input_bids_files(dwi_json_files, "DWI_JSON",
                                                self.bids_directory, self.subjects, self.sessions)

        # Inputs from fmap/ folder
        # ========================
        # Magnitude1 file:
        fmap_magnitude_files = self.bids_layout.get(type='magnitude1', return_type='file', extensions=['.nii|.nii.gz'],
                                                    subject=participant_labels, session=session_labels)
        error_message += check_input_bids_files(fmap_magnitude_files, "FMAP_MAGNITUDE1_NII",
                                                self.bids_directory, self.subjects, self.sessions)

        # PhaseDiff file:
        fmap_phasediff_files = self.bids_layout.get(type='phasediff', return_type='file', extensions=['.nii|.nii.gz'],
                                                    subject=participant_labels, session=session_labels)
        error_message += check_input_bids_files(fmap_phasediff_files, "FMAP_PHASEDIFF_NII",
                                                self.bids_directory, self.subjects, self.sessions)

        # PhaseDiff JSON file:
        fmap_phasediff_json_files = self.bids_layout.get(type='phasediff', return_type='file', extensions=['.json'],
                                                              subject=participant_labels, session=session_labels)
        error_message += check_input_bids_files(fmap_phasediff_json_files, "FMAP_PHASEDIFF_JSON",
                                                self.bids_directory, self.subjects, self.sessions)

        if error_message:
            raise ClinicaBIDSError(error_message)

        if len(dwi_files) == 0:
            import sys
            cprint('%s\nEither all the images were already run by the pipeline or no image was found to run the pipeline. '
                   'The program will now exit.%s' % (Fore.BLUE, Fore.RESET))
            sys.exit(0)
        else:
            cprint('Found %s image(s) in BIDS dataset' % len(self.subjects))

        read_node = npe.Node(name="ReadingFiles",
                             iterables=[
                                 ('dwi', dwi_files),
                                 ('bvec', bvec_files),
                                 ('bval', bval_files),
                                 ('dwi_json', dwi_json_files),
                                 ('fmap_magnitude', fmap_magnitude_files),
                                 ('fmap_phasediff', fmap_phasediff_files),
                                 ('fmap_phasediff_json', fmap_phasediff_json_files),
                             ],
                             synchronize=True,
                             interface=nutil.IdentityInterface(
                                 fields=self.get_input_fields()))
        self.connect([
            (read_node, self.input_node, [('dwi', 'dwi'),
                                          ('bvec', 'bvec'),
                                          ('bval', 'bval'),
                                          ('dwi_json', 'dwi_json'),
                                          ('fmap_magnitude', 'fmap_magnitude'),
                                          ('fmap_phasediff', 'fmap_phasediff'),
                                          ('fmap_phasediff_json', 'fmap_phasediff_json')]),
        ])

    def build_output_node(self):
        """Build and connect an output node to the pipeline."""
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio
        from clinica.utils.io import fix_join
        from . import dwi_preprocessing_using_fmap_utils as utils

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
            (self.input_node, container_path,    [('dwi', 'dwi_filename')]),
            (self.input_node,  rename_into_caps, [('dwi',             'in_bids_dwi')]),
            (self.output_node, rename_into_caps, [('preproc_dwi',     'fname_dwi'),
                                                  ('preproc_bval',    'fname_bval'),
                                                  ('preproc_bvec',    'fname_bvec'),
                                                  ('b0_mask',         'fname_brainmask'),
                                                  ('calibrated_fmap', 'fname_calibrated_fmap')]),
            (container_path, write_results,      [(('container', fix_join, 'dwi'), 'container')]),
            (rename_into_caps, write_results,    [('out_caps_dwi',             'preprocessing.@preproc_dwi'),
                                                  ('out_caps_bval',            'preprocessing.@preproc_bval'),
                                                  ('out_caps_bvec',            'preprocessing.@preproc_bvec'),
                                                  ('out_caps_brainmask',       'preprocessing.@b0_mask'),
                                                  ('out_caps_calibrated_fmap', 'preprocessing.@calibrated_fmap')]),
            (self.output_node, write_results, [('dilate_b0_mask', 'preprocessing.@dilate_b0_mask')]),
        ])

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline."""
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.utility as niu
        import nipype.interfaces.fsl as fsl
        import nipype.interfaces.ants as ants
        import nipype.interfaces.mrtrix3 as mrtrix3

        from clinica.lib.nipype.interfaces.fsl.epi import Eddy

        from clinica.utils.dwi import generate_acq_file, generate_index_file, compute_average_b0
        from clinica.utils.fmap import remove_filename_extension

        from .dwi_preprocessing_using_fmap_workflows import (prepare_phasediff_fmap, ants_bias_correction)
        from .dwi_preprocessing_using_fmap_utils import (init_input_node,
                                                         get_grad_fsl,
                                                         print_end_pipeline)

        # Nodes creation
        # ==============
        # Initialize input parameters and print begin message
        init_node = npe.Node(interface=nutil.Function(
            input_names=self.get_input_fields(),
            output_names=['image_id',
                          'dwi', 'bvec', 'bval', 'total_readout_time', 'phase_encoding_direction',
                          'fmap_magnitude', 'fmap_phasediff', 'delta_echo_time'],
            function=init_input_node),
            name='0-InitNode')

        # Generate (bvec, bval) tuple for MRtrix interfaces
        get_grad_fsl = npe.Node(nutil.Function(
            input_names=['bval', 'bvec'],
            output_names=['grad_fsl'],
            function=get_grad_fsl),
            name='0-GetFslGrad')

        # Compute brain mask from uncorrected dataset
        pre_mask_b0 = npe.Node(mrtrix3.BrainMask(),
                               name='0-PreMaskB0')
        pre_mask_b0.inputs.out_file = 'brainmask.nii.gz'  # On default, .mif file is generated

        # Bias field correction of the magnitude image
        bias_mag_fmap = npe.Node(ants.N4BiasFieldCorrection(dimension=3), name='FMap-0-N4MagnitudeFmap')
        # Brain extraction of the magnitude image
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

        # Generate <image_id>_acq.txt for eddy
        gen_acq_txt = npe.Node(nutil.Function(
            input_names=['in_dwi', 'fsl_phase_encoding_direction', 'total_readout_time', 'image_id'],
            output_names=['out_acq'],
            function=generate_acq_file),
            name='0-GenerateAcqFile')

        # Generate <image_id>_index.txt for eddy
        gen_index_txt = npe.Node(nutil.Function(
            input_names=['in_bval', 'low_bval', 'image_id'],
            output_names=['out_index'],
            function=generate_index_file),
            name='0-GenerateIndexFile')
        gen_index_txt.inputs.low_bval = self._low_bval

        # Remove ".nii.gz" from fieldmap filename for eddy --field
        rm_extension = npe.Node(interface=nutil.Function(
            input_names=['in_file'],
            output_names=['file_without_extension'],
            function=remove_filename_extension),
            name="0-RemoveFNameExtension")

        dilate = npe.Node(fsl.maths.MathsCommand(
            nan2zeros=True, args='-kernel sphere 5 -dilM'),
            name='dilate')

        # Run FSL eddy
        eddy = npe.Node(name='2b-Eddy',
                        interface=Eddy())
        eddy.inputs.repol = True
        if self._use_cuda_8_0:
            eddy.inputs.use_cuda8_0 = self._use_cuda_8_0
        if self._use_cuda_9_1:
            eddy.inputs.use_cuda9_1 = self._use_cuda_9_1
        if self._seed_fsl_eddy:
            eddy.inputs.initrand = self._seed_fsl_eddy

        # Bias correction
        bias = ants_bias_correction(name='RemoveBias')

        # Compute average b0 (for brain mask extraction)
        compute_average_b0 = npe.Node(niu.Function(input_names=['in_dwi', 'in_bval'],
                                                   output_names=['out_b0_average'],
                                                   function=compute_average_b0),
                                      name='ComputeB0Average')
        compute_average_b0.inputs.low_bval = 5.0

        # Compute b0 mask on corrected DWI
        mask_b0 = npe.Node(fsl.BET(mask=True, robust=True),
                           name='mask_b0')

        # Print end message
        print_end_message = npe.Node(
            interface=nutil.Function(
                input_names=['image_id', 'final_file'],
                function=print_end_pipeline),
            name='99-WriteEndMessage')

        # Connection
        # ==========
        self.connect([
            # Initialize input parameters and print begin message
            (self.input_node, init_node, [('dwi', 'dwi'),
                                          ('bvec', 'bvec'),
                                          ('bval', 'bval'),
                                          ('dwi_json', 'dwi_json'),
                                          ('fmap_magnitude', 'fmap_magnitude'),
                                          ('fmap_phasediff', 'fmap_phasediff'),
                                          ('fmap_phasediff_json', 'fmap_phasediff_json')]),
            # Generate (bvec, bval) tuple for MRtrix interfaces
            (init_node, get_grad_fsl, [('bval', 'bval'),
                                       ('bvec', 'bvec')]),
            (get_grad_fsl, pre_mask_b0, [('grad_fsl',  'grad_fsl')]),
            (init_node,    pre_mask_b0, [('dwi',  'in_file')]),
            # Bias field correction of the magnitude image
            (init_node, bias_mag_fmap, [('fmap_magnitude', 'input_image')]),
            # Brain extraction of the magnitude image
            (bias_mag_fmap, bet_mag_fmap, [('output_image', 'in_file')]),
            # Calibration of the FMap
            (bet_mag_fmap, calibrate_fmap, [('mask_file', 'input_node.fmap_mask')]),
            (init_node, calibrate_fmap, [('fmap_phasediff', 'input_node.fmap_phasediff'),
                                         ('fmap_magnitude', 'input_node.fmap_magnitude'),
                                         ('delta_echo_time', 'input_node.delta_echo_time')]),
            # FMap <-> B0 registration

            (pre_mask_b0,    res_fmap, [('out_file', 'in_b0')]),
            (calibrate_fmap, res_fmap, [('output_node.calibrated_fmap', 'in_fmap')]),
            (res_fmap, smoothing, [('out_resampled_fmap', 'in_file')]),

            # Remove ".nii.gz" from fieldmap filename for eddy --field
            (smoothing, rm_extension, [('out_file', 'in_file')]),
            # Generate <image_id>_acq.txt for eddy
            (init_node, gen_acq_txt, [('dwi', 'in_dwi'),
                                      ('total_readout_time', 'total_readout_time'),
                                      ('phase_encoding_direction', 'fsl_phase_encoding_direction'),
                                      ('image_id', 'image_id')]),
            # Generate <image_id>_index.txt for eddy
            (init_node, gen_index_txt, [('bval', 'in_bval'),
                                        ('image_id', 'image_id')]),
            # Run FSL eddy
            (init_node,     eddy, [('dwi', 'in_file'),
                                   ('bval', 'in_bval'),
                                   ('bvec', 'in_bvec'),
                                   ('image_id', 'out_base')]),
            (gen_acq_txt,   eddy, [('out_acq', 'in_acqp')]),
            (gen_index_txt, eddy, [('out_index', 'in_index')]),
            (rm_extension,  eddy, [('file_without_extension', 'field')]),
            (pre_mask_b0, dilate, [('out_file', 'in_file')]),
            #            (dilate, eddy, [('out_file', 'in_mask')]),
            (pre_mask_b0,   eddy, [('out_file', 'in_mask')]),
            # Bias correction
            (pre_mask_b0, bias, [('out_file', 'input_node.mask')]),
            (eddy, bias, [('out_corrected', 'input_node.dwi'),
                          ('out_rotated_bvecs', 'input_node.bvec')]),
            (init_node, bias, [('bval', 'input_node.bval')]),
            # Compute average b0 (for brain mask extraction)
            (init_node, compute_average_b0, [('bval', 'in_bval')]),
            (bias, compute_average_b0, [('output_node.bias_corrected_dwi', 'in_dwi')]),
            # Compute b0 mask on corrected DWI
            (compute_average_b0, mask_b0, [('out_b0_average', 'in_file')]),
            # Print end message
            (init_node, print_end_message, [('image_id', 'image_id')]),
            (mask_b0,   print_end_message, [('mask_file', 'final_file')]),
            # Output node
            (init_node, self.output_node, [('bval', 'preproc_bval')]),
            (eddy,      self.output_node, [('out_rotated_bvecs',   'preproc_bvec')]),
            (bias,      self.output_node, [('output_node.bias_corrected_dwi', 'preproc_dwi')]),
            (mask_b0,   self.output_node, [('mask_file',  'b0_mask')]),
            (res_fmap,  self.output_node, [('out_resampled_fmap',  'calibrated_fmap')]),

            (dilate,    self.output_node, [('out_file',            'dilate_b0_mask')]),
        ])
