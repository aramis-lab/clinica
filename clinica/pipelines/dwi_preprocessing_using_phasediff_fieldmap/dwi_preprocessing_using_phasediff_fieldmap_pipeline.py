# coding: utf8

import clinica.pipelines.engine as cpe

# Use hash instead of parameters for iterables folder names
# Otherwise path will be too long and generate OSError
from nipype import config
cfg = dict(execution={'parameterize_dirs': False})
config.update_config(cfg)


class DwiPreprocessingUsingPhaseDiffFieldmap(cpe.Pipeline):
    """DWI Preprocessing using phase difference fieldmap.

    Args:
        input_dir(str): Input directory in a BIDS hierarchy.
        output_dir(str): Output directory in a CAPS hierarchy.
        subjects_sessions_list(str): The Subjects-Sessions list file (in .tsv
            format).

    Returns:
        A clinica pipeline object containing the DwiPreprocessingUsingPhaseDiffFieldmap pipeline.

    Raises:

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
                      'effective_echo_spacing', 'phase_encoding_direction',
                      'fmap_magnitude', 'fmap_phasediff',
                      'delta_echo_time']

        return input_list

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """
        output_list = ['preproc_dwi', 'preproc_bvec', 'preproc_bval',
                       'b0_mask']

        return output_list

    def build_input_node(self):
        """Build and connect an input node to the pipeline.
        """
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio
        from colorama import Fore
        from clinica.utils.stream import cprint
        from clinica.utils.io import check_input_bids_file, extract_metadata_from_json
        from clinica.utils.dwi import check_dwi_volume
        from clinica.utils.epi import bids_dir_to_fsl_dir

        list_dwi_files = []
        list_bvec_files = []
        list_bval_files = []
        list_effective_echo_spacings = []
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
            [eff_echo_spacing, enc_direction] = extract_metadata_from_json(dwi_json[0], ['EffectiveEchoSpacing', 'PhaseEncodingDirection'])
            list_effective_echo_spacings.append(eff_echo_spacing)
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

            cprint('From JSON files: EffectiveEchoSpacing = %s, PhaseEncodingDirection = %s, '
                   'EchoTime1 = %s, EchoTime2 = %s (DeltaEchoTime = %s)' %
                   (eff_echo_spacing, enc_direction, echo_time_1, echo_time_2, abs(echo_time_2 - echo_time_1)))

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
                                 ('effective_echo_spacing', list_effective_echo_spacings),
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
            (read_node, self.input_node, [('effective_echo_spacing', 'effective_echo_spacing')]),
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
        import clinica.pipelines.dwi_preprocessing_using_phasediff_fieldmap.dwi_preprocessing_using_phasediff_fieldmap_utils as utils

        # Find container path from DWI filename
        # =====================================
        container_path = npe.Node(nutil.Function(
            input_names=['dwi_filename'],
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
            (self.input_node, container_path,    [('dwi', 'dwi_filename')]),  # noqa
            (self.input_node,  rename_into_caps, [('dwi',          'in_bids_dwi')]),  # noqa
            (self.output_node, rename_into_caps, [('preproc_dwi',  'fname_dwi'),  # noqa
                                                  ('preproc_bval', 'fname_bval'),  # noqa
                                                  ('preproc_bvec', 'fname_bvec'),  # noqa
                                                  ('b0_mask',      'fname_brainmask')]),  # noqa
            (container_path, write_results,      [(('container', fix_join, 'dwi'), 'container')]),  # noqa
            (rename_into_caps, write_results,    [('out_caps_dwi',       'preprocessing.@preproc_dwi'),  # noqa
                                                  ('out_caps_bval',      'preprocessing.@preproc_bval'),  # noqa
                                                  ('out_caps_bvec',      'preprocessing.@preproc_bvec'),  # noqa
                                                  ('out_caps_brainmask', 'preprocessing.@b0_mask')])  # noqa
        ])

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline.
        """
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.fsl as fsl
        from nipype.workflows.dmri.fsl.utils import apply_all_corrections

        from clinica.utils.dwi import prepare_reference_b0
        from clinica.workflows.dwi_preprocessing import ecc_pipeline
        from clinica.workflows.dwi_preprocessing import hmc_pipeline
        from clinica.workflows.dwi_preprocessing import remove_bias

        import clinica.pipelines.dwi_preprocessing_using_phasediff_fieldmap.dwi_preprocessing_using_phasediff_fieldmap_workflows as workflows

        # Nodes creation
        # ==============
        # Prepare b0 image for further corrections
        prepare_b0 = npe.Node(name="PrepareB0", interface=nutil.Function(
            input_names=['in_dwi', 'in_bval', 'in_bvec', 'low_bval', 'working_directory'],
            output_names=['out_reference_b0', 'out_b0_dwi_merge',
                          'out_updated_bval', 'out_updated_bvec'],
            function=prepare_reference_b0))
        prepare_b0.inputs.low_bval = self._low_bval
        prepare_b0.inputs.working_directory = self.base_dir
        # Mask b0 for computations purposes
        mask_b0_pre = npe.Node(fsl.BET(frac=0.3, mask=True, robust=True),
                               name='PreMaskB0')
        # Head-motion correction
        hmc = hmc_pipeline(name='HeadMotionCorrection')
        # Eddy-currents correction
        ecc = ecc_pipeline(name='EddyCurrentCorrection')
        # Susceptibility distortion correction using fmap
        sdc = workflows.susceptibility_distortion_correction_using_phasediff_fmap(
            register_fmap_on_b0=True,
            name='SusceptibilityDistortionCorrection'
        )
        # Apply all corrections
        unwarp = apply_all_corrections(name='ApplyAllCorrections')
        # Remove bias correction
        bias = remove_bias(name='RemoveBias')
        #
        # Connection
        # ==========
        self.connect([
            # Preliminary step (possible computation of a mean b0)
            (self.input_node, prepare_b0, [('dwi',  'in_dwi'),  # noqa
                                           ('bval', 'in_bval'),  # noqa
                                           ('bvec', 'in_bvec')]),  # noqa
            # Mask b0 before corrections
            (prepare_b0, mask_b0_pre, [('out_reference_b0', 'in_file')]),  # noqa
            # Head-motion correction
            (prepare_b0,  hmc, [('out_b0_dwi_merge', 'inputnode.in_file'),  # noqa
                               ('out_updated_bval',  'inputnode.in_bval'),  # noqa
                               ('out_updated_bvec',  'inputnode.in_bvec')]),  # noqa
            (mask_b0_pre, hmc, [('mask_file',        'inputnode.in_mask')]),  # noqa
            # Eddy-current correction
            (hmc,         ecc, [('outputnode.out_xfms', 'inputnode.in_xfms')]),  # noqa
            (prepare_b0,  ecc, [('out_b0_dwi_merge',    'inputnode.in_file')]),  # noqa
            (prepare_b0,  ecc, [('out_updated_bval',    'inputnode.in_bval')]),  # noqa
            (mask_b0_pre, ecc, [('mask_file',           'inputnode.in_mask')]),  # noqa
            # Magnetic susceptibility correction
            (ecc,             sdc, [('outputnode.out_file',      'inputnode.in_dwi')]),  # noqa
            (mask_b0_pre,     sdc, [('mask_file',                'inputnode.in_mask')]),  # noqa
            (self.input_node, sdc, [('fmap_phasediff',           'inputnode.in_fmap_phasediff')]),  # noqa
            (self.input_node, sdc, [('fmap_magnitude',           'inputnode.in_fmap_magnitude')]),  # noqa
            (self.input_node, sdc, [('delta_echo_time',          'inputnode.delta_echo_time')]),  # noqa
            (self.input_node, sdc, [('effective_echo_spacing',   'inputnode.effective_echo_spacing')]),  # noqa
            (self.input_node, sdc, [('phase_encoding_direction', 'inputnode.fsl_phase_encoding_direction')]),  # noqa
            # Apply all corrections
            (prepare_b0, unwarp, [('out_b0_dwi_merge',    'inputnode.in_dwi')]),  # noqa
            (hmc,        unwarp, [('outputnode.out_xfms', 'inputnode.in_hmc')]),  # noqa
            (ecc,        unwarp, [('outputnode.out_xfms', 'inputnode.in_ecc')]),  # noqa
            (sdc,        unwarp, [('outputnode.out_warp', 'inputnode.in_sdc')]),  # noqa
            # Bias correction
            (unwarp, bias, [('outputnode.out_file', 'inputnode.in_file')]),
            # Outputnode
            (bias,       self.output_node, [('outputnode.out_file', 'preproc_dwi')]),  # noqa
            (hmc,        self.output_node, [('outputnode.out_bvec', 'preproc_bvec')]),  # noqa
            (prepare_b0, self.output_node, [('out_updated_bval',    'preproc_bval')]),  # noqa
            (bias,       self.output_node, [('outputnode.b0_mask',  'b0_mask')])   # noqa
        ])
