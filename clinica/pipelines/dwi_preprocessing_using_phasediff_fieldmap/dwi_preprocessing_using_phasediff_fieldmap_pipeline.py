# coding: utf8

import clinica.pipelines.engine as cpe
from nipype import config

__author__ = ["Alexandre Routier"]
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Nipype", "Junhao Wen"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__status__ = "Development"


# Use hash instead of parameters for iterables folder names
# Otherwise path will be too long and generate OSError
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
        A clinica pipeline object containing the DWIPreprocessingUsingPhaseDiffFieldmap pipeline.

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
        input_list = ['dwi', 'bvec', 'bval', 'delta_echo_time',
                      'effective_echo_spacing', 'phase_encoding_direction',
                      'fmap_magnitude', 'fmap_phasediff']

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
        from clinica.utils.dwi import check_dwi_volume
        from clinica.utils.inputs import clinica_file_reader
        import clinica.pipelines.dwi_preprocessing_using_phasediff_fieldmap.dwi_preprocessing_using_phasediff_fieldmap_utils as utils
        from clinica.utils.exceptions import ClinicaBIDSError
        from clinica.utils.stream import cprint

        all_errors = []

        # DWI preprocessing NIfTI
        dwi_bids, err_msg = clinica_file_reader(self.subjects,
                                                self.sessions,
                                                self.bids_directory,
                                                {'pattern': 'dwi/*dwi.nii*',
                                                 'description': 'DWI NifTI'})
        if err_msg:
            all_errors.append(err_msg)

        # DWI json
        dwi_json, err_msg = clinica_file_reader(self.subjects,
                                                self.sessions,
                                                self.bids_directory,
                                                {'pattern': 'dwi/*_dwi.json',
                                                 'description': 'DWI json file'})

        # Create list_eff_echo_spacings and list_enc_directions
        list_eff_echo_spacings = []
        list_enc_directions = []
        for json in dwi_json:
            [eff_echo_spacing, enc_direction] = utils.parameters_from_dwi_metadata(json)
            list_eff_echo_spacings.append(eff_echo_spacing)
            list_enc_directions.append(enc_direction)

        # bval files
        bval_files, err_msg = clinica_file_reader(self.subjects,
                                                  self.sessions,
                                                  self.bids_directory,
                                                  {'pattern': 'dwi/*_dwi.bval',
                                                   'description': 'bval files'})
        if err_msg:
            all_errors.append(err_msg)

        # bvec files
        bvec_files, err_msg = clinica_file_reader(self.subjects,
                                                  self.sessions,
                                                  self.bids_directory,
                                                  {'pattern': 'dwi/*_dwi.bvec',
                                                   'description': 'bvec files'})
        if err_msg:
            all_errors.append(err_msg)

        for (dwi, bvec, bval) in zip(dwi_bids, bvec_files, bval_files):
            check_dwi_volume(in_dwi=dwi, in_bvec=bvec, in_bval=bval)

        # Phasediff json
        fmap_phasediff_json, err_msg = clinica_file_reader(self.subjects,
                                                           self.sessions,
                                                           self.bids_directory,
                                                           {'pattern': 'fmap/*_phasediff.json',
                                                            'description': 'phasediff json file'})
        # Then deduce delta echo times
        list_delta_echo_times = [utils.delta_echo_time_from_bids_fmap(json_phasediff)
                                 for json_phasediff in fmap_phasediff_json]

        if err_msg:
            all_errors.append(err_msg)

        # Phasediff nifti
        phasediff_nifti, err_msg = clinica_file_reader(self.subjects,
                                                       self.sessions,
                                                       self.bids_directory,
                                                       {'pattern': 'fmap/*_phasediff.nii*',
                                                        'description': 'phasediff json file'})

        if err_msg:
            all_errors.append(err_msg)

        # Magnitude1
        magnitude1, err_msg = clinica_file_reader(self.subjects,
                                                  self.sessions,
                                                  self.bids_directory,
                                                  {'pattern': 'fmap/*_magnitude1.nii*',
                                                   'description': 'magnitude file'})

        if len(all_errors) > 0:
            error_message = 'Clinica faced error(s) while trying to read files in your CAPS directory.\n'
            for msg in all_errors:
                error_message += msg
            raise ClinicaBIDSError(error_message)

        cprint("List JSON parameters for DWI Preprocessing:")
        cprint("- PhaseEncodingDirections")
        cprint(list_enc_directions)
        cprint("- EffectiveEchoSpacing")
        cprint(list_eff_echo_spacings)
        cprint("- DeltaEchoTime")
        cprint(list_delta_echo_times)

        read_input_node = npe.Node(name="LoadingCLIArguments",
                                   interface=nutil.IdentityInterface(
                                       fields=self.get_input_fields(),
                                       mandatory_inputs=True),
                                   iterables=[('dwi', dwi_bids),
                                              ('bvec', bvec_files),
                                              ('bval', bval_files),
                                              ('delta_echo_time', list_delta_echo_times),
                                              ('effective_echo_spacing', list_eff_echo_spacings),
                                              ('phase_encoding_direction', list_enc_directions),
                                              ('fmap_magnitude', magnitude1),
                                              ('fmap_phasediff', phasediff_nifti)],
                                   synchronize=True)

        self.connect([
            (read_input_node, self.input_node, [('fmap_magnitude', 'fmap_magnitude')]),
            (read_input_node, self.input_node, [('fmap_phasediff', 'fmap_phasediff')]),
            (read_input_node, self.input_node, [('dwi', 'dwi')]),
            (read_input_node, self.input_node, [('bval', 'bval')]),
            (read_input_node, self.input_node, [('bvec', 'bvec')]),
            (read_input_node, self.input_node, [('delta_echo_time', 'delta_echo_time')]),
            (read_input_node, self.input_node, [('effective_echo_spacing', 'effective_echo_spacing')]),
            (read_input_node, self.input_node, [('phase_encoding_direction', 'phase_encoding_direction')])
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
            input_names=['in_dwi', 'in_bval', 'in_bvec', 'low_bval'],
            output_names=['out_reference_b0', 'out_b0_dwi_merge',
                          'out_updated_bval', 'out_updated_bvec'],
            function=prepare_reference_b0))
        prepare_b0.inputs.low_bval = self._low_bval
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
            (self.input_node, sdc, [('phase_encoding_direction', 'inputnode.phase_encoding_direction')]),  # noqa
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
