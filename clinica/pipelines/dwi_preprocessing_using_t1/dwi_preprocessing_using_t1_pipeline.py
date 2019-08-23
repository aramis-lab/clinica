# coding: utf8

import clinica.pipelines.engine as cpe

__author__ = ["Junhao Wen", "Thomas Jacquemont", "Alexandre Routier"]
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Nipype"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__status__ = "Development"


# Use hash instead of parameters for iterables folder names
# Otherwise path will be too long and generate OSError
from nipype import config
cfg = dict(execution={'parameterize_dirs': False})
config.update_config(cfg)


class DwiPreprocessingUsingT1(cpe.Pipeline):
    """DWI Preprocessing using T1 image for susceptibility distortion step.

    Todo:
        [X] Detect c3d_affine_tool dependency
        [X] Refactor input_node (cf AM's comments on https://github.com/aramis-lab/clinica/pull/8)
        [X] Read data from JSON
              [X] - PhaseEncodingDirection
              [X] - TotalReadoutTime
        [ ] Interfaces
              [X] - antsRegistrationSyNQuick.sh
              [ ] - antsApplyTransforms? (optional, no output message when run)
              [ ] - CreateJacobianDeterminantImage? (optional, no output message when run)
              [ ] - WarpImageMultiTransform? (optional, no output message when run)
        [X] Add CLI flag for eddy_cuda8.0 / eddy_cuda9.1
        [ ] Update space_required_by_pipeline.csv info
        [ ] CI
              [/] Make FSL eddy reproducible
              [ ] Data CI
        [ ] Wiki page

    Ideas for improvement:
        [ ] Replace prepare_reference_b0 function by a first run of FSL eddy
        [ ] Replace B0-T1w registration by FA-T1w registration

    Questions:
        [ ] Why eddy --flm is set to linear instead of default quadratic?

    Warnings:
        - Do not use this pipeline if you have fieldmap data in your dataset.

    Returns:
        A clinica pipeline object containing the DwiPreprocessingUsingT1 pipeline.
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
            low_bval(optional[int]): Define the b0 volumes as all volume bval <= low_bval (default: 5)
            use_cuda_8_0(optional[bool]): Use CUDA 8.0 implementation of FSL eddy (default: False).
            use_cuda_9_1(optional[bool]): Use CUDA 9.1 implementation of FSL eddy (default: False).
            seed_fsl_eddy(optional[int]): Initialize the random number generator for FSL eddy (default: None).

        Raise:
            ClinicaBIDSError: If low_bval is not ranging from 0 to 100.
            ClinicaBIDSError: If CUDA 8.0 and CUDA 9.1 are chosen.
       """
        from colorama import Fore
        from clinica.utils.exceptions import ClinicaException

        super(DwiPreprocessingUsingT1, self).__init__(
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
        """Check dependencies that can not be listed in the `info.json` file."""
        pass

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipelines.

        Note:
            The list of inputs of the DwiPreprocessingUsingT1 pipeline is:
                * t1w (str): Path of the T1 weighted image in BIDS format
                * dwi (str): Path of the diffusion weighted image in BIDS format
                * bvec (str): Path of the bvec file in BIDS format
                * bval (str): Path of the bval file in BIDS format
                * dwi_json (str): Path to the DWI JSON file in BIDS format and containing
                    TotalReadoutTime and PhaseEncodingDirection metadata (see BIDS specifications)

        Returns:
            A list of (string) input fields name.
        """
        input_list = ['t1w', 'dwi', 'bvec', 'bval', 'dwi_json']
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
                                 ('t1w', t1w_files),
                                 ('dwi', dwi_files),
                                 ('bvec', bvec_files),
                                 ('bval', bval_files),
                                 ('dwi_json', dwi_json_files),
                             ],
                             synchronize=True,
                             interface=nutil.IdentityInterface(
                                 fields=self.get_input_fields()))
        self.connect([
            (read_node, self.input_node, [('t1w', 't1w'),
                                          ('dwi', 'dwi'),
                                          ('bvec', 'bvec'),
                                          ('bval', 'bval'),
                                          ('dwi_json', 'dwi_json')]),
        ])

    def build_output_node(self):
        """Build and connect an output node to the pipelines."""
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio
        from clinica.utils.io import fix_join
        import clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_utils as utils

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
            (self.input_node, container_path,    [('dwi',          'bids_dwi_filename')]),  # noqa
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
        """Build and connect the core nodes of the pipelines."""
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.fsl as fsl

        from clinica.utils.dwi import prepare_reference_b0
        from .dwi_preprocessing_using_t1_workflows import (eddy_fsl_pipeline, epi_pipeline, remove_bias)
        from .dwi_preprocessing_using_t1_utils import init_input_node, print_end_pipeline

        # Nodes creation
        # ==============
        # Initialize input parameters and print begin message
        init_node = npe.Node(interface=nutil.Function(
            input_names=self.get_input_fields(),
            output_names=['image_id',
                          't1w', 'dwi', 'bvec', 'bval', 'total_readout_time', 'phase_encoding_direction'],
            function=init_input_node),
            name='0-InitNode')

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
        # Head-motion correction + Eddy-currents correction
        eddy_fsl = eddy_fsl_pipeline(low_bval=self._low_bval,
                                     use_cuda_8_0=self._use_cuda_8_0,
                                     use_cuda_9_1=self._use_cuda_9_1,
                                     seed_fsl_eddy=self._seed_fsl_eddy)
        # Susceptibility distortion correction using T1w image
        sdc = epi_pipeline(name='SusceptibilityDistortionCorrection')
        # Remove bias correction
        bias = remove_bias(name='RemoveBias')

        # Print end message
        print_end_message = npe.Node(
            interface=nutil.Function(
                input_names=['image_id', 'final_file'],
                function=print_end_pipeline),
            name='WriteEndMessage')

        # Connection
        # ==========
        self.connect([
            # Initialize input parameters and print begin message
            (self.input_node, init_node, [('t1w', 't1w'),
                                          ('dwi', 'dwi'),
                                          ('bvec', 'bvec'),
                                          ('bval', 'bval'),
                                          ('dwi_json', 'dwi_json')]),
            # Preliminary step (possible computation of a mean b0):
            (init_node, prepare_b0, [('dwi',  'in_dwi'),
                                     ('bval', 'in_bval'),
                                     ('bvec', 'in_bvec')]),
            # Mask b0 before corrections
            (prepare_b0, mask_b0_pre, [('out_reference_b0', 'in_file')]),
            # Head-motion correction + eddy current correction
            (init_node, eddy_fsl, [('total_readout_time', 'inputnode.total_readout_time'),
                                   ('phase_encoding_direction', 'inputnode.phase_encoding_direction')]),
            (prepare_b0, eddy_fsl, [('out_b0_dwi_merge', 'inputnode.in_file'),
                                    ('out_updated_bval', 'inputnode.in_bval'),
                                    ('out_updated_bvec', 'inputnode.in_bvec'),
                                    ('out_reference_b0', 'inputnode.ref_b0')]),
            (mask_b0_pre, eddy_fsl, [('mask_file', 'inputnode.in_mask')]),
            # Magnetic susceptibility correction
            (init_node, sdc, [('t1w', 'inputnode.T1')]),
            (eddy_fsl,  sdc, [('outputnode.out_corrected', 'inputnode.DWI')]),
            (eddy_fsl,  sdc, [('outputnode.out_rotated_bvecs', 'inputnode.bvec')]),
            # Bias correction
            (sdc, bias, [('outputnode.DWIs_epicorrected', 'inputnode.in_file')]),
            # Print end message
            (init_node, print_end_message, [('image_id', 'image_id')]),  # noqa
            (bias,      print_end_message, [('outputnode.out_file', 'final_file')]),  # noqa
            # Output node
            (bias,       self.output_node, [('outputnode.out_file', 'preproc_dwi')]),  # noqa
            (sdc,        self.output_node, [('outputnode.out_bvec', 'preproc_bvec')]),  # noqa
            (prepare_b0, self.output_node, [('out_updated_bval',    'preproc_bval')]),  # noqa
            (bias,       self.output_node, [('outputnode.b0_mask',  'b0_mask')])   # noqa
        ])
