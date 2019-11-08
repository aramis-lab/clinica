# coding: utf8

import clinica.pipelines.engine as cpe
from nipype import config

__author__ = ["Thomas Jacquemont", "Alexandre Routier"]
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Nipype"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__status__ = "Development"


# Use hash instead of parameters for iterables folder names
# Otherwise path will be too long and generate OSError
cfg = dict(execution={'parameterize_dirs': False})
config.update_config(cfg)


class DwiPreprocessingUsingT1(cpe.Pipeline):
    """DWI Preprocessing using T1 image for susceptibility distortion step.

    Warnings:
        - Do not use this pipelines if you have fieldmap data in your dataset.

    Args:
        input_dir(str): Input directory in a BIDS hierarchy.
        output_dir(str): Output directory in a CAPS hierarchy.
        subjects_sessions_list(str): The Subjects-Sessions list file (in .tsv
            format).

    Returns:
        A clinica pipeline object containing the DWIPreprocessingUsingT1 pipeline.

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

        super(DwiPreprocessingUsingT1, self).__init__(
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
            warnings.warn('Warning: The low_bval parameter is ('
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
        input_list = ['T1w', 'dwi', 'bvec', 'bval']
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
        from clinica.utils.inputs import clinica_file_reader
        from clinica.utils.stream import cprint
        from clinica.utils.dwi import check_dwi_volume
        from clinica.utils.exceptions import ClinicaBIDSError, ClinicaException
        import clinica.utils.input_files as input_files

        all_errors = []
        try:
            t1w_files = clinica_file_reader(self.subjects,
                                            self.sessions,
                                            self.bids_directory,
                                            input_files.T1W_NII)
        except ClinicaException as e:
            all_errors.append(e)
        try:
            dwi_files = clinica_file_reader(self.subjects,
                                            self.sessions,
                                            self.bids_directory,
                                            input_files.DWI_NII)
        except ClinicaException as e:
            all_errors.append(e)

        # bval files
        try:
            bval_files = clinica_file_reader(self.subjects,
                                             self.sessions,
                                             self.bids_directory,
                                             input_files.DWI_BVAL)
        except ClinicaException as e:
            all_errors.append(e)

        # bvec files
        try:
            bvec_files = clinica_file_reader(self.subjects,
                                             self.sessions,
                                             self.bids_directory,
                                             input_files.DWI_BVEC)
        except ClinicaException as e:
            all_errors.append(e)

        if len(all_errors) > 0:
            error_message = 'Clinica faced error(s) while trying to read files in your BIDS directory.\n'
            for msg in all_errors:
                error_message += str(msg)
            raise ClinicaBIDSError(error_message)

        # Perform the check after potential issue while reading inputs
        for (dwi, bvec, bval) in zip(dwi_files, bvec_files, bval_files):
            check_dwi_volume(in_dwi=dwi, in_bvec=bvec, in_bval=bval)

        read_input_node = npe.Node(name="LoadingCLIArguments",
                                   interface=nutil.IdentityInterface(
                                       fields=self.get_input_fields(),
                                       mandatory_inputs=True),
                                   iterables=[('T1w', t1w_files),
                                              ('dwi', dwi_files),
                                              ('bvec', bvec_files),
                                              ('bval', bval_files)],
                                   synchronize=True)

        self.connect([
            (read_input_node, self.input_node, [('T1w', 'T1w')]),
            (read_input_node, self.input_node, [('dwi', 'dwi')]),
            (read_input_node, self.input_node, [('bvec', 'bvec')]),
            (read_input_node, self.input_node, [('bval', 'bval')])
        ])

    def build_output_node(self):
        """Build and connect an output node to the pipelines.
        """
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio
        from os.path import join
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
            (self.input_node, container_path,    [('dwi',                    'bids_dwi_filename')]),  # noqa
            (self.input_node,  rename_into_caps, [('dwi',                          'in_bids_dwi')]),  # noqa
            (self.output_node, rename_into_caps, [('preproc_dwi',                      'fname_dwi'),  # noqa
                                                  ('preproc_bval',                    'fname_bval'),  # noqa
                                                  ('preproc_bvec',                    'fname_bvec'),  # noqa
                                                  ('b0_mask',                  'fname_brainmask')]),  # noqa
            (container_path, write_results,      [(('container', fix_join, 'dwi'), 'container')]),  # noqa
            (rename_into_caps, write_results,    [('out_caps_dwi',    'preprocessing.@preproc_dwi'),  # noqa
                                                  ('out_caps_bval',  'preprocessing.@preproc_bval'),  # noqa
                                                  ('out_caps_bvec',  'preprocessing.@preproc_bvec'),  # noqa
                                                  ('out_caps_brainmask', 'preprocessing.@b0_mask')])  # noqa
        ])

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipelines.
        """
        import clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_workflows as workflows
        import clinica.pipelines.dwi_preprocessing_using_t1.dwi_preprocessing_using_t1_utils as utils

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.fsl as fsl

        from clinica.utils.dwi import prepare_reference_b0
        from clinica.workflows.dwi_preprocessing import ecc_pipeline
        from clinica.workflows.dwi_preprocessing import hmc_pipeline
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
        # Head-motion correction
        hmc = hmc_pipeline(name='HeadMotionCorrection')
        # Eddy-currents correction
        ecc = ecc_pipeline(name='EddyCurrentCorrection')
        # Susceptibility distortion correction using T1w image
        sdc = workflows.susceptibility_distortion_correction_using_t1(
            name='SusceptibilityDistortionCorrection')
        # Remove bias correction
        bias = remove_bias(name='RemoveBias')
        # Apply all corrections
        aac = workflows.apply_all_corrections_using_ants(
            name='ApplyAllCorrections')

        print_begin_message = npe.Node(
            interface=nutil.Function(
                input_names=['in_bids_or_caps_file'],
                function=utils.print_begin_pipeline),
            name='Write-Begin_Message')

        print_end_message = npe.Node(
            interface=nutil.Function(
                input_names=['in_bids_or_caps_file', 'final_file'],
                function=utils.print_end_pipeline),
            name='Write-End_Message')

        # Connection
        # ==========

        self.connect([
            # Print begin message
            (self.input_node, print_begin_message, [('dwi', 'in_bids_or_caps_file')]),  # noqa
            # Preliminary step (possible computation of a mean b0):
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
            (self.input_node, sdc, [('T1w',                 'inputnode.in_t1')]),  # noqa
            (ecc,             sdc, [('outputnode.out_file', 'inputnode.in_dwi')]),  # noqa
            (hmc,             sdc, [('outputnode.out_bvec', 'inputnode.in_bvec')]),  # noqa
            # Apply all corrections
            (prepare_b0,      aac, [('out_b0_dwi_merge',    'inputnode.in_dwi')]),  # noqa
            (hmc,             aac, [('outputnode.out_xfms', 'inputnode.in_hmc')]),  # noqa
            (ecc,             aac, [('outputnode.out_xfms', 'inputnode.in_ecc')]),  # noqa
            (sdc,             aac, [('outputnode.out_warp', 'inputnode.in_sdc_syb')]),  # noqa
            (self.input_node, aac, [('T1w',                 'inputnode.in_t1')]),  # noqa
            # Bias correction
            (aac, bias, [('outputnode.out_file', 'inputnode.in_file')]),
            # Outputnode:
            (bias,       self.output_node, [('outputnode.out_file', 'preproc_dwi')]),  # noqa
            (sdc,        self.output_node, [('outputnode.out_bvec', 'preproc_bvec')]),  # noqa
            (prepare_b0, self.output_node, [('out_updated_bval',    'preproc_bval')]),  # noqa
            (bias,       self.output_node, [('outputnode.b0_mask',  'b0_mask')]),  # noqa
            # Print end message
            (self.input_node, print_end_message, [('dwi',         'in_bids_or_caps_file')]),  # noqa
            (bias,            print_end_message, [('outputnode.out_file', 'final_file')]),  # noqa
        ])
