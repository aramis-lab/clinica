# coding: utf8

from nipype import config

import clinica.pipelines.engine as cpe

# Use hash instead of parameters for iterables folder names
# Otherwise path will be too long and generate OSError
cfg = dict(execution={'parameterize_dirs': False})
config.update_config(cfg)


class DwiPreprocessingUsingT1(cpe.Pipeline):
    """DWI Preprocessing using T1 image for susceptibility distortion step.

    Ideas for improvement:
        - Replace prepare_reference_b0 function by a first run of FSL eddy
        - Replace B0-T1w registration by FA-T1w registration
        - Replace os.system(cmd) by call of Nipype interface
        - Use promising sdcflows workflows

    Warnings:
        - Do not use this pipeline if you have fieldmap data in your dataset.

    Note:
        Some reading regarding the reproducibility of FSL eddy command:
        https://www.jiscmail.ac.uk/cgi-bin/webadmin?A2=fsl;1ccf038f.1608

    Returns:
        A clinica pipeline object containing the DWIPreprocessingUsingT1 pipeline.
    """
    def check_pipeline_parameters(self):
        """Check pipeline parameters."""
        from colorama import Fore
        from clinica.utils.exceptions import ClinicaException
        from clinica.utils.stream import cprint

        if self.parameters['low_bval'] < 0:
            raise ValueError('%sThe low_bval is equals to %s: it should be zero or close to zero.%s' %
                             (Fore.RED, self.parameters['low_bval'], Fore.RESET))

        if self.parameters['low_bval'] > 100:
            cprint('%sWarning: The low_bval parameter is %s: it should be close to zero.%s' %
                   (Fore.YELLOW, self.parameters['low_bval'], Fore.RESET))

        if 'use_cuda_8_0' not in self.parameters.keys():
            self.parameters['use_cuda_8_0'] = False
        if 'use_cuda_9_1' not in self.parameters.keys():
            self.parameters['use_cuda_9_1'] = False
        if self.parameters['use_cuda_8_0'] and self.parameters['use_cuda_9_1']:
            raise ClinicaException('\n%s[Error] You must choose between CUDA 8.0 or CUDA 9.1, not both.%s' %
                                   (Fore.RED, Fore.RESET))

        if self.parameters['initrand'] not in self.parameters.keys():
            self.parameters['initrand'] = None

    def check_custom_dependencies(self):
        """Check dependencies that can not be listed in the `info.json` file."""
        pass

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Note:
            The list of inputs of the DwiPreprocessingUsingPhaseDiffFieldmap pipeline is:
                * t1w (str): Path of the T1w image in BIDS format
                * dwi (str): Path of the diffusion weighted image in BIDS format
                * dwi_json (str): Path of the DWI JSON file in BIDS format and containing
                    TotalReadoutTime and PhaseEncodingDirection metadata (see BIDS specifications)
                * bvec (str): Path of the bvec file in BIDS format
                * bval (str): Path of the bval file in BIDS format

        Returns:
            A list of (string) input fields name.
        """
        input_list = ['t1w', 'dwi', 'dwi_json', 'bvec', 'bval']

        return input_list

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Note:
            The list of outputs of the DwiPreprocessingUsingPhaseDiffFieldmap pipeline is:
                * preproc_dwi (str): Path of the preprocessed DWI
                * preproc_bvec (str): Path of the preprocessed bvec
                * preproc_bval (str): Path of the preprocessed bval
                * b0_mask (str): Path of the b0 brainmask

        Returns:
            A list of (string) output fields name.
        """
        output_list = ['preproc_dwi', 'preproc_bvec', 'preproc_bval', 'b0_mask']

        return output_list

    def build_input_node(self):
        """Build and connect an input node to the pipeline."""
        import os

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        from clinica.utils.inputs import clinica_list_of_files_reader
        from clinica.utils.input_files import T1W_NII, DWI_JSON, DWI_NII, DWI_BVEC, DWI_BVAL
        from clinica.utils.filemanip import save_participants_sessions
        from clinica.utils.stream import cprint
        from clinica.utils.ux import print_images_to_process

        list_bids_files = clinica_list_of_files_reader(
            self.subjects,
            self.sessions,
            self.bids_directory,
            [T1W_NII, DWI_JSON, DWI_NII, DWI_BVEC, DWI_BVAL],
            raise_exception=True)

        # Save subjects to process in <WD>/<Pipeline.name>/participants.tsv
        folder_participants_tsv = os.path.join(self.base_dir, self.name)
        save_participants_sessions(self.subjects, self.sessions, folder_participants_tsv)

        if len(self.subjects):
            print_images_to_process(self.subjects, self.sessions)
            cprint('List available in %s' % os.path.join(folder_participants_tsv, 'participants.tsv'))
            cprint('Computational time will depend of the number of volumes in your DWI dataset and the use of CUDA.')

        read_node = npe.Node(name="ReadingFiles",
                             iterables=[
                                 ('t1w', list_bids_files[0]),
                                 ('dwi_json', list_bids_files[1]),
                                 ('dwi', list_bids_files[2]),
                                 ('bvec', list_bids_files[3]),
                                 ('bval', list_bids_files[4]),
                             ],
                             synchronize=True,
                             interface=nutil.IdentityInterface(
                                 fields=self.get_input_fields()))
        self.connect([
            (read_node, self.input_node, [('t1w', 't1w'),
                                          ('dwi', 'dwi'),
                                          ('dwi_json', 'dwi_json'),
                                          ('bvec', 'bvec'),
                                          ('bval', 'bval')]),
        ])

    def build_output_node(self):
        """Build and connect an output node to the pipeline."""
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio
        from clinica.utils.nipype import fix_join, container_from_filename
        from .dwi_preprocessing_using_t1_utils import rename_into_caps

        # Find container path from DWI filename
        # =====================================
        container_path = npe.Node(nutil.Function(
            input_names=['bids_or_caps_filename'],
            output_names=['container'],
            function=container_from_filename),
            name='container_path')

        rename_into_caps = npe.Node(nutil.Function(
            input_names=['in_bids_dwi', 'fname_dwi', 'fname_bval',
                         'fname_bvec', 'fname_brainmask'],
            output_names=['out_caps_dwi', 'out_caps_bval', 'out_caps_bvec',
                          'out_caps_brainmask'],
            function=rename_into_caps),
            name='rename_into_caps')

        # Writing results into CAPS
        # =========================
        write_results = npe.Node(name='write_results',
                                 interface=nio.DataSink())
        write_results.inputs.base_directory = self.caps_directory
        write_results.inputs.parameterization = False

        self.connect([
            (self.input_node, container_path,    [('dwi',          'bids_or_caps_filename')]),
            (self.input_node,  rename_into_caps, [('dwi',          'in_bids_dwi')]),
            (self.output_node, rename_into_caps, [('preproc_dwi',  'fname_dwi'),
                                                  ('preproc_bval', 'fname_bval'),
                                                  ('preproc_bvec', 'fname_bvec'),
                                                  ('b0_mask',      'fname_brainmask')]),
            (container_path, write_results,      [(('container', fix_join, 'dwi'), 'container')]),
            (rename_into_caps, write_results,    [('out_caps_dwi',       'preprocessing.@preproc_dwi'),
                                                  ('out_caps_bval',      'preprocessing.@preproc_bval'),
                                                  ('out_caps_bvec',      'preprocessing.@preproc_bvec'),
                                                  ('out_caps_brainmask', 'preprocessing.@b0_mask')])
        ])

    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline."""
        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import nipype.interfaces.fsl as fsl

        from .dwi_preprocessing_using_t1_workflows import eddy_fsl_pipeline, epi_pipeline, remove_bias
        from .dwi_preprocessing_using_t1_utils import init_input_node, prepare_reference_b0, print_end_pipeline

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
        prepare_b0.inputs.low_bval = self.parameters['low_bval']
        prepare_b0.inputs.working_directory = self.base_dir
        # Mask b0 for computations purposes
        mask_b0_pre = npe.Node(fsl.BET(frac=0.3, mask=True, robust=True),
                               name='PreMaskB0')
        # Head-motion correction + Eddy-currents correction
        eddy_fsl = eddy_fsl_pipeline(low_bval=self.parameters['low_bval'],
                                     use_cuda_8_0=self.parameters['use_cuda_8_0'],
                                     use_cuda_9_1=self.parameters['use_cuda_9_1'],
                                     seed_fsl_eddy=self.parameters['initrand'])
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
            (init_node, print_end_message, [('image_id', 'image_id')]),
            (bias,      print_end_message, [('outputnode.out_file', 'final_file')]),
            # Output node
            (bias,       self.output_node, [('outputnode.out_file', 'preproc_dwi')]),
            (sdc,        self.output_node, [('outputnode.out_bvec', 'preproc_bvec')]),
            (prepare_b0, self.output_node, [('out_updated_bval',    'preproc_bval')]),
            (bias,       self.output_node, [('outputnode.b0_mask',  'b0_mask')])
        ])
