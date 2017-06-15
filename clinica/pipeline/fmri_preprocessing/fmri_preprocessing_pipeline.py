"""fMRI Preprocessing - Clinica Pipeline.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details: https://gitlab.icm-institute.org/aramis/clinica/wikis/docs/InteractingWithClinica.
"""

# WARNING: Don't put any import statement here except if it's absolutly
# necessary. Put it *inside* the different methods.
# Otherwise it will slow down the dynamic loading of the pipelines list by the
# command line tool.
import clinica.pipeline.engine as cpe


class fMRIPreprocessing(cpe.Pipeline):
    """Create fMRI preprocessing pipeline object.

    Warnings:
        - The Fieldmap node is still under revision as a pull request
        - The RealingUnwarp node is still under revision as a pull request

    Todos:
        - [x] Replace reg_node target image by the brain only using c1 + c2 + c3 dilated-eroded-filled.
        - [x] Develop SPM Realign and Unwarp wrapper and integrate it.
        - [x] Develop SPM Fieldmap Calculation Tool wrapper and integrate it.
        - [ ] Replace standard DataGrabber by a BIDS tree finder.
        - [ ] Add support of gzipped nifti inputs.

    Args:
        input_dir: A BIDS directory.
        output_dir: An empty output directory where CAPS structured data will be written.
        subjects_sessions_list: The Subjects-Sessions list file (in .tsv format).

    Returns:
        A nipype workflow object containing the full fMRI preprocessing pipeline.

    Raises:
        IOError:

    Example:
        >>> from pipelines.fmri_preprocessing import fMRIPreprocessing
        >>> pipeline = fMRIPreprocessing('~/MYDATASET_BIDS', '~/MYDATASET_CAPS')
        >>> pipeline.parameters = {
        >>>     'num_slices': 45,
        >>>     'time_repetition': 2.4,
        >>>     'echo_times': [5.19, 7.65],
        >>>     'blip_direction': 1,
        >>>     'total_readout_time': 15.6799
        >>> }
        >>> pipeline.run()
    """


    def get_input_fields(self):
        """Specify the list of possible inputs of this pipeline.

        Returns:
            A list of (string) input fields name.
        """

        return ['magnitude1', 'phasediff', 'bold', 'T1w']


    def get_output_fields(self):
        """Specify the list of possible outputs of this pipeline.

        Returns:
            A list of (string) output fields name.
        """

        return ['t1_brain_mask', 'mc_params', 'native_fmri', 't1_fmri', 'mni_fmri', 'mni_smoothed_fmri']


    def build_input_node(self):
        """Build and connect an input node to the pipeline.
        """

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        # Reading BIDS
        # ============
        read_node = npe.Node(name="ReadingBIDS",
                             interface=nutil.IdentityInterface(fields=self.get_input_fields(),
                                                               mandatory_inputs=True))
        read_node.inputs.magnitude1 = self.bids_layout.get(return_type='file', type='magnitude1', extensions='nii')
        read_node.inputs.phasediff = self.bids_layout.get(return_type='file', type='phasediff', extensions='nii')
        read_node.inputs.bold = self.bids_layout.get(return_type='file', type='bold', extensions='nii')
        read_node.inputs.T1w = self.bids_layout.get(return_type='file', run='[1]', type='T1w', extensions='nii')

        self.connect([
            # Reading BIDS
            (read_node, self.input_node, [('magnitude1', 'magnitude1')]),
            (read_node, self.input_node, [('phasediff', 'phasediff')]),
            (read_node, self.input_node, [('bold', 'bold')]),
            (read_node, self.input_node, [('T1w', 'T1w')]),
        ])


    def build_output_node(self):
        """Build and connect an output node to the pipeline.
        """

        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio

        # Writing CAPS
        # ============
        write_node = npe.MapNode(name='WritingCAPS',
                                 iterfield=['container'] + self.get_output_fields(),
                                 interface=nio.DataSink(infields=self.get_output_fields()))
        write_node.inputs.base_directory = self.caps_directory
        write_node.inputs.parameterization = False
        write_node.inputs.container = ['subjects/' + self.subjects[i] + '/' + self.sessions[i] +
                                       '/fmri/preprocessing' for i in range(len(self.subjects))]
        write_node.inputs.remove_dest_dir = True
        write_node.inputs.regexp_substitutions = [
            (r't1_brain_mask/c3(.+)_maths_dil_ero_thresh_fillh\.nii\.gz$', r'\1_brainmask.nii.gz'),
            (r'mc_params/rp_a(.+)\.txt$', r'\1_motionparams.txt'),
            (r'native_fmri/ua(.+)\.nii$', r'\1_space-native.nii'),
            (r't1_fmri/rua(.+)\.nii$', r'\1_space-t1.nii'),
            (r'mni_fmri/wrua(.+)\.nii$', r'\1_space-mni.nii'),
            (r'mni_smoothed_fmri/swrua(.+)\.nii$', r'\1_space-mni_smoothed.nii'),
            # I don't know why it's adding this empty folder, so I remove it:
            (r'trait_added', r''),
        ]

        self.connect([
            # Writing CAPS
            (self.output_node, write_node, [('t1_brain_mask', 't1_brain_mask')]),
            (self.output_node, write_node, [('mc_params', 'mc_params')]),
            (self.output_node, write_node, [('native_fmri', 'native_fmri')]),
            (self.output_node, write_node, [('t1_fmri', 't1_fmri')]),
            (self.output_node, write_node, [('mni_fmri', 'mni_fmri')]),
            (self.output_node, write_node, [('mni_smoothed_fmri', 'mni_smoothed_fmri')]),
        ])


    def build_core_nodes(self):
        """Build and connect the core nodes of the pipeline.
        """

        import fmri_preprocessing_utils as utils
        import nipype.interfaces.utility as nutil
        import nipype.interfaces.spm as spm
        import nipype.pipeline.engine as npe
        from clinica.utils.io import zip_nii, unzip_nii

        # FieldMap calculation
        # ====================
        fm_node = npe.MapNode(name="FieldMapCalculation",
                              iterfield=['phase', 'magnitude', 'epi'],
                              interface=utils.FieldMap())
        fm_node.inputs.et = self.parameters['echo_times']
        fm_node.inputs.blipdir = self.parameters['blip_direction']
        fm_node.inputs.tert = self.parameters['total_readout_time']

        # Slice timing correction
        # =======================
        st_node = npe.Node(name="SliceTimingCorrection",
                           interface=spm.SliceTiming())
        st_node.inputs.time_repetition = self.parameters['time_repetition']
        st_node.inputs.slice_order = range(1, self.parameters['num_slices'] + 1)
        st_node.inputs.num_slices = self.parameters['num_slices']
        st_node.inputs.ref_slice = self.parameters['num_slices'] / 2
        st_node.inputs.time_acquisition = self.parameters['time_repetition'] - self.parameters['time_repetition'] \
                                                                               / float(self.parameters['num_slices'])

        # Motion correction and unwarping
        # ===============================
        mc_node = npe.MapNode(name="MotionCorrectionUnwarping",
                              iterfield=["scans", "pmscan"],
                              interface=utils.RealignUnwarp())
        mc_node.inputs.register_to_mean = True
        mc_node.inputs.write_mask = False

        # Brain extraction
        # ================
        bet_node = utils.BrainExtractionWorkflow(name="BrainExtraction")

        # Registration
        # ============
        reg_node = npe.MapNode(interface=spm.Coregister(),
                               iterfield=["apply_to_files", "source", "target"],
                               name="Registration")

        # Normalization
        # =============
        norm_node = npe.MapNode(interface=spm.Normalize12(),
                                iterfield=['image_to_align', 'apply_to_files'],
                                name='Normalization')

        # Smoothing
        # =========
        smooth_node = npe.MapNode(interface=spm.Smooth(),
                                  iterfield=['in_files'],
                                  name='Smoothing')

        # Zipping
        # =======
        zip_node = npe.MapNode(name='Zipping',
                               iterfield=['in_file'],
                               interface=nutil.Function(input_names=['in_file'],
                                                        output_names=['out_file'],
                                                        function=zip_nii))

        # Connection
        # ==========
        self.connect([
            # FieldMap calculation
            (self.input_node, fm_node, [('magnitude1', 'magnitude')]),
            (self.input_node, fm_node, [('phasediff', 'phase')]),
            (self.input_node, fm_node, [('bold', 'epi')]),
            # Brain extraction
            (self.input_node, bet_node, [('T1w', 'Segmentation.data')]),
            (self.input_node, bet_node, [('T1w', 'ApplyMask.in_file')]),
            # Slice timing correction
            (self.input_node, st_node, [('bold', 'in_files')]),
            # Motion correction and unwarping
            (st_node, mc_node, [('timecorrected_files', 'scans')]),
            (fm_node, mc_node, [('vdm', 'pmscan')]),
            # Registration
            (mc_node, reg_node, [('mean_image', 'source')]),
            (mc_node, reg_node, [('runwarped_files', 'apply_to_files')]),
            (bet_node, reg_node, [('ApplyMask.out_file', 'target')]),
            # Normalization
            (self.input_node, norm_node, [('T1w', 'image_to_align')]),
            (reg_node, norm_node, [('coregistered_files', 'apply_to_files')]),
            # Smoothing
            (norm_node, smooth_node, [('normalized_files', 'in_files')]),
            # Returning output
            (bet_node, self.output_node, [('Fill.out_file', 't1_brain_mask')]),
            (mc_node, self.output_node, [('realignment_parameters', 'mc_params')]),
            (mc_node, self.output_node, [('runwarped_files', 'native_fmri')]),
            (reg_node, self.output_node, [('coregistered_files', 't1_fmri')]),
            (norm_node, self.output_node, [('normalized_files', 'mni_fmri')]),
            (smooth_node, self.output_node, [('smoothed_files', 'mni_smoothed_fmri')]),
        ])