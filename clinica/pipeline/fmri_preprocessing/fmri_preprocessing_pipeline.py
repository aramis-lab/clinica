"""fMRI Preprocessing - Clinica Pipeline.
This file has been generated automatically by the `clinica generate template`
command line tool. See here for more details: https://gitlab.icm-institute.org/aramis/clinica/wikis/docs/InteractingWithClinica.
"""

# WARNING: Don't put any import statement here except if it's absolutly
# necessary. Put it *inside* the different methods.
# Otherwise it will slow down the dynamic loading of the pipelines list by the
# command line tool.
import clinica.pipeline.engine as cpe

__author__ = "Jeremy Guillon"
__copyright__ = "Copyright 2016,2017 The Aramis Lab Team"
__credits__ = ["Jeremy Guillon", "Romain Valabregue"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jeremy Guillon"
__email__ = "jeremy.guillon@inria.fr"
__status__ = "Development"


class fMRIPreprocessing(cpe.Pipeline):
    """Create fMRI preprocessing pipeline object.

    Warnings:
        - The Fieldmap node is still under revision as a pull request
        - The RealingUnwarp node is still under revision as a pull request

    Todos:
        - [x] Replace reg_node target image by the brain only using c1 + c2 + c3 dilated-eroded-filled.
        - [x] Develop SPM Realign and Unwarp wrapper and integrate it.
        - [x] Develop SPM Fieldmap Calculation Tool wrapper and integrate it.
        - [x] Replace standard DataGrabber by a BIDS tree finder.
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
        >>> from clinica.pipeline.fmri_preprocessing.fmri_preprocessing_pipeline import fMRIPreprocessing
        >>> pipeline = fMRIPreprocessing('~/MYDATASET_BIDS', '~/MYDATASET_CAPS')
        >>> pipeline.parameters = {
        >>>     'num_slices' : 45,
        >>>     'time_repetition' : 2.4,
        >>>     'echo_times' : [5.19, 7.65],
        >>>     'blip_direction' : 1,
        >>>     'total_readout_time' : 15.6799,
        >>>     'full_width_at_half_maximum' : [8, 8, 8],
        >>>     't1_native_space' : False
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

        if 't1_native_space' not in self.parameters or self.parameters[
            't1_native_space']:
            return ['t1_brain_mask', 'mc_params', 'native_fmri', 't1_fmri',
                    'mni_fmri', 'mni_smoothed_fmri']
        else:
            return ['t1_brain_mask', 'mc_params', 'native_fmri', 'mni_fmri',
                    'mni_smoothed_fmri']

    def build_input_node(self):
        """Build and connect an input node to the pipeline.
        """

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe

        # Reading BIDS
        # ============
        read_node = npe.Node(name="ReadingBIDS",
                             interface=nutil.IdentityInterface(
                                 fields=self.get_input_fields(),
                                 mandatory_inputs=True))
        # I remove the 'sub-' prefix that is not considered by the pybids'
        # layout object.
        subject_regex = '|'.join(s[4:] for s in self.subjects)
        read_node.inputs.magnitude1 = self.bids_layout.get(return_type='file',
                                                           type='magnitude1',
                                                           extensions='nii',
                                                           subject=subject_regex)
        read_node.inputs.phasediff = self.bids_layout.get(return_type='file',
                                                          type='phasediff',
                                                          extensions='nii',
                                                          subject=subject_regex)
        read_node.inputs.bold = self.bids_layout.get(return_type='file',
                                                     type='bold',
                                                     extensions='nii',
                                                     subject=subject_regex)
        read_node.inputs.T1w = self.bids_layout.get(return_type='file',
                                                    run='[1]', type='T1w',
                                                    extensions='nii',
                                                    subject=subject_regex)

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
                                 iterfield=['container']
                                           + self.get_output_fields(),
                                 interface=nio.DataSink(
                                     infields=self.get_output_fields()))
        write_node.inputs.base_directory = self.caps_directory
        write_node.inputs.parameterization = False
        write_node.inputs.container = [
            'subjects/' + self.subjects[i] + '/' + self.sessions[i] +
            '/fmri/preprocessing' for i in range(len(self.subjects))]
        write_node.inputs.remove_dest_dir = True
        write_node.inputs.regexp_substitutions = [
            (r't1_brain_mask/c3(.+)_maths_dil_ero_thresh_fillh\.nii\.gz$',
             r'\1_brainmask.nii.gz'),
            (r'mc_params/rp_a(.+)\.txt$', r'\1_confounds.tsv'),
            (r'native_fmri/ua(.+)\.nii.gz$', r'\1_space-meanBOLD.nii.gz'),
            (r't1_fmri/rua(.+)\.nii.gz$', r'\1_space-T1w.nii.gz'),
            (r'mni_fmri/wrua(.+)\.nii.gz$', r'\1_space-Ixi549Space.nii.gz'),
            (r'mni_smoothed_fmri/swrua(.+)\.nii.gz$',
             r'\1_space-Ixi549Space_fwhm-' + '-'.join(map(str, self.parameters[
                 'full_width_at_half_maximum'])) + '.nii.gz'),
            # I don't know why it's adding this empty folder, so I remove it:
            (r'trait_added', r''),
        ]

        # fMRI images in the subject's T1 native space are large, we add it
        # only if specified:
        if self.parameters['t1_native_space']:
            self.connect([
                (self.output_node, write_node, [('t1_fmri', 't1_fmri')]),
            ])

        self.connect([
            # Writing CAPS
            (self.output_node, write_node,
             [('t1_brain_mask', 't1_brain_mask')]),
            (self.output_node, write_node,
             [('mc_params', 'mc_params')]),
            (self.output_node, write_node,
             [('native_fmri', 'native_fmri')]),
            (self.output_node, write_node,
             [('mni_fmri', 'mni_fmri')]),
            (self.output_node, write_node,
             [('mni_smoothed_fmri', 'mni_smoothed_fmri')]),
        ])

    def check_custom_dependencies(self):
        pass

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
        st_node.inputs.time_acquisition = self.parameters['time_repetition'] - \
                                          self.parameters['time_repetition'] \
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
        smooth_node.inputs.fwhm = self.parameters['full_width_at_half_maximum']

        # Zipping
        # =======
        zip_node = npe.MapNode(name='Zipping',
                               iterfield=['in_file'],
                               interface=nutil.Function(input_names=['in_file'],
                                                        output_names=[
                                                            'out_file'],
                                                        function=zip_nii))
        zip_bet_node = zip_node.clone('ZippingBET')
        zip_mc_node = zip_node.clone('ZippingMC')
        zip_reg_node = zip_node.clone('ZippingRegistration')
        zip_norm_node = zip_node.clone('ZippingNormalization')
        zip_smooth_node = zip_node.clone('ZippingSmoothing')

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
            (bet_node, zip_bet_node, [('Fill.out_file', 'in_file')]),
            (mc_node, zip_mc_node, [('runwarped_files', 'in_file')]),
            (reg_node, zip_reg_node, [('coregistered_files', 'in_file')]),
            (norm_node, zip_norm_node, [('normalized_files', 'in_file')]),
            (smooth_node, zip_smooth_node, [('smoothed_files', 'in_file')]),
            # Returning output
            (zip_bet_node, self.output_node, [('out_file', 't1_brain_mask')]),
            (mc_node, self.output_node,
             [('realignment_parameters', 'mc_params')]),
            (zip_mc_node, self.output_node, [('out_file', 'native_fmri')]),
            (zip_reg_node, self.output_node, [('out_file', 't1_fmri')]),
            (zip_norm_node, self.output_node, [('out_file', 'mni_fmri')]),
            (zip_smooth_node, self.output_node, [('out_file',
                                                  'mni_smoothed_fmri')]),
        ])
