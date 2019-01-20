# coding: utf8

# WARNING: Don't put any import statement here except if it's absolutly
# necessary. Put it *inside* the different methods.
# Otherwise it will slow down the dynamic loading of the pipelines list by the
# command line tool.
import clinica.pipelines.engine as cpe

__author__ = "Jeremy Guillon"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Jeremy Guillon", "Romain Valabregue"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jeremy Guillon"
__email__ = "jeremy.guillon@inria.fr"
__status__ = "Development"


class fMRIPreprocessing(cpe.Pipeline):
    """Create fMRI preprocessing pipelines object.

    Warnings:
        - The Fieldmap node is still under revision as a pull request
        - The RealingUnwarp node is still under revision as a pull request

    Todos:
        - [x] Don't read inputs if not needed (i.e. --unwarp or no)
        - [x] Read parameters from sidecar `*.json` files.
        - [x] Add support of gzipped nifti inputs.
        - [x] Replace reg_node target image by the brain only using c1 + c2 +
        c3 dilated-eroded-filled.
        - [x] Develop SPM Realign and Unwarp wrapper and integrate it.
        - [x] Develop SPM Fieldmap Calculation Tool wrapper and integrate it.
        - [x] Replace standard DataGrabber by a BIDS tree finder.
        - [x] Export only gzipped nifti files.

    Args:
        input_dir: A BIDS directory.
        output_dir: An empty output directory where CAPS structured data will
        be written.
        subjects_sessions_list: The Subjects-Sessions list file (in .tsv
        format).

    Returns:
        A nipype workflow object containing the full fMRI preprocessing
        pipelines.

    Raises:
        IOError:

    Example:
        >>> from clinica.pipelines.fmri_preprocessing
        .fmri_preprocessing_pipeline import fMRIPreprocessing
        >>> pipelines = fMRIPreprocessing('~/MYDATASET_BIDS',
        '~/MYDATASET_CAPS')
        >>> pipelines.parameters = {
        >>>     'num_slices' : 45,
        >>>     'time_repetition' : 2.4,
        >>>     'echo_times' : [5.19, 7.65],
        >>>     'blip_direction' : 1,
        >>>     'total_readout_time' : 15.6799,
        >>>     'full_width_at_half_maximum' : [8, 8, 8],
        >>>     't1_native_space' : False
        >>> }
        >>> pipelines.run()
    """

    def get_input_fields(self):
        """Specify the list of possible inputs of this pipelines.

        Returns:
            A list of (string) input fields name.
        """

        if ('unwarping' in self.parameters) and self.parameters['unwarping']:
            return ['et', 'blipdir', 'tert', 'time_repetition', 'num_slices',
                    'magnitude1', 'slice_order', 'ref_slice',
                    'time_acquisition', 'phasediff', 'bold', 'T1w']
        else:
            return ['time_repetition', 'num_slices', 'slice_order', 'ref_slice',
                    'time_acquisition', 'bold', 'T1w']

    def get_output_fields(self):
        """Specify the list of possible outputs of this pipelines.

        Returns:
            A list of (string) output fields name.
        """

        if ('t1_native_space' in self.parameters) and self.parameters[
                't1_native_space']:
            return ['t1_brain_mask', 'mc_params', 'native_fmri', 't1_fmri',
                    'mni_fmri', 'mni_smoothed_fmri']
        else:
            return ['t1_brain_mask', 'mc_params', 'native_fmri', 'mni_fmri',
                    'mni_smoothed_fmri']

    def build_input_node(self):
        """Build and connect an input node to the pipelines.

        References:
            https://lcni.uoregon.edu/kb-articles/kb-0003

        """

        import nipype.interfaces.utility as nutil
        import nipype.pipeline.engine as npe
        import json
        import numpy as np
        from clinica.utils.stream import cprint

        # Reading BIDS files
        # ==================
        read_node = npe.Node(name="ReadingBIDS",
                             interface=nutil.IdentityInterface(
                                     fields=self.get_input_fields(),
                                     mandatory_inputs=True))
        # I remove the 'sub-' prefix that is not considered by the pybids'
        # layout object.
        subject_regex = '|'.join(s[4:] for s in self.subjects)

        if ('unwarping' in self.parameters) and self.parameters['unwarping']:
            read_node.inputs.magnitude1 = self.bids_layout.get(
                    return_type='file',
                    type='magnitude1',
                    extensions='nii.gz',
                    subject=subject_regex)
            read_node.inputs.phasediff = self.bids_layout.get(
                    return_type='file',
                    type='phasediff',
                    extensions='nii.gz',
                    subject=subject_regex)
        read_node.inputs.bold = self.bids_layout.get(return_type='file',
                                                     type='bold',
                                                     extensions='nii.gz',
                                                     subject=subject_regex)
        read_node.inputs.T1w = self.bids_layout.get(return_type='file',
                                                    type='T1w',
                                                    extensions='nii.gz',
                                                    subject=subject_regex)

        # Reading BIDS json
        # =================

        read_node.inputs.et = []
        read_node.inputs.blipdir = []
        read_node.inputs.tert = []
        read_node.inputs.time_repetition = []
        read_node.inputs.num_slices = []
        read_node.inputs.slice_order = []
        read_node.inputs.ref_slice = []
        read_node.inputs.time_acquisition = []

        for i in range(len(self.subjects)):
            cprint('Loading subject "{sub}"...'.format(sub=self.subjects[i]))

            if self.parameters['unwarping']:
                # From phasediff json file
                phasediff_json = self.bids_layout.get(return_type='file',
                                                      type='phasediff',
                                                      extensions='json',
                                                      subject=self.subjects[i][
                                                              4:])
                with open(phasediff_json[0]) as json_file:
                    data = json.load(json_file)
                    # SPM echo times
                    read_node.inputs.et.append([data['EchoTime1'],
                                                data['EchoTime2']])
                    # SPM blip direction
                    # TODO: Verifiy that it is the correct way to get the
                    # blipdir
                    blipdir_raw = data['PhaseEncodingDirection']
                    if len(blipdir_raw) > 1 and blipdir_raw[1] == '-':
                        read_node.inputs.blipdir.append(-1)
                    else:
                        read_node.inputs.blipdir.append(1)

            # From func json file
            func_json = self.bids_layout.get(return_type='file',
                                             type='bold',
                                             extensions='json',
                                             subject=self.subjects[i][4:])
            with open(func_json[0]) as json_file:
                data = json.load(json_file)
                # SPM Total readout time
                read_node.inputs.tert.append(
                        1 / data['BandwidthPerPixelPhaseEncode'])
                # SPM Repetition time
                read_node.inputs.time_repetition.append(data['RepetitionTime'])
                # Number of slices
                slice_timing = data['SliceTiming']
                read_node.inputs.num_slices.append(len(slice_timing))
                # Slice order
                slice_order = np.argsort(slice_timing) + 1
                read_node.inputs.slice_order.append(slice_order.tolist())
                read_node.inputs.ref_slice.append(np.argmin(slice_timing) + 1)
                read_node.inputs.time_acquisition.append(
                        data['RepetitionTime'] - data['RepetitionTime']
                        / float(len(slice_timing)))

            cprint(read_node.inputs)

        if ('unwarping' in self.parameters) and self.parameters['unwarping']:
            self.connect([
                # Reading BIDS json
                (read_node, self.input_node, [('et', 'et')]),
                (read_node, self.input_node, [('blipdir', 'blipdir')]),
                (read_node, self.input_node, [('tert', 'tert')]),
                # Reading BIDS files
                (read_node, self.input_node, [('phasediff', 'phasediff')]),
                (read_node, self.input_node, [('magnitude1', 'magnitude1')]),
            ])

        self.connect([
            # Reading BIDS json
            (read_node, self.input_node,
             [('time_repetition', 'time_repetition')]),
            (read_node, self.input_node, [('num_slices', 'num_slices')]),
            (read_node, self.input_node, [('slice_order', 'slice_order')]),
            (read_node, self.input_node, [('ref_slice', 'ref_slice')]),
            (read_node, self.input_node,
             [('time_acquisition', 'time_acquisition')]),
            # Reading BIDS files
            (read_node, self.input_node, [('bold', 'bold')]),
            (read_node, self.input_node, [('T1w', 'T1w')]),
        ])

    def build_output_node(self):
        """Build and connect an output node to the pipelines.
        """

        import nipype.pipeline.engine as npe
        import nipype.interfaces.io as nio

        # Writing CAPS
        # ============
        write_node = npe.MapNode(name='WritingCAPS',
                                 iterfield=['container'] + self.get_output_fields(),
                                 interface=nio.DataSink(
                                         infields=self.get_output_fields()))
        write_node.inputs.base_directory = self.caps_directory
        write_node.inputs.parameterization = False
        write_node.inputs.container = [
            'subjects/' + self.subjects[i] + '/' + self.sessions[i] +
            '/fmri/preprocessing' for i in range(len(self.subjects))]
        write_node.inputs.remove_dest_dir = True
        if ('freesurfer_brain_mask' in self.parameters) and \
                not (self.parameters['freesurfer_brain_mask']):
            write_node.inputs.regexp_substitutions = [
                (r't1_brain_mask/c3(.+)_maths_dil_ero_thresh_fillh\.nii\.gz$',
                 r'\1_brainmask.nii.gz'),
                (r'mc_params/rp_a(.+)\.txt$', r'\1_motion.tsv'),
                (r'native_fmri/[u|r]a(.+)\.nii.gz$',
                 r'\1_space-meanBOLD_preproc.nii.gz'),
                (r't1_fmri/r[u|r]a(.+)\.nii.gz$',
                 r'\1_space-T1w_preproc.nii.gz'),
                (r'mni_fmri/wr[u|r]a(.+)\.nii.gz$',
                 r'\1_space-Ixi549Space_preproc.nii.gz'),
                (r'mni_smoothed_fmri/swr[u|r]a(.+)\.nii.gz$',
                 r'\1_space-Ixi549Space_fwhm-' + 'x'.join(
                     map(str, self.parameters[
                         'full_width_at_half_maximum'])) + '_preproc.nii.gz'),
                # I don't know why it's adding this empty folder, so I remove
                # it:
                (r'trait_added', r''),
            ]
        else:
            write_node.inputs.regexp_substitutions = [
                (r't1_brain_mask/(.+)\.nii\.gz$', r'\1_brainmask.nii.gz'),
                (r'mc_params/rp_a(.+)\.txt$', r'\1_motion.tsv'),
                (r'native_fmri/[u|r]a(.+)\.nii.gz$',
                 r'\1_space-meanBOLD_preproc.nii.gz'),
                (r't1_fmri/r[u|r]a(.+)\.nii.gz$',
                 r'\1_space-T1w_preproc.nii.gz'),
                (r'mni_fmri/wr[u|r]a(.+)\.nii.gz$',
                 r'\1_space-Ixi549Space_preproc.nii.gz'),
                (r'mni_smoothed_fmri/swr[u|r]a(.+)\.nii.gz$',
                 r'\1_space-Ixi549Space_fwhm-' + 'x'.join(
                     map(str, self.parameters[
                         'full_width_at_half_maximum'])) +
                 '_preproc.nii.gz'),
                # I don't know why it's adding this empty folder, so I remove
                # it:
                (r'trait_added', r''),
            ]

        # fMRI images in the subject's T1 native space are large, we add it
        # only if specified:
        if ('t1_native_space' in self.parameters) and self.parameters[
                't1_native_space']:
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
        """Build and connect the core nodes of the pipelines.
        """

        import clinica.pipelines.fmri_preprocessing.fmri_preprocessing_utils as utils
        import nipype.interfaces.utility as nutil
        import nipype.interfaces.spm as spm
        import nipype.pipeline.engine as npe
        from clinica.utils.io import zip_nii, unzip_nii

        # Zipping
        # =======
        unzip_node = npe.MapNode(name='Unzipping',
                                 iterfield=['in_file'],
                                 interface=nutil.Function(
                                         input_names=['in_file'],
                                         output_names=[
                                             'out_file'],
                                         function=unzip_nii))

        unzip_T1w = unzip_node.clone('UnzippingT1w')
        unzip_phasediff = unzip_node.clone('UnzippingPhasediff')
        unzip_bold = unzip_node.clone('UnzippingBold')
        unzip_magnitude1 = unzip_node.clone('UnzippingMagnitude1')

        # FieldMap calculation
        # ====================
        if self.parameters['unwarping']:
            fm_node = npe.MapNode(name="FieldMapCalculation",
                                  iterfield=['phase', 'magnitude', 'epi',
                                             'et', 'blipdir', 'tert'],
                                  interface=spm.FieldMap())

        # Slice timing correction
        # =======================
        st_node = npe.MapNode(name="SliceTimingCorrection",
                              iterfield=['in_files', 'time_repetition',
                                         'slice_order', 'num_slices',
                                         'ref_slice', 'time_acquisition'],
                              interface=spm.SliceTiming())

        # Motion correction and unwarping
        # ===============================

        if self.parameters['unwarping']:
            mc_node = npe.MapNode(name="MotionCorrectionUnwarping",
                                  iterfield=["scans", "pmscan"],
                                  interface=spm.RealignUnwarp())
            mc_node.inputs.register_to_mean = True
            mc_node.inputs.reslice_mask = False
        else:
            mc_node = npe.MapNode(name="MotionCorrection",
                                  iterfield=["in_files"],
                                  interface=spm.Realign())
            mc_node.inputs.register_to_mean = True

        # Brain extraction
        # ================
        import os.path as path
        from nipype.interfaces.freesurfer import MRIConvert
        if self.parameters['freesurfer_brain_mask']:
            brain_masks = [path.join(self.caps_directory, 'subjects',
                                     self.subjects[i], self.sessions[i],
                                     't1/freesurfer_cross_sectional',
                                     self.subjects[i] + '_' + self.sessions[i],
                                     'mri/brain.mgz')
                           for i in range(len(self.subjects))]
            conv_brain_masks = [str(self.subjects[i] + '_' + self.sessions[i] +
                                    '.nii')
                                for i in range(len(self.subjects))]
            bet_node = npe.MapNode(interface=MRIConvert(),
                                   iterfield=["in_file", "out_file"],
                                   name="BrainConversion")
            bet_node.inputs.in_file = brain_masks
            bet_node.inputs.out_file = conv_brain_masks
            bet_node.inputs.out_type = 'nii'
        else:
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

        # Connections
        # ===========

        if self.parameters['freesurfer_brain_mask']:
            self.connect([
                # Brain extraction
                (bet_node, reg_node, [('out_file', 'target')]),
                (bet_node, zip_bet_node, [('out_file', 'in_file')]),
            ])
        else:
            self.connect([
                # Brain extraction
                (unzip_T1w, bet_node, [('out_file', 'Segmentation.data')]),
                (unzip_T1w, bet_node, [('out_file', 'ApplyMask.in_file')]),
                (bet_node, reg_node, [('ApplyMask.out_file', 'target')]),
                (bet_node, zip_bet_node, [('Fill.out_file', 'in_file')]),
            ])

        if self.parameters['unwarping']:
            self.connect([
                # FieldMap calculation
                (self.input_node, fm_node, [('et', 'et')]),
                (self.input_node, fm_node, [('blipdir', 'blipdir')]),
                (self.input_node, fm_node, [('tert', 'tert')]),
                (self.input_node, unzip_phasediff, [('phasediff', 'in_file')]),
                (self.input_node, unzip_magnitude1,
                 [('magnitude1', 'in_file')]),
                (unzip_magnitude1, fm_node, [('out_file', 'magnitude')]),
                (unzip_phasediff, fm_node, [('out_file', 'phase')]),
                (unzip_bold, fm_node, [('out_file', 'epi')]),
                # Motion correction and unwarping
                (st_node, mc_node, [('timecorrected_files', 'scans')]),
                (fm_node, mc_node, [('vdm', 'pmscan')]),
                (mc_node, reg_node, [('realigned_unwarped_files', 'apply_to_files')]),
                (mc_node, zip_mc_node, [('realigned_unwarped_files', 'in_file')]),
            ])
        else:
            self.connect([
                # Motion correction and unwarping
                (st_node, mc_node, [('timecorrected_files', 'in_files')]),
                (mc_node, reg_node, [('realigned_files', 'apply_to_files')]),
                (mc_node, zip_mc_node, [('realigned_files', 'in_file')]),
            ])
        self.connect([
            # Unzipping
            (self.input_node, unzip_T1w, [('T1w', 'in_file')]),
            (self.input_node, unzip_bold, [('bold', 'in_file')]),
            # Slice timing correction
            (unzip_bold, st_node, [('out_file', 'in_files')]),
            (self.input_node, st_node,
             [('time_repetition', 'time_repetition')]),
            (self.input_node, st_node, [('num_slices', 'num_slices')]),
            (self.input_node, st_node, [('slice_order', 'slice_order')]),
            (self.input_node, st_node, [('ref_slice', 'ref_slice')]),
            (self.input_node, st_node,
             [('time_acquisition', 'time_acquisition')]),
            # Registration
            (mc_node, reg_node, [('mean_image', 'source')]),
            # Normalization
            (unzip_T1w, norm_node, [('out_file', 'image_to_align')]),
            (reg_node, norm_node, [('coregistered_files', 'apply_to_files')]),
            # Smoothing
            (norm_node, smooth_node, [('normalized_files', 'in_files')]),
            # Zipping
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
