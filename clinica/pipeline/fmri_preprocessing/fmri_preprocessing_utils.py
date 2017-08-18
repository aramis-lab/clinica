"""

"""


from __future__ import print_function, division, unicode_literals, absolute_import
from nipype.utils.filemanip import (fname_presuffix, filename_to_list, list_to_filename)
from nipype.interfaces.base import (OutputMultiPath, TraitedSpec, isdefined, traits, File)
from nipype.interfaces.spm.base import (SPMCommand, scans_for_fname, func_is_3d, scans_for_fnames, SPMCommandInputSpec)
import nipype.pipeline.engine as npe
import nipype.interfaces.spm as spm
import nipype.interfaces.fsl as fsl


class RealignUnwarpInputSpec(SPMCommandInputSpec):
    scans = traits.Either(traits.List(File(exists=True)), File(exists=True),
        mandatory=True, copyfile=True,
        field='data.scans',
        desc='scans for this session')
    pmscan = File(mandatory=True, exists=True, copyfile=False,
        field='data.pmscan',
        desc='phase map (vdm* file)')
    # Estimation Options
    quality = traits.Range(low=0.0, high=1.0,
        field='eoptions.quality',
        desc='0.1 = fast, 1.0 = precise')
    fwhm = traits.Range(low=0.0,
        field='eoptions.fwhm',
        desc='gaussian smoothing kernel width')
    separation = traits.Range(low=0.0,
        field='eoptions.sep',
        desc='sampling separation in mm')
    register_to_mean = traits.Bool(
        field='eoptions.rtm',
        desc='Indicate whether realignment is done to the mean image')
    weight_img = File(exists=True,
        field='eoptions.weight',
        desc='filename of weighting image')
    interp = traits.Range(low=0, high=7,
        field='eoptions.interp',
        desc='degree of b-spline used for interpolation')
    wrap = traits.List(traits.Int(), minlen=3, maxlen=3,
        field='eoptions.wrap',
        desc='Check if interpolation should wrap in [x,y,z]')
    # Unwarp Options
    basfcn = traits.List(traits.Int(), minlen=2, maxlen=2,
        field='uweoptions.basfcn',
        desc='basic functions')
    regorder = traits.Range(low=0, high=3,
        field='uweoptions.regorder',
        desc='regularization')
    regfactor = traits.Range(low=0,
        field='uweoptions.lambda',
        desc='regularization factor')
    jm = traits.Bool(
        field='uweoptions.jm',
        desc='Jacobian deformations')
    fot = traits.List(traits.Int(), minlen=2, maxlen=2,
        field='uweoptions.fot',
        desc='first-order effects')
    sot = traits.List(traits.Int(), minlen=2, maxlen=2,
        field='uweoptions.sot',
        desc='second-order effects')
    uwfwhm = traits.Range(low=0,
        field='uweoptions.uwfwhm',
        desc='smoothing for unwarp')
    rem = traits.Bool(
        field='uweoptions.rem',
        desc='')
    noi = traits.Range(low=0,
        field='uweoptions.noi',
        desc='maximum number of iterations')
    expround = traits.Enum('Average', 'First', 'Last',
        field='uweoptions.expround',
        desc='one of: Average, First, Last')
    # Reslice Options
    write_which = traits.ListInt([2, 1],
        minlen=2, maxlen=2, usedefault=True,
        field='uwroptions.uwwhich',
        desc='determines which images to reslice')
    write_interp = traits.Range(low=0, high=7,
        field='uwroptions.rinterp',
        desc='degree of b-spline used for interpolation')
    write_wrap = traits.List(traits.Int(), minlen=3, maxlen=3,
        field='uwroptions.wrap',
        desc='Check if interpolation should wrap in [x,y,z]')
    write_mask = traits.Bool(
        field='uwroptions.mask',
        desc='True/False mask output image')
    out_prefix = traits.String('u', usedefault=True,
        field='uwroptions.prefix',
        desc='realigned output prefix')


class RealignUnwarpOutputSpec(TraitedSpec):
    mean_image = File(exists=True, desc='Mean image file from the realignment')
    modified_scans = traits.Either(
        traits.List(File(exists=True)), File(exists=True),
        desc=('Copies of all files passed to '
              'scans. Headers will have '
              'been modified to align all '
              'images with the first, or '
              'optionally to first do that, '
              'extract a mean image, and '
              're-align to that mean image.'))
    runwarped_files = traits.Either(
        traits.List(File(exists=True)), File(exists=True),
        desc='These will be the resliced and unwarped files.')
    realignment_parameters = OutputMultiPath(File(exists=True),
                                             desc=('Estimated translation and '
                                                   'rotation parameters'))


class RealignUnwarp(SPMCommand):
    """Use spm_realign for estimating within modality rigid body alignment

    http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf#page=25

    Examples
    --------

    >>> import nipype.interfaces.spm as spm
    >>> runwarp = spm.RealignUnwarp()
    >>> runwarp.inputs.scans = 'functional.nii'
    >>> runwarp.inputs.pmscan = 'functional.nii'
    >>> runwarp.inputs.register_to_mean = True
    >>> runwarp.run() # doctest: +SKIP

    """

    input_spec = RealignUnwarpInputSpec
    output_spec = RealignUnwarpOutputSpec

    _jobtype = 'spatial'
    _jobname = 'realignunwarp'

    def _format_arg(self, opt, spec, val):
        """Convert input to appropriate format for spm
        """
        if opt == 'scans' or opt == 'pmscan':
            return scans_for_fnames(filename_to_list(val))
        return super(RealignUnwarp, self)._format_arg(opt, spec, val)

    def _list_outputs(self):
        outputs = self._outputs().get()
        resliced_all = self.inputs.write_which[0] > 0
        resliced_mean = self.inputs.write_which[1] > 0

        if isdefined(self.inputs.scans):
            outputs['realignment_parameters'] = []
        for imgf in filename_to_list(self.inputs.scans):
            if isinstance(imgf, list):
                tmp_imgf = imgf[0]
            else:
                tmp_imgf = imgf
            outputs['realignment_parameters'].append(
                fname_presuffix(tmp_imgf, prefix='rp_', suffix='.txt',
                                use_ext=False))
            if not isinstance(imgf, list) and func_is_3d(imgf):
                break

        outputs['modified_scans'] = self.inputs.scans


        first_image = filename_to_list(self.inputs.scans)[0]
        if resliced_mean:
            outputs['mean_image'] = fname_presuffix(first_image,
                                                    prefix=('mean'+self.inputs.out_prefix))

        if resliced_all:
            outputs['runwarped_files'] = []
            for idx, imgf in enumerate(
                    filename_to_list(self.inputs.scans)):
                runwarped_run = []
                if isinstance(imgf, list):
                    for i, inner_imgf in enumerate(filename_to_list(imgf)):
                        newfile = fname_presuffix(
                            inner_imgf, prefix=self.inputs.out_prefix)
                        runwarped_run.append(newfile)
                else:
                    runwarped_run = fname_presuffix(
                        imgf, prefix=self.inputs.out_prefix)
                outputs['runwarped_files'].append(runwarped_run)
        return outputs


class FieldMapInputSpec(SPMCommandInputSpec):
    jobtype = traits.Enum('calculatevdm', 'applyvdm', usedefault=True,
        desc='one of: calculatevdm, applyvdm')
    phase = File(mandatory=True, exists=True, copyfile=False,
        field='subj.data.presubphasemag.phase',
        desc='presubstracted phase file')
    magnitude = File(mandatory=True, exists=True, copyfile=False,
        field='subj.data.presubphasemag.magnitude',
        desc='presubstracted magnitude file')
    et = traits.List(traits.Float(), minlen=2, maxlen=2, mandatory=True,
        field='subj.defaults.defaultsval.et',
        desc='short and long echo times')
    maskbrain = traits.Bool(True, usedefault=True,
        field='subj.defaults.defaultsval.maskbrain',
        desc='masking or no masking of the brain')
    blipdir = traits.Enum(1, -1, mandatory=True,
        field='subj.defaults.defaultsval.blipdir',
        desc='polarity of the phase-encode blips')
    tert = traits.Float(mandatory=True,
        field='subj.defaults.defaultsval.tert',
        desc='total EPI readout time')
    epifm = traits.Bool(False, usedefault=True,
        field='subj.defaults.defaultsval.epifm',
        desc='epi-based field map');
    ajm = traits.Bool(False, usedefault=True,
        field='subj.defaults.defaultsval.ajm',
        desc='jacobian modulation');
    # Unwarping defaults parameters
    method = traits.Enum('Mark3D', 'Mark2D', 'Huttonish', usedefault=True,
        desc='One of: Mark3D, Mark2D, Huttonish',
        field='subj.defaults.defaultsval.uflags.method');
    fwhm = traits.Range(low=0, value=10, usedefault=True,
        field='subj.defaults.defaultsval.uflags.fwhm',
        desc='gaussian smoothing kernel width');
    pad = traits.Range(low=0, value=0, usedefault=True,
        field='subj.defaults.defaultsval.uflags.pad',
        desc='padding kernel width');
    ws = traits.Bool(True, usedefault=True,
        field='subj.defaults.defaultsval.uflags.ws',
        desc='weighted smoothing');
    # Brain mask defaults parameters
    template = traits.File(copyfile=False, exists=True,
        field='subj.defaults.defaultsval.mflags.template',
        desc='template image for brain masking');
    fwhm = traits.Range(low=0, value=5, usedefault=True,
        field='subj.defaults.defaultsval.mflags.fwhm',
        desc='gaussian smoothing kernel width');
    nerode = traits.Range(low=0, value=2, usedefault=True,
        field='subj.defaults.defaultsval.mflags.nerode',
        desc='number of erosions');
    ndilate = traits.Range(low=0, value=4, usedefault=True,
        field='subj.defaults.defaultsval.mflags.ndilate',
        desc='number of erosions');
    thresh = traits.Float(0.5, usedefault=True,
        field='subj.defaults.defaultsval.mflags.thresh',
        desc='threshold used to create brain mask from segmented data');
    reg = traits.Float(0.02, usedefault=True,
        field='subj.defaults.defaultsval.mflags.reg',
        desc='regularization value used in the segmentation');
    # EPI unwarping for quality check
    epi = traits.File(copyfile=False, exists=True, mandatory=True,
        field='subj.session.epi',
        desc='EPI to unwarp');
    matchvdm = traits.Bool(True, usedefault=True,
        field='subj.matchvdm',
        desc='match VDM to EPI');
    sessname = traits.String('_run-', usedefault=True,
        field='subj.sessname',
        desc='VDM filename extension');
    writeunwarped = traits.Bool(False, usedefault=True,
        field='subj.writeunwarped',
        desc='write unwarped EPI');
    anat = traits.File(copyfile=False, exists=True,
        field='subj.anat',
        desc='anatomical image for comparison');
    matchanat = traits.Bool(True, usedefault=True,
        field='subj.matchanat',
        desc='match anatomical image to EPI');


class FieldMapOutputSpec(TraitedSpec):
    vdm = File(exists=True, desc='voxel difference map')


class FieldMap(SPMCommand):
    """Use spm to calculate fieldmap vdm.

    http://www.fil.ion.ucl.ac.uk/spm/doc/manual.pdf#page=19

    To do
    -----
    Deal with real/imag magnitude images and with the two phase files case.

    Examples
    --------
    >>> from nipype.interfaces.spm import FieldMap
    >>> fm = FieldMap()
    >>> fm.inputs.phase = 'phasediff.nii'
    >>> fm.inputs.magnitude = 'magnitude1.nii'
    >>> fm.inputs.et = [5.19, 7.65]
    >>> fm.inputs.blipdir = 1
    >>> fm.inputs.tert = 15.6
    >>> fm.inputs.epi = 'bold.nii'
    >>> fm.run() # doctest: +SKIP

    """

    input_spec = FieldMapInputSpec
    output_spec = FieldMapOutputSpec
    _jobtype = 'tools'
    _jobname = 'fieldmap'

    def _format_arg(self, opt, spec, val):
        """Convert input to appropriate format for spm
        """
        if opt == 'phase' or opt == 'magnitude' or opt == 'anat':
            return scans_for_fname(filename_to_list(val))
        if opt == 'epi' or opt == 'magnitude':
            return scans_for_fname(filename_to_list(val))

        return super(FieldMap, self)._format_arg(opt, spec, val)

    def _parse_inputs(self):
        """validate spm fieldmap options if set to None ignore
        """
        einputs = super(FieldMap, self)._parse_inputs()
        jobtype = self.inputs.jobtype
        return [{'%s' % (jobtype): einputs[0]}]

    def _list_outputs(self):
        outputs = self._outputs().get()
        jobtype = self.inputs.jobtype
        if jobtype == "calculatevdm":
            outputs['vdm'] = []
            for phase in filename_to_list(self.inputs.phase):
                outputs['vdm'].append(fname_presuffix(phase, prefix='vdm5_sc'))
            outputs['vdm'] = list_to_filename(outputs['vdm'])

        return outputs


class BrainExtractionWorkflow(npe.Workflow):

    def __init__(self, name, base_dir=None):

        super(BrainExtractionWorkflow, self).__init__(name, base_dir)

        # Segmentation
        # ============
        seg_node = npe.MapNode(name="Segmentation",
                               iterfield="data",
                               interface=spm.Segment())
        seg_node.inputs.gm_output_type = [False, False, True]
        seg_node.inputs.wm_output_type = [False, False, True]
        seg_node.inputs.csf_output_type = [False, False, True]
        add1_node = npe.MapNode(name="AddGMWM",
                                iterfield=["in_file", "operand_file"],
                                interface=fsl.BinaryMaths())
        add1_node.inputs.operation = 'add'
        add2_node = npe.MapNode(name="AddGMWMCSF",
                                iterfield=["in_file", "operand_file"],
                                interface=fsl.BinaryMaths())
        add2_node.inputs.operation = 'add'
        dil_node = npe.MapNode(name="Dilate",
                               iterfield="in_file",
                               interface=fsl.DilateImage())
        dil_node.inputs.operation = 'mean'
        ero_node = npe.MapNode(name="Erode",
                               iterfield="in_file",
                               interface=fsl.ErodeImage())
        thre_node = npe.MapNode(name="Threshold",
                                iterfield="in_file",
                                interface=fsl.Threshold())
        thre_node.inputs.thresh = 0.5
        fill_node = npe.MapNode(name="Fill",
                                iterfield="in_file",
                                interface=fsl.UnaryMaths())
        fill_node.inputs.operation = 'fillh'
        mask_node = npe.MapNode(name="ApplyMask",
                                iterfield=["in_file", "mask_file"],
                                interface=fsl.ApplyMask())
        mask_node.inputs.output_type = str("NIFTI")

        self.connect([
            (seg_node   , add1_node   , [('native_gm_image' , 'in_file'     )]),
            (seg_node   , add1_node   , [('native_wm_image' , 'operand_file')]),
            (seg_node   , add2_node   , [('native_csf_image', 'in_file'     )]),
            (add1_node  , add2_node   , [('out_file'        , 'operand_file')]),
            (add2_node  , dil_node    , [('out_file'        , 'in_file'     )]),
            (dil_node   , ero_node    , [('out_file'        , 'in_file'     )]),
            (ero_node   , thre_node   , [('out_file'        , 'in_file'     )]),
            (thre_node  , fill_node   , [('out_file'        , 'in_file'     )]),
            (fill_node  , mask_node   , [('out_file'        , 'mask_file'   )]),
        ])
