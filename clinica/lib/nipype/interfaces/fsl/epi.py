# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""The fsl module provides classes for interfacing with the `FSL
<http://www.fmrib.ox.ac.uk/fsl/index.html>`_ command line tools.  This
was written to work with FSL version 5.0.4.
"""
from __future__ import print_function, division, unicode_literals, \
    absolute_import
from builtins import str

import os

from nipype.interfaces.base import (traits, TraitedSpec, File, isdefined)
from nipype.interfaces.fsl.base import FSLCommand, FSLCommandInputSpec


class EddyInputSpec(FSLCommandInputSpec):
    in_file = File(
        exists=True,
        mandatory=True,
        argstr='--imain=%s',
        desc=('File containing all the images to estimate '
              'distortions for'))
    in_mask = File(
        exists=True,
        mandatory=True,
        argstr='--mask=%s',
        desc='Mask to indicate brain')
    in_index = File(
        exists=True,
        mandatory=True,
        argstr='--index=%s',
        desc=('File containing indices for all volumes in --imain '
              'into --acqp and --topup'))
    in_acqp = File(
        exists=True,
        mandatory=True,
        argstr='--acqp=%s',
        desc='File containing acquisition parameters')
    in_bvec = File(
        exists=True,
        mandatory=True,
        argstr='--bvecs=%s',
        desc=('File containing the b-vectors for all volumes in '
              '--imain'))
    in_bval = File(
        exists=True,
        mandatory=True,
        argstr='--bvals=%s',
        desc=('File containing the b-values for all volumes in '
              '--imain'))
    out_base = traits.Str(
        'eddy_corrected',
        argstr='--out=%s',
        usedefault=True,
        desc=('basename for output (warped) image'))
    session = File(
        exists=True,
        argstr='--session=%s',
        desc=('File containing session indices for all volumes in '
              '--imain'))
    in_topup_fieldcoef = File(
        exists=True,
        argstr="--topup=%s",
        requires=['in_topup_movpar'],
        desc=('topup file containing the field '
              'coefficients'))
    in_topup_movpar = File(
        exists=True,
        requires=['in_topup_fieldcoef'],
        desc='topup movpar.txt file')

    flm = traits.Enum(
        'linear',
        'quadratic',
        'cubic',
        argstr='--flm=%s',
        desc='First level EC model')

    slm = traits.Enum(
        'none',
        'linear',
        'quadratic',
        argstr='--slm=%s',
        desc='Second level EC model')

    fep = traits.Bool(
        False, argstr='--fep', desc='Fill empty planes in x- or y-directions')

    interp = traits.Enum(
        'spline',
        'trilinear',
        argstr='--interp=%s',
        desc='Interpolation model for estimation step')

    nvoxhp = traits.Int(
        1000, usedefault=True,
        argstr='--nvoxhp=%s',
        desc=('# of voxels used to estimate the '
              'hyperparameters'))

    fudge_factor = traits.Float(
        10.0, usedefault=True,
        argstr='--ff=%s',
        desc=('Fudge factor for hyperparameter '
              'error variance'))

    dont_sep_offs_move = traits.Bool(
        False,
        argstr='--dont_sep_offs_move',
        desc=('Do NOT attempt to separate '
              'field offset from subject '
              'movement'))

    dont_peas = traits.Bool(
        False,
        argstr='--dont_peas',
        desc="Do NOT perform a post-eddy alignment of "
        "shells")

    fwhm = traits.Float(
        desc=('FWHM for conditioning filter when estimating '
              'the parameters'),
        argstr='--fwhm=%s')

    niter = traits.Int(5, usedefault=True,
                       argstr='--niter=%s', desc='Number of iterations')

    method = traits.Enum(
        'jac',
        'lsr',
        argstr='--resamp=%s',
        desc=('Final resampling method (jacobian/least '
              'squares)'))

    repol = traits.Bool(
        False, argstr='--repol', desc='Detect and replace outlier slices')

    num_threads = traits.Int(
        1,
        usedefault=True,
        nohash=True,
        desc="Number of openmp threads to use")

    is_shelled = traits.Bool(
        False,
        argstr='--data_is_shelled',
        desc="Override internal check to ensure that "
        "date are acquired on a set of b-value "
        "shells")

    field = traits.Str(
        argstr='--field=%s',
        desc="NonTOPUP fieldmap scaled in Hz - filename has "
        "to be provided without an extension. TOPUP is "
        "strongly recommended")

    field_mat = File(
        exists=True,
        argstr='--field_mat=%s',
        desc="Matrix that specifies the relative locations of "
        "the field specified by --field and first volume "
        "in file --imain")

    use_cuda7_5 = traits.Bool(
        False,
        desc="Run eddy using CUDA 7.5 (for FSL < 6.0)",
    )

    use_cuda8_0 = traits.Bool(
        False,
        desc="Run eddy using CUDA 8.0 (for FSL >= 6.0)",
        xor=['use_cuda9_1']
    )

    use_cuda9_1 = traits.Bool(
        False,
        desc="Run eddy using CUDA 9.1 (for FSL >= 6.0)",
        xor=['use_cuda8_0']
    )

    initrand = traits.Int(
        False,
        argstr='--initrand=%s',
        desc='Initialize the random number generator'
    )


class EddyOutputSpec(TraitedSpec):
    out_corrected = File(
        exists=True, desc='4D image file containing all the corrected volumes')
    out_parameter = File(
        exists=True,
        desc=('text file with parameters definining the field and'
              'movement for each scan'))
    out_rotated_bvecs = File(
        exists=True, desc='File containing rotated b-values for all volumes')
    out_movement_rms = File(
        exists=True, desc='Summary of the "total movement" in each volume')
    out_restricted_movement_rms = File(
        exists=True,
        desc=('Summary of the "total movement" in each volume '
              'disregarding translation in the PE direction'))
    out_shell_alignment_parameters = File(
        exists=True,
        desc=('File containing rigid body movement parameters '
              'between the different shells as estimated by a '
              'post-hoc mutual information based registration'))
    out_outlier_report = File(
        exists=True,
        desc=('Text-file with a plain language report on what '
              'outlier slices eddy has found'))


class Eddy(FSLCommand):
    """
    Interface for FSL eddy, a tool for estimating and correcting eddy
    currents induced distortions. `User guide
    <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Eddy/UsersGuide>`_ and
    `more info regarding acqp file
    <http://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy/Faq#How_do_I_know_what_to_put_into_my_--acqp_file>`_.

    Examples
    --------

    >>> from nipype.interfaces.fsl import Eddy
    >>> eddy = Eddy()
    >>> eddy.inputs.in_file = 'epi.nii'
    >>> eddy.inputs.in_mask  = 'epi_mask.nii'
    >>> eddy.inputs.in_index = 'epi_index.txt'
    >>> eddy.inputs.in_acqp  = 'epi_acqp.txt'
    >>> eddy.inputs.in_bvec  = 'bvecs.scheme'
    >>> eddy.inputs.in_bval  = 'bvals.scheme'
    >>> eddy.inputs.use_cuda8_0 = True
    >>> eddy.cmdline # doctest: +ELLIPSIS
    'eddy_cuda8.0 --ff=10.0 --acqp=epi_acqp.txt --bvals=bvals.scheme \
--bvecs=bvecs.scheme --imain=epi.nii --index=epi_index.txt \
--mask=epi_mask.nii --niter=5 --nvoxhp=1000 --out=.../eddy_corrected'
    >>> eddy.inputs.use_cuda8_0 = False
    >>> eddy.cmdline # doctest: +ELLIPSIS
    'eddy_openmp --ff=10.0 --acqp=epi_acqp.txt --bvals=bvals.scheme \
--bvecs=bvecs.scheme --imain=epi.nii --index=epi_index.txt \
--mask=epi_mask.nii --niter=5 --nvoxhp=1000 --out=.../eddy_corrected'
    >>> res = eddy.run() # doctest: +SKIP

    """
    _cmd = 'eddy_openmp'
    input_spec = EddyInputSpec
    output_spec = EddyOutputSpec

    _num_threads = 1

    def __init__(self, **inputs):
        super(Eddy, self).__init__(**inputs)
        self.inputs.on_trait_change(self._num_threads_update, 'num_threads')
        if not isdefined(self.inputs.num_threads):
            self.inputs.num_threads = self._num_threads
        else:
            self._num_threads_update()
        self.inputs.on_trait_change(self._use_cuda7_5, 'use_cuda7_5')
        self.inputs.on_trait_change(self._use_cuda8_0, 'use_cuda8_0')
        self.inputs.on_trait_change(self._use_cuda9_1, 'use_cuda9_1')
        if isdefined(self.inputs.use_cuda7_5):
            self._use_cuda7_5()
        if isdefined(self.inputs.use_cuda8_0):
            self._use_cuda8_0()
        if isdefined(self.inputs.use_cuda9_1):
            self._use_cuda9_1()

    def _num_threads_update(self):
        self._num_threads = self.inputs.num_threads
        if not isdefined(self.inputs.num_threads):
            if 'OMP_NUM_THREADS' in self.inputs.environ:
                del self.inputs.environ['OMP_NUM_THREADS']
        else:
            self.inputs.environ['OMP_NUM_THREADS'] = str(
                self.inputs.num_threads)

    def _use_cuda7_5(self):
        self._cmd = 'eddy_cuda' if self.inputs.use_cuda7_5 else 'eddy_openmp'

    def _use_cuda8_0(self):
        self._cmd = 'eddy_cuda8.0' if self.inputs.use_cuda8_0 else 'eddy_openmp'

    def _use_cuda9_1(self):
        self._cmd = 'eddy_cuda9.1' if self.inputs.use_cuda9_1 else 'eddy_openmp'

    def _run_interface(self, runtime):
        # If 'eddy_openmp' is missing, use 'eddy'
        FSLDIR = os.getenv('FSLDIR', '')
        cmd = self._cmd
        if all((FSLDIR != '', cmd == 'eddy_openmp',
                not os.path.exists(os.path.join(FSLDIR, 'bin', cmd)))):
            self._cmd = 'eddy'
        runtime = super(Eddy, self)._run_interface(runtime)

        # Restore command to avoid side-effects
        self._cmd = cmd
        return runtime

    def _format_arg(self, name, spec, value):
        if name == 'in_topup_fieldcoef':
            return spec.argstr % value.split('_fieldcoef')[0]
        if name == 'out_base':
            return spec.argstr % os.path.abspath(value)
        return super(Eddy, self)._format_arg(name, spec, value)

    def _list_outputs(self):
        outputs = self.output_spec().get()
        outputs['out_corrected'] = os.path.abspath(
            '%s.nii.gz' % self.inputs.out_base)
        outputs['out_parameter'] = os.path.abspath(
            '%s.eddy_parameters' % self.inputs.out_base)

        # File generation might depend on the version of EDDY
        out_rotated_bvecs = os.path.abspath(
            '%s.eddy_rotated_bvecs' % self.inputs.out_base)
        out_movement_rms = os.path.abspath(
            '%s.eddy_movement_rms' % self.inputs.out_base)
        out_restricted_movement_rms = os.path.abspath(
            '%s.eddy_restricted_movement_rms' % self.inputs.out_base)
        out_shell_alignment_parameters = os.path.abspath(
            '%s.eddy_post_eddy_shell_alignment_parameters' %
            self.inputs.out_base)
        out_outlier_report = os.path.abspath(
            '%s.eddy_outlier_report' % self.inputs.out_base)

        if os.path.exists(out_rotated_bvecs):
            outputs['out_rotated_bvecs'] = out_rotated_bvecs
        if os.path.exists(out_movement_rms):
            outputs['out_movement_rms'] = out_movement_rms
        if os.path.exists(out_restricted_movement_rms):
            outputs['out_restricted_movement_rms'] = \
                out_restricted_movement_rms
        if os.path.exists(out_shell_alignment_parameters):
            outputs['out_shell_alignment_parameters'] = \
                out_shell_alignment_parameters
        if os.path.exists(out_outlier_report):
            outputs['out_outlier_report'] = out_outlier_report

        return outputs
