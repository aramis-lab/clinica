# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
# -*- coding: utf-8 -*-
from __future__ import (print_function, division, unicode_literals,
                        absolute_import)

import os.path as op

from nipype.interfaces.base import (CommandLineInputSpec, CommandLine, traits, TraitedSpec,
                                    File, InputMultiPath, isdefined)
from nipype.interfaces.mrtrix3.base import MRTrix3BaseInputSpec, MRTrix3Base


class TensorMetricsInputSpec(CommandLineInputSpec):
    in_file = File(
        exists=True,
        argstr='%s',
        mandatory=True,
        position=-1,
        desc='input DTI image')

    out_fa = File(argstr='-fa %s', desc='output FA file')
    out_adc = File(argstr='-adc %s', desc='output ADC file')
    out_ad = File(argstr='-ad %s', desc='output AD file')
    out_rd = File(argstr='-rd %s', desc='output RD file')
    out_evec = File(
        argstr='-vector %s', desc='output selected eigenvector(s) file')
    out_eval = File(
        argstr='-value %s', desc='output selected eigenvalue(s) file')
    component = traits.List(
        [1],
        usedefault=True,
        argstr='-num %s',
        sep=',',
        desc=('specify the desired eigenvalue/eigenvector(s). Note that '
              'several eigenvalues can be specified as a number sequence'))
    in_mask = File(
        exists=True,
        argstr='-mask %s',
        desc=('only perform computation within the specified binary'
              ' brain mask image'))
    modulate = traits.Enum(
        'FA',
        'none',
        'eval',
        argstr='-modulate %s',
        desc=('how to modulate the magnitude of the'
              ' eigenvectors'))


class TensorMetricsOutputSpec(TraitedSpec):
    out_fa = File(desc='output FA file')
    out_adc = File(desc='output ADC file')
    out_ad = File(desc='output AD file')
    out_rd = File(desc='output RD file')
    out_evec = File(desc='output selected eigenvector(s) file')
    out_eval = File(desc='output selected eigenvalue(s) file')


class TensorMetrics(CommandLine):
    """
    Compute metrics from tensors


    Example
    -------

    >>> import nipype.interfaces.mrtrix3 as mrt
    >>> comp = mrt.TensorMetrics()
    >>> comp.inputs.in_file = 'dti.mif'
    >>> comp.inputs.out_fa = 'fa.mif'
    >>> comp.cmdline                               # doctest: +ELLIPSIS
    'tensor2metric -num 1 -fa fa.mif dti.mif'
    >>> comp.run()                                 # doctest: +SKIP
    """

    _cmd = 'tensor2metric'
    input_spec = TensorMetricsInputSpec
    output_spec = TensorMetricsOutputSpec

    def _list_outputs(self):
        outputs = self.output_spec().get()

        for k in list(outputs.keys()):
            if isdefined(getattr(self.inputs, k)):
                outputs[k] = op.abspath(getattr(self.inputs, k))

        return outputs
