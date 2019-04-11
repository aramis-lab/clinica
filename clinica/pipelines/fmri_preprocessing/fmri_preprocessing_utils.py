# coding: utf8

from __future__ import print_function, division, unicode_literals, \
    absolute_import
from nipype.utils.filemanip import (fname_presuffix, filename_to_list,
                                    list_to_filename)
from nipype.interfaces.base import (OutputMultiPath, TraitedSpec, isdefined,
                                    traits, File)
from nipype.interfaces.spm.base import (SPMCommand, scans_for_fname, func_is_3d,
                                        scans_for_fnames, SPMCommandInputSpec)
import nipype.pipeline.engine as npe
import nipype.interfaces.spm as spm
import nipype.interfaces.fsl as fsl

import os

if 'SPMSTANDALONE_HOME' in os.environ:
    if 'MCR_HOME' in os.environ:
        matlab_cmd = os.path.join(os.environ['SPMSTANDALONE_HOME'],
                                  'run_spm12.sh') \
                     + ' ' + os.environ['MCR_HOME'] \
                     + ' script'
        spm.SPMCommand.set_mlab_paths(matlab_cmd=matlab_cmd, use_mcr=True)


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
            (seg_node, add1_node, [('native_gm_image', 'in_file')]),
            (seg_node, add1_node, [('native_wm_image', 'operand_file')]),
            (seg_node, add2_node, [('native_csf_image', 'in_file')]),
            (add1_node, add2_node, [('out_file', 'operand_file')]),
            (add2_node, dil_node, [('out_file', 'in_file')]),
            (dil_node, ero_node, [('out_file', 'in_file')]),
            (ero_node, thre_node, [('out_file', 'in_file')]),
            (thre_node, fill_node, [('out_file', 'in_file')]),
            (fill_node, mask_node, [('out_file', 'mask_file')]),
        ])
