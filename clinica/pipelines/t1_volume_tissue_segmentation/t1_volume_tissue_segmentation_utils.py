"""Utils for the implementation of t1-volume-tissue-segmentation pipeline."""

import os
from os import PathLike

import numpy as np
from nipype.interfaces.base import (
    File,
    InputMultiPath,
    OutputMultiPath,
    TraitedSpec,
    traits,
)
from nipype.interfaces.spm.base import SPMCommand, SPMCommandInputSpec
from nipype.utils.filemanip import filename_to_list, list_to_filename


def init_input_node(t1w):
    """Extract "sub-<participant_id>_ses-<session_label>" from <t1w> and print begin message."""
    from clinica.utils.filemanip import get_subject_id
    from clinica.utils.ux import print_begin_image

    subject_id = get_subject_id(t1w)
    print_begin_image(subject_id)

    return subject_id, t1w


def print_end_pipeline(subject_id, final_file):
    """Display end message for <subject_id> when <final_file> is connected."""
    from clinica.utils.ux import print_end_image

    print_end_image(subject_id)


def zip_list_files(class_images, zip_files=False):
    from clinica.utils.filemanip import zip_nii

    if zip_files:
        return [zip_nii(tissue, True) for tissue in class_images]

    return [tissue for tissue in class_images]


def get_tissue_tuples(
    tissue_map: PathLike,
    tissue_classes: list,
    dartel_tissues: list,
    save_warped_unmodulated: bool,
    save_warped_modulated: bool,
) -> list:
    """Extract list of tuples, one for each tissue class.
    Parameters
    ---------
    tissue_map : PathLike
        Path to tissue maps
    tissue_classes: list
        Classes of images to obtain from segmentation. Ex: [1,2,3] is GM, WM and CSF
    dartel_tissues: list
        Classes of images to save for DARTEL template calculation. Ex: [1] is only GM'
    save_warped_unmodulated : bool
        Save warped unmodulated images for tissues specified in --tissue_classes
    save_warped_modulated: bool
        Save warped modulated images for tissues specified in --tissue_classes
    Returns
    -------
    list
        List of tuples according to NewSegment input for tissues
    Notes
    -----
     The returned list contains tissue probability map (4D), 1-based index to frame
     - number of gaussians
     - which maps to save [Native, DARTEL] - a tuple of two boolean values
     - which maps to save [Unmodulated, Modulated] - a tuple of two boolean values
    """
    tissues = []

    for i in range(1, 7):
        n_gaussians = 2

        if i == 4 or i == 5:
            n_gaussians = i - 1

        native_space = i in tissue_classes
        dartel_input = i in dartel_tissues
        warped_unmodulated = (i in tissue_classes) and save_warped_unmodulated
        warped_modulated = (i in tissue_classes) and save_warped_modulated

        tissues.append(
            (
                (tissue_map, i),
                n_gaussians,
                (native_space, dartel_input),
                (warped_unmodulated, warped_modulated),
            )
        )

    return tissues


class ApplySegmentationDeformationInput(SPMCommandInputSpec):

    deformation_field = File(
        exists=True,
        mandatory=True,
        field="comp{1}.def",
        desc="SPM Segmentation deformation file",
    )
    in_files = InputMultiPath(
        File(exists=True),
        mandatory=True,
        field="out{1}.pull.fnames",
        desc="Files on which deformation field is applied",
    )
    interpolation = traits.Range(
        low=0,
        high=7,
        field="out{1}.pull.interp",
        desc="degree of b-spline used for interpolation",
    )
    mask = traits.Int(
        0, usedefault=True, field="out{1}.pull.mask", desc="image masking"
    )
    fwhm = traits.List(
        traits.Float(0),
        field="out{1}.pull.fwhm",
        minlen=3,
        maxlen=3,
        desc="3-element list (opt)",
    )


class ApplySegmentationDeformationOutput(TraitedSpec):
    out_files = OutputMultiPath(File(exists=True), desc="Transformed files")


class ApplySegmentationDeformation(SPMCommand):
    """Uses SPM to apply a deformation field obtained from Segmentation routine to a given file

    Examples
    --------

    >>> import clinica.pipelines.t1_volume_tissue_segmentation.t1_volume_tissue_segmentation_utils as seg_utils
    >>> inv = seg_utils.ApplySegmentationDeformation()
    >>> inv.inputs.in_files = 'T1w.nii'
    >>> inv.inputs.deformation = 'y_T1w.nii'
    >>> inv.run() # doctest: +SKIP
    """

    input_spec = ApplySegmentationDeformationInput
    output_spec = ApplySegmentationDeformationOutput

    _jobtype = "util"
    _jobname = "defs"

    def _format_arg(self, opt, spec, val):
        """Convert input to appropriate format for spm"""
        if opt == "deformation_field":
            return np.array([list_to_filename(val)], dtype=object)
        if opt == "in_files":
            return np.array(filename_to_list(val), dtype=object)
        return val

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs["out_files"] = []
        for filename in self.inputs.in_files:
            _, fname = os.path.split(filename)
            outputs["out_files"].append(os.path.realpath("w%s" % fname))
        return outputs
