# coding: utf8

"""
Utils for the implementation of t1-volume-tissue-segmentation pipeline
"""

import os
import numpy as np
from nipype.interfaces.spm.base import SPMCommandInputSpec, SPMCommand
from nipype.interfaces.base import TraitedSpec, OutputMultiPath, InputMultiPath, File, traits
from nipype.utils.filemanip import filename_to_list, list_to_filename

__author__ = "Jorge Samper-Gonzalez"
__copyright__ = "Copyright 2016-2019 The Aramis Lab Team"
__credits__ = ["Jorge Samper-Gonzalez"]
__license__ = "See LICENSE.txt file"
__version__ = "0.1.0"
__maintainer__ = "Jorge Samper-Gonzalez"
__email__ = "jorge.samper-gonzalez@inria.fr"
__status__ = "Development"


def select_bids_images(subjects, sessions, image_type, bids_layout):
    """
    read the subjects_visits_tsv
    """
    if len(subjects) != len(sessions):
        raise RuntimeError("Subjects list and sessions list must have the same length.")

    return [select_image(subjects[i], sessions[i], image_type, bids_layout) for i in range(len(subjects))]


def select_image(participant_id, session_id, image_type, bids_layout):
    """
    starting from the tsv with the participant_id and the session_id look for the corresponding image in the bids directory
    """
    import warnings

    if participant_id.startswith('sub-'):
        participant_id = participant_id[4:]
    if session_id.startswith('ses-'):
        session_id = session_id[4:]

    selected_images = bids_layout.get(subject=participant_id, session=session_id, type=image_type,
                                      return_type='file', extensions='nii.gz')
    if len(selected_images) == 0:
        selected_images = bids_layout.get(subject=participant_id, session=session_id, type=image_type,
                                          return_type='file', extensions='nii')
    if len(selected_images) == 0:
        raise RuntimeError('No ' + image_type + ' images were found for participant ' + participant_id
                           + ' and session ' + session_id)
    if len(selected_images) > 1:
        warnings.warn('Several ' + image_type + ' images were found for participant ' + participant_id
                      + ' and session ' + session_id, RuntimeWarning)
    return selected_images[0]


def group_nested_images_by_subject(class_images, zip_files=False):
    """

    """
    from clinica.utils.io import zip_nii

    if zip_files:
        return [zip_nii([s for tissue in subject for s in tissue], True) for subject in class_images]

    return [[s for tissue in subject for s in tissue] for subject in class_images]


def get_tissue_tuples(tissue_map, tissue_classes, dartel_tissues, save_warped_unmodulated, save_warped_modulated):
    """
    Method to obtain the list of tuples, one for each tissue class, with the following fields:
         - tissue probability map (4D), 1-based index to frame
         - number of gaussians
         - which maps to save [Native, DARTEL] - a tuple of two boolean values
         - which maps to save [Unmodulated, Modulated] - a tuple of two boolean values

    :param tissue_map: Path to tissue maps
    :param tissue_classes: Classes of images to obtain from segmentation. Ex: [1,2,3] is GM, WM and CSF
    :param dartel_tissues: Classes of images to save for DARTEL template calculation. Ex: [1] is only GM'
    :param save_warped_unmodulated: Save warped unmodulated images for tissues specified in --tissue_classes
    :param save_warped_modulated: Save warped modulated images for tissues specified in --tissue_classes
    :return: List of tuples according to NewSegment input por tissues
    """
    tissues = []

    for i in range(1, 7):
        n_gaussians = 2

        if i == 4 or i == 5:
            n_gaussians = i - 1

        native_space = False
        dartel_input = False
        warped_unmodulated = False
        warped_modulated = False

        if i in tissue_classes:
            native_space = True
            if save_warped_unmodulated:
                warped_unmodulated = True
            if save_warped_modulated:
                warped_modulated = True

        if i in dartel_tissues:
            dartel_input = True

        tissues.append(((tissue_map, i),
                        n_gaussians,
                        (native_space, dartel_input),
                        (warped_unmodulated, warped_modulated)))
    return tissues


class ApplySegmentationDeformationInput(SPMCommandInputSpec):

    deformation_field = File(
        exists=True,  mandatory=True,
        field='comp{1}.def',
        desc='SPM Segmentation deformation file')
    in_files = InputMultiPath(
        File(exists=True), mandatory=True,
        field='out{1}.pull.fnames',
        desc='Files on which deformation field is applied')
    interpolation = traits.Range(
        low=0, high=7,
        field='out{1}.pull.interp',
        desc='degree of b-spline used for interpolation')
    mask = traits.Int(
        0,
        usedefault=True,
        field='out{1}.pull.mask',
        desc='image masking')
    fwhm = traits.List(
        traits.Float(0),
        field='out{1}.pull.fwhm',
        minlen=3, maxlen=3,
        desc='3-element list (opt)')


class ApplySegmentationDeformationOutput(TraitedSpec):
    out_files = OutputMultiPath(File(exists=True),
                                desc='Transformed files')


class ApplySegmentationDeformation(SPMCommand):
    """ Uses spm to apply a deformation field obtained from Segmentation routine to a given file

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

    _jobtype = 'util'
    _jobname = 'defs'

    def _format_arg(self, opt, spec, val):
        """Convert input to appropriate format for spm
        """
        if opt == 'deformation_field':
            return np.array([list_to_filename(val)], dtype=object)
        if opt == 'in_files':
            return np.array(filename_to_list(val), dtype=object)
        return val

    def _list_outputs(self):
        outputs = self._outputs().get()
        outputs['out_files'] = []
        for filename in self.inputs.in_files:
            _, fname = os.path.split(filename)
            outputs['out_files'].append(os.path.realpath('w%s' % fname))
        return outputs
