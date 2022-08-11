import os
from multiprocessing.dummy import Array
from os import PathLike
from pathlib import PosixPath
from typing import Union

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
from pydra.mark import annotate, task


@task
@annotate({"return": {"out_file": PathLike}})
def zip_nii(in_var: Union[PathLike, list], same_dir: bool = False) -> PathLike:
    """Zips a file or a list of files

    Parameters
    ---------
    in_var : Union[PathLike, list]
        file or list of files to zip

    same_dir : bool
        True if we want to zip in the origin directory, False if in cwd

    Returns
    ------
    PathLike
        path of the output

    """

    import gzip
    import os
    import shutil

    from nipype.utils.filemanip import split_filename
    from traits.trait_base import _Undefined

    if (in_var is None) or isinstance(in_var, _Undefined):
        return None

    if not isinstance(in_var, str):  # type(in_file) is list:
        return [zip_nii(f, same_dir) for f in in_var]

    orig_dir, base, ext = split_filename(str(in_var))

    # Already compressed

    if os.path.splitext("ext") == ".gz":
        return in_var

    # Not compressed

    out_file = os.abspath(
        os.join(orig_dir if same_dir else os.getcwd(), base + ext + ".gz")
    )

    with open(in_var, "rb") as f_in, gzip.open(out_file, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)

    return out_file


@task
@annotate({"return": {"out_status": None}})
def check_volume_location_in_world_coordinate_system(
    nifti_list: list,
    bids_dir: PathLike,
    modality: str = "t1w",
    skip_question: bool = False,
) -> None:
    """Check if images are centered around the origin of the world coordinate

    Parameters
    ----
    nifti_list: list
        list of path to nifti files

    bids_dir: str
        path to bids directory associated with this check

    modality: str
        the modality of the image

    skip_question: bool
        if True, assume answer is yes

    Returns
    -------
    None

    Warns
    ------
    If volume is not centered on origin of the world coordinate system

    Notes
    -----
    the NIfTI file list provided in argument are approximately centered around the origin of the
    world coordinates. Otherwise, issues may arise with further processing such as SPM segmentation. When not centered,
    we warn the user of the problem propose to exit clinica to run clinica iotools center-nifti or to continue with the execution
    of the pipeline

    """

    import sys
    from os.path import abspath, basename

    import click
    import numpy as np

    nifti_list = []
    if isinstance(nifti_list, PosixPath):
        if not is_centered(nifti_list):
            list_non_centered_files.append(nifti_list)
    elif isinstance(nifti_list, list):
        list_non_centered_files = [file for file in nifti_list if not is_centered(file)]

    if len(list_non_centered_files) > 0:
        centers = [
            get_world_coordinate_of_center(file) for file in list_non_centered_files
        ]
        l2_norm = [np.linalg.norm(center, ord=2) for center in centers]

        # File column width : 3 spaces more than the longest string to display
        file_width = 3 + max(len(basename(file)) for file in list_non_centered_files)

        # Center column width (with a fixed minimum size) : 3 spaces more than the longest string to display
        center_width = max(
            len("Coordinate of center") + 3,
            3 + max(len(str(center)) for center in centers),
        )

        warning_message = (
            f"It appears that {str(len(list_non_centered_files))} files "
            "have a center way out of the origin of the world coordinate system. SPM has a high "
            "probability to fail on these files (for co-registration or segmentation):\n\n"
        )
        warning_message += (
            "%-" + str(file_width) + "s%-" + str(center_width) + "s%-s"
        ) % ("File", "Coordinate of center", "Distance to origin")
        # 18 is the length of the string 'Distance to origin'
        warning_message += "\n" + "-" * (file_width + center_width + 18) + "\n"
        for file, center, l2 in zip(list_non_centered_files, centers, l2_norm):
            warning_message += (
                "%-" + str(file_width) + "s%-" + str(center_width) + "s%-25.2f\n"
            ) % (basename(file), str(center), l2)

        cmd_line = f"`clinica iotools center-nifti {abspath(bids_dir)} {abspath(bids_dir)}_centered --modality {modality}`"

        warning_message += (
            "\nIf you are trying to launch the t1-freesurfer pipeline, you can ignore this message "
            "if you do not want to run the pet-surface pipeline afterward."
        )

        warning_message += (
            "\nClinica provides a tool to counter this problem by replacing the center of the volume"
            " at the origin of the world coordinates.\nUse the following command line to correct the "
            f"header of the faulty NIFTI volumes in a new folder:\n{cmd_line}"
            "You will find more information on the command by typing "
            "clinica iotools center-nifti in the console."
        )

        click.echo(warning_message)

        if not skip_question:
            if not click.confirm("Do you still want to launch the pipeline?"):
                click.echo("Clinica will now exit...")
                sys.exit(0)
    return


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
