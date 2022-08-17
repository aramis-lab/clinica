import os
from os import PathLike
from typing import Union

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

from clinica.pydra.t1_volume.utils import get_world_coordinate_of_center, is_centered


@task
@annotate({"return": {"out_status": None}})
def check_volume_location_in_world_coordinate_system(
    nifti_input: Union[PathLike, list],
    bids_dir: PathLike,
    modality: str = "t1w",
    skip_question: bool = False,
) -> bool:
    """Check if images are centered around the origin of the world coordinate

    Parameters
    ----
    nifti_input: Union[list, PathLike]
        list of path to nifti files or path

    bids_dir: str
        path to bids directory associated with this check

    modality: str
        the modality of the image

    skip_question: bool
        if True, assume answer is yes

    Returns
    -------
    bool
        True if they are centered, False otherwise

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

    flag_centered = True
    list_non_centered_files = []

    if isinstance(nifti_input, list):
        list_non_centered_files = [
            file for file in nifti_input if not is_centered(file)
        ]
    elif os.path.exists(nifti_input):
        if not is_centered(nifti_input):
            list_non_centered_files.append(nifti_input)

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

        flag_centered = False

        if not skip_question:
            if not click.confirm("Do you still want to launch the pipeline?"):
                click.echo("Clinica will now exit...")
                sys.exit(0)

    return flag_centered


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
    """Uses SPM to apply a deformation field obtained from Segmentation routine to a given file"""

    input_spec = ApplySegmentationDeformationInput
    output_spec = ApplySegmentationDeformationOutput

    _jobtype = "util"
    _jobname = "defs"

    def _format_arg(self, opt, spec, val):
        """Convert input to appropriate format for spm"""

        import numpy as np

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
