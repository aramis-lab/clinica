import typing as ty

import nipype.interfaces.spm as spm
import pydra
from nipype.algorithms.misc import Gunzip
from pydra import Workflow
from pydra.tasks.nipype1.utils import Nipype1Task

from clinica.pipelines.t1_volume_tissue_segmentation.t1_volume_tissue_segmentation_utils import (
    ApplySegmentationDeformation,
)
from clinica.pydra.engine import clinica_io
from clinica.pydra.t1_volume.utils import initialize_tissues_spm_segment


@clinica_io
def t1volume_tissue_segmentation(
    name: str = "t1-volume-tissue-segmentation", parameters: dict = {}
) -> Workflow:
    """Core workflow for tissue segmentation, bias correction and spatial normalization

    Parameters
    ----------
    name : str
        name of pipeline. Default="t1-volume-tissue-segmentation".

    parameters : dict
        dictionary of pipeline parameters

    Returns
    -------
    Workflow
        pydra workflow for core functionalities
    """

    from clinica.pydra.t1_volume.tasks import (
        task_volume_location_in_world_coordinate_system,
    )

    input_spec = pydra.specs.SpecInfo(
        name="Input",
        fields=[("T1w", str, {"mandatory": True})],
        bases=(pydra.specs.BaseSpec,),
    )

    workflow = Workflow(
        name,
        input_spec=input_spec,
    )

    workflow.split("T1w")

    workflow.add(
        task_volume_location_in_world_coordinate_system(
            name="check_location_world",
            nifti_input=workflow.lzin.T1w,
            skip_question=parameters["skip_question"],
        )
    )

    workflow.add(
        Nipype1Task(name="unzipT1w", interface=Gunzip(), in_file=workflow.lzin.T1w)
    )

    spm_segment = Nipype1Task(
        name="SpmSegmentation",
        interface=spm.NewSegment(),
    )

    spm_segment.inputs.write_deformation_fields = [True, True]

    tissue_tuples = initialize_tissues_spm_segment(parameters)
    spm_segment.inputs.tissues = tissue_tuples

    spm_segment.inputs.channel_files = workflow.unzipT1w.lzout.out_file

    workflow.add(spm_segment)

    t1_to_mni = Nipype1Task(
        name="T1wToMNI",
        interface=ApplySegmentationDeformation(),
    )

    t1_to_mni.inputs.in_files = workflow.unzipT1w.lzout.out_file
    t1_to_mni.inputs.deformation_field = (
        workflow.SpmSegmentation.lzout.forward_deformation_field
    )

    workflow.add(t1_to_mni)

    workflow.set_output([("test_out", workflow.T1wToMNI.lzout.out_files)])

    return workflow
