from pathlib import Path

import nipype.interfaces.io as nio
import nipype.interfaces.spm as spm
import nipype.interfaces.utility as nutil
import nipype.pipeline.engine as npe
from nipype.algorithms.misc import Gunzip
from pydra import Submitter, Workflow
from pydra.tasks.nipype1.utils import Nipype1Task

from clinica.pydra.engine import clinica_io
from clinica.pydra.t1_volume.t1_volume_utils import (
    ApplySegmentationDeformation,
    initialize_parameters,
)

in_dir = Path("/Users/omar.elrifai/workspace/experimentations/pydra/IN/")
nifti_list = [in_dir / "t1w.nii.gz"]
test_nii = in_dir / "t1w.nii"


@clinica_io
def t1volume_tissue_segmentation(name: str = "t1volume") -> Workflow:
    """Workflow for tissue segmentation, bias correction and spatial normalization"""

    from clinica.pydra.t1_volume.t1_volume_tasks import (
        check_volume_location_in_world_coordinate_system,
    )

    workflow = Workflow(
        name,
        input_spec=["T1w"],
    )

    workflow.add(
        check_volume_location_in_world_coordinate_system(
            name="check_location_world",
            nifti_list=workflow.lzin.T1w,
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

    tissue_tuples = initialize_parameters()
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
