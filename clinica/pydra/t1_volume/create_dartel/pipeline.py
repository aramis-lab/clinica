import nipype.interfaces.spm as spm
from nipype.algorithms.misc import Gunzip
from pydra import Workflow
from pydra.tasks.nipype1.utils import Nipype1Task

from clinica.pydra.engine import clinica_io


@clinica_io
def t1volume_create_dartel(
    name: str = "t1volume-create-dartel", parameters: dict = {}
) -> Workflow:
    """Workflow for tissue segmentation, bias correction and spatial normalization

    Parameters
    ----------
    name : str
        name of pipeline
    Returns
    -------
    Workflow
        pydra workflow for core functionalities
    """

    workflow = Workflow(
        name,
        input_spec=["T1w"],
    )

    workflow.add(
        Nipype1Task(name="task_unzip", interface=Gunzip(), in_file=workflow.lzin.T1w)
        .split("in_file")
        .combine("in_file")
    )

    dartel_template_task = Nipype1Task(
        name="dartel_template",
        interface=spm.DARTEL(),
    )

    dartel_template_task.inputs.image_files = workflow.task_unzip.lzout.out_file

    workflow.add(dartel_template_task)

    workflow.set_output(
        [
            ("test_out_template", workflow.dartel_template.lzout.final_template_file),
            ("test_out_unzipped", workflow.task_unzip.lzout.out_file),
        ]
    )

    return workflow
