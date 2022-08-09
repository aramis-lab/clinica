from os import PathLike, path
from pathlib import Path, PurePath

import pytest
from pydra import Workflow
from pydra.mark import annotate, task

CURRENT_DIR = path.dirname(path.realpath(__file__))


@task
@annotate({"return": {"smoothed_image": PurePath}})
def smooth_image(input_image: PathLike) -> PurePath:

    from nilearn.image import smooth_img

    smoothed_image = smooth_img(input_image, fwhm=3)

    return smoothed_image


@pytest.fixture
def bids_query():
    return {"T1w": {"suffix": "T1w", "extension": [".nii.gz"]}}


@pytest.fixture
def pipeline():
    return Workflow(
        name="pipeline_workflow",
        input_spec=["input_dir"],
        input_dir=Path(CURRENT_DIR) / "data",
    )


@pytest.fixture
def core_workflow():
    wf = Workflow("smoothing_t1w", input_spec=["T1w"])
    wf.add(
        smooth_image(
            name="smooth_image",
            interface=smooth_image,
            input_image=wf.lzin.T1w,
        )
    )
    wf.set_output([("smoothed_image", wf.smooth_image.lzout.smoothed_image)])
    return wf


def test_build_input_workflow(pipeline: Workflow, core_workflow: Workflow):
    from clinica.pydra.engine import build_input_workflow

    assert len(pipeline.nodes) == 0

    _ = build_input_workflow(pipeline, core_workflow)

    assert len(pipeline.nodes) == 1

    assert "input_workflow" in [x.name for x in pipeline.nodes]

    return


def test_add_input_task(bids_query: dict):
    from clinica.pydra.engine import add_input_task

    input_workflow = Workflow(name="input_workflow", input_spec=["input_dir"])

    input_workflow.inputs.input_dir = Path(CURRENT_DIR) / "data"

    add_input_task(input_workflow, bids_query)

    assert "bids_reader_task" in [x.name for x in input_workflow.nodes]
    return


def test_build_output_workflow(pipeline: Workflow, core_workflow: Workflow):
    from clinica.pydra.engine import build_output_workflow

    pipeline.add(core_workflow)

    decorated_pipeline = build_output_workflow(
        pipeline, core_workflow, Path(CURRENT_DIR) / "data"
    )
    assert "bids_writer_task_smoothed_image" in [
        x.name for x in decorated_pipeline.nodes
    ]
    return
