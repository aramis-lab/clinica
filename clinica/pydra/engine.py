import functools
from os import PathLike
from typing import Callable

from pydra import Workflow
from pydra.engine.core import TaskBase

import clinica.pydra.engine_utils as pu
from clinica.pydra.interfaces import bids_reader, bids_writer, caps_reader


def clinica_io(func):
    """Decorator to add BIDS reader/writer to any Pydra workflow."""

    @functools.wraps(func)
    def run_wrapper(
        name: str, parameters: dict, input_dir: PathLike, output_dir: PathLike
    ) -> Workflow:

        core_workflow = func(parameters=parameters)

        pipeline = Workflow(
            name=name,
            input_spec=["bids_dir", "caps_dir"],
            bids_dir=input_dir,
            caps_dir=output_dir,
        )

        build_input_workflow(pipeline, core_workflow)

        # TODO: define condition on split if multiple fields
        pipeline.add(core_workflow)

        build_output_workflow(pipeline, core_workflow, output_dir)

        return pipeline

    return run_wrapper


def add_input_reading_task(
    pipeline: Workflow,
    core_workflow: Workflow,
    query_maker: Callable,
    reader: TaskBase,
) -> Workflow:
    """Configure and add the reading tasks of input workflow.

    Parameters
    ----------
    pipeline : Workflow
        The main Workflow to which the readers should be added.

    core_inputs : dict
        The inputs specified by the core workflow. This defines
        what the reader workflow should read and how it should
        connect to the core Workflow.

    query_maker : Callable
        Function responsible for parsing the core_inputs into a
        proper query.

    reader : TaskBase
        Task responsible for reading data.

    Returns
    -------
    pipeline : Workflow
        The main Workflow with readers added to it.
    """
    query = query_maker(pu.list_workflow_inputs(core_workflow))
    if len(query) == 0:
        return pipeline
    input_dir = "bids_dir" if "bids" in reader.__name__ else "caps_dir"
    input_workflow = Workflow(
        name=f"input_workflow_{reader.__name__}",
        input_spec=["input_dir"],
    )
    try:
        input_workflow.inputs.input_dir = getattr(pipeline.lzin, input_dir)
    except AttributeError:
        raise AttributeError(
            f"Workflow has no {input_dir} input. Please verify your input specifications."
        )
    input_workflow = add_input_task(
        input_workflow,
        reader(query=query, input_dir=input_workflow.lzin.input_dir),
    )
    pipeline.add(input_workflow)

    # connect input workflow to core workflow
    for field, _ in query.items():
        read_data = getattr(input_workflow.lzout, field)
        setattr(core_workflow.inputs, field, read_data)

    return pipeline


add_input_reading_task_bids = functools.partial(
    add_input_reading_task,
    query_maker=pu.bids_query,
    reader=bids_reader,
)


add_input_reading_task_caps = functools.partial(
    add_input_reading_task,
    query_maker=pu.caps_query,
    reader=caps_reader,
)


def build_input_workflow(pipeline: Workflow, core_workflow: Workflow) -> Workflow:
    """Setup for an input workflow.

    For now, the input workflow is responsible for:

        - reading BIDS data
        - reading CAPS data

    Parameters
    ----------
    pipeline :  Workflow
         The high level workflow containing (input -> core -> output).

    core_workflow : Workflow
         The functional workflow.

    Returns
    -------
    pipeline : Workflow
        The pipeline with the input workflow.
    """
    pipeline = add_input_reading_task_bids(pipeline, core_workflow)
    pipeline = add_input_reading_task_caps(pipeline, core_workflow)
    return pipeline


def add_input_task(input_workflow: Workflow, task: TaskBase) -> Workflow:
    """Add a task to the input workflow and define the workflow outputs.

    Parameters
    ----------
    input_workflow : Workflow
        The high level workflow containing (input -> core -> output).

    task : TaskBase
        The task to be added to the input_workflow.

    Returns
    -------
    input_workflow : Workflow
         The input workflow with the task added to it.
    """
    input_workflow.add(task)
    input_workflow.set_output(
        [
            (
                field,
                getattr(getattr(input_workflow, f"{task.name}").lzout, field),
            )
            for field in task.output_names
        ]
    )
    return input_workflow


def build_output_workflow(
    pipeline: Workflow, core_workflow: Workflow, output_dir: PathLike
) -> Workflow:
    """Example of an output workflow.

    Parameters
    ----------
    pipeline : Workflow
        the resulting workflow consisting of (input/core/output)
    core_workflow : Workflow
        contains the core interfaces
    output_dir : PathLike
        Path of the directory to be written to

    Returns
    -------
    Workflow
        The pipeline with the output workflow.
    """

    output_attrs = []

    for i, field in enumerate(pu.list_out_fields(core_workflow)):

        pipeline.add(bids_writer(name="bids_writer_task_" + str(field)))

        writer_task = getattr(pipeline, "bids_writer_task_" + str(field))
        writer_task_inputs = getattr(writer_task, "inputs")

        output_data = getattr(core_workflow.lzout, field)

        setattr(writer_task_inputs, "output_file", output_data)
        setattr(writer_task_inputs, "output_dir", output_dir)

        output_attrs.append((field, output_data))

    pipeline.set_output(output_attrs)

    return pipeline
