import functools
from os import PathLike

from pydra import Workflow

import clinica.pydra.engine_utils as pu
from clinica.pydra.interfaces import bids_reader, bids_writer


def clinica_io(func):
    """Decorator to add BIDS reader/writer to any Pydra workflow."""

    @functools.wraps(func)
    def run_wrapper(name: str, input_dir: PathLike, output_dir: PathLike) -> Workflow:

        core_workflow = func()

        pipeline = Workflow(
            name=name,
            input_spec=["input_dir"],
            input_dir=input_dir,
        )

        split_key = build_input_workflow(pipeline, core_workflow)

        # TODO: define condition on split if multiple fields
        pipeline.add(core_workflow.split(split_key))

        build_output_workflow(pipeline, core_workflow, output_dir)

        return pipeline

    return run_wrapper


def build_input_workflow(pipeline: Workflow, core_workflow: Workflow) -> str:
    """Setup for an input workflow to read BIDS data.

    Parameters
    ----------
    pipeline :  Workflow
         the high level workflow containing (input -> core -> output)
    core_workflow : Workflow
         the functional workflow

    Returns
    ------
    str
        the field to split on.
    """

    field = ""

    list_core_inputs = pu.list_in_fields(core_workflow)
    query_dict = pu.bids_query(list_core_inputs)

    input_workflow = Workflow(name="input_workflow", input_spec=["input_dir"])

    input_workflow.inputs.input_dir = pipeline.lzin.input_dir

    input_workflow = add_input_task(input_workflow, query_dict)

    pipeline.add(input_workflow)

    # connect input workflow to core workflow

    for field in list_core_inputs:
        read_data = getattr(input_workflow.lzout, field)
        setattr(core_workflow.inputs, field, read_data)

    return field


def add_input_task(input_workflow: Workflow, query_bids: dict) -> Workflow:
    """Construct and parameterize the input workflow.

    Parameters
    ----------
    input_workflow : Workflow
        The high level workflow containing (input -> core -> output)
    query_bids : dict
        The dictionary containing the information needed to query the BIDS folder

    Returns
    -------
    Workflow
        An BIDS reader based workflow
    """

    input_workflow.add(
        bids_reader(query_bids=query_bids, input_dir=input_workflow.lzin.input_dir)
    )

    data_keys = list(query_bids.keys())

    input_workflow.set_output(
        [
            (field, getattr(input_workflow.bids_reader_task.lzout, field))
            for field in data_keys
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
        The output workflow.
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
