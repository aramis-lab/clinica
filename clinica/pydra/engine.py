import functools
from os import PathLike
from typing import Callable

from pydra import Workflow

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
            input_spec=["input_dir"],
            input_dir=input_dir,
        )

        build_input_workflow(pipeline, core_workflow)

        # TODO: define condition on split if multiple fields
        pipeline.add(core_workflow)

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
    -------
    Workflow
        The pipeline with the input workflow.
    """

    field = ""

    list_core_inputs = pu.list_in_fields(core_workflow)
    bids_query_dict = pu.bids_query(list_core_inputs)

    caps_dict_inputs = pu.list_dict_in_fields(core_workflow)
    caps_query_dict = pu.caps_query(caps_dict_inputs)

    input_workflow = Workflow(name="input_workflow", input_spec=["input_dir"])

    input_workflow.inputs.input_dir = pipeline.lzin.input_dir

    for reader, query in zip(
        [bids_reader, caps_reader], [bids_query_dict, caps_query_dict]
    ):
        input_workflow = add_input_task(input_workflow, reader, query)
        pipeline.add(input_workflow)

    # connect input workflow to core workflow

    for field in list_core_inputs:
        read_data = getattr(input_workflow.lzout, field)
        setattr(core_workflow.inputs, field, read_data)

    return pipeline


def add_input_task(
    input_workflow: Workflow,
    reader: Callable,
    query: dict,
) -> Workflow:
    """Construct and parameterize the input workflow with a BIDS/CAPS query.

    Parameters
    ----------
    input_workflow : Workflow
        The high level workflow containing (input -> core -> output).
    reader : The reader to use (bids_reader or caps_reader)
    query : dict
        The dictionary containing the information needed to query the BIDS/CAPS folder.

    Returns
    -------
    Workflow
        An BIDS/CAPS reader based workflow.
    """

    input_workflow.add(reader(query=query, input_dir=input_workflow.lzin.input_dir))

    data_keys = list(query.keys())

    input_workflow.set_output(
        [
            (
                field,
                getattr(
                    getattr(input_workflow, f"{reader.__name__}_task").lzout, field
                ),
            )
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
