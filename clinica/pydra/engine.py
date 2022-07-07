import functools
from os import PathLike
from pathlib import PurePath

import clinica.pydra.engine_utils as pu
from pydra import Workflow
from clinica.pydra.interfaces import bids_reader, bids_writer

import typing as ty


def clinica_io(func):
    """
    Prepend an input workflow (BIDS reader)
    and append and output_workflow (BIDS writer)
    to any pydra (core) workflow
    """

    @functools.wraps(func)
    def run_wrapper(name: str, input_dir: PathLike, output_dir: PathLike, **kwargs):

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


def build_input_workflow(pipeline: Workflow, core_workflow: Workflow) -> ty.Tuple:
    """
    Setup for an input workflow to read BIDS data

    :pipeline: the high level workflow containing (input -> core -> output)
    :core_workflow: the functional workflow
    :return: the field to split on.
    """

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


def add_input_task(input_workflow: Workflow, query_dict) -> Workflow:
    """
    Construct and parameterize the input workflow

    :param pipeline: the high level workflow containing (input -> core -> output)
    :return: an BIDS reader based workflow
    """

    input_workflow.add(
        bids_reader(query_bids=query_dict, input_dir=input_workflow.lzin.input_dir)
    )

    data_keys = pu.list_keys(query_dict)

    input_workflow.set_output(
        [
            (field, getattr(input_workflow.bids_reader_task.lzout, field))
            for field in data_keys
        ]
    )
    return input_workflow


def build_output_workflow(
    pipeline: Workflow, core_workflow: Workflow, output_dir: PurePath
) -> None:
    """Example of an output workflow.

    :param name: The name of the workflow.
    :return: The output workflow.
    """

    output_attrs = []

    for i, field in enumerate(pu.list_out_fields(core_workflow)):

        pipeline.add(bids_writer(name="bids_writer_task_" + str(i)))

        writer_task = getattr(pipeline, "bids_writer_task_" + str(i))
        writer_task_inputs = getattr(writer_task, "inputs")

        output_data = getattr(core_workflow.lzout, field)

        setattr(writer_task_inputs, "output_file", output_data)
        setattr(writer_task_inputs, "output_dir", output_dir)

        output_attrs.append((field, output_data))

    pipeline.set_output(output_attrs)

    return pipeline
