from pathlib import PurePath
import nipype.interfaces.io as nio
from pydra import Submitter, Workflow

import clinica.pipelines.engine_pydra_utils as pu
from clinica.pipelines.interfaces import bids_writer, bids_reader


class Pipeline:
    def __init__(
        self, name: str, pipeline_in_dir: PurePath, pipeline_out_dir: PurePath
    ):
        """
        :pipeline_in_dir: a BIDS compliant input directory
        """
        self.pipeline_in_dir = pipeline_in_dir
        self.pipeline_out_dir = pipeline_out_dir
        self.name = name

    def build_workflow(self, core_workflow: Workflow):
        """
        Construct Clinica compliant pipeline
        :core_workflow: a pydra workflow with the core interfaces
        :query_bids: dict with keys which are a string and with a dict value containing (datatype/suffix/extension)
        """

        self.workflow = Workflow(
            name=self.name,
            input_spec=["input_dir"],
            input_dir=self.pipeline_in_dir,
        )

        # setup the input workflow
        list_core_inputs = pu.list_in_fields(core_workflow)

        self.input_workflow = pu.bids_query(list_core_inputs)
        self.input_workflow.inputs.input_dir = self.workflow.lzin.input_dir

        self.workflow.add(self.input_workflow)

        # connect input workflow to core workflow
        for field in list_core_inputs:
            bids_data = getattr(self.input_workflow.lzout, field)
            setattr(core_workflow.inputs, field, bids_data)

        # @TODO: define condition on split if multiple fields

        self.workflow.add(core_workflow.split(field))

        # add a bids_writer task for each core output
        for i, field in enumerate(pu.list_out_fields(core_workflow)):

            self.workflow.add(bids_writer(name="bids_writer_task_" + str(i)))
            output_data = getattr(core_workflow.lzout, field)
            writer_task = getattr(self.workflow, "bids_writer_task_" + str(i))
            inputs_writer_task = getattr(writer_task, "inputs")

            setattr(inputs_writer_task, "output_file", output_data)

        # @TODO: update with fields from output_workflow output_spec
        self.workflow.set_output(
            [("core_vars", getattr(core_workflow.lzout, "t1w_cropped_file"))]
        )

        return self.workflow

    @property
    def input_workflow(self) -> Workflow:
        return self._input_workflow

    @input_workflow.setter
    def input_workflow(self, query_bids: dict, name: str = "input"):

        """Generic Input workflow
        :name: The name of the workflow.
        :return: The input workflow.
        """
        self._input_workflow = Workflow(name=name, input_spec=["input_dir"])

        self._input_workflow.add(
            bids_reader(
                query_bids=query_bids, input_dir=self._input_workflow.lzin.input_dir
            )
        )

        data_keys = pu.list_keys(query_bids)

        self._input_workflow.set_output(
            [
                (field, getattr(self._input_workflow.bids_reader_task.lzout, field))
                for field in data_keys
            ]
        )

    def run(self):
        with Submitter(plugin="cf") as submitter:
            submitter(self.workflow)

        results = self.workflow.result(return_inputs=True)
        # @TODO: decide where to store results
        print(results)
