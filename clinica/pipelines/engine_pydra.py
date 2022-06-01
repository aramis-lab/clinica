from pydra import Submitter, Workflow

from clinica.pipelines.interfaces import bids_reader_task, bids_writer_task


class Pipeline:
    def __init__(self, pipeline_in_dir, pipeline_out_dir):
        """
        :pipeline_in_dir: a BIDS compliant input directory
        """
        self.pipeline_in_dir = pipeline_in_dir
        self.pipeline_out_dir = pipeline_out_dir

    def build_workflow(self, wf_core, query_bids):
        """
        Construct Clinica compliant pipeline
        :wf_core: a pydra workflow with the core interfaces
        :query_bids: dict with keys which are a string and with a dict value containing (datatype/suffix/extension)
        """

        self.workflow = Workflow(
            name="t1_linear",
            input_spec=["input_dir"],
            input_dir=self.pipeline_in_dir,
        )

        self.input_workflow = query_bids
        self.input_workflow.inputs.input_dir = self.workflow.lzin.input_dir
        self.workflow.add(self.input_workflow)

        for out_field in self.input_workflow.output_spec.fields:

            field_name = out_field[0]
            wf_core.inputs.t1w_file = getattr(self.input_workflow.lzout, field_name)

        # @TODO: define condition on split if multiple fields and make generic
        self.workflow.add(wf_core.split("t1w_file"))

        # @TODO: change to make generic
        self.output_workflow = "output"
        self.output_workflow.inputs.output_dir = self.pipeline_out_dir
        self.output_workflow.inputs.input_file = wf_core.lzout.t1w_cropped_file

        self.workflow.add(self.output_workflow)

        # @TODO: change output spec to make generic
        self.workflow.set_output(
            [
                ("t1w_cropped_file", self.output_workflow.lzout.output_file),
                ("xfm_file", wf_core.lzout.xfm_file),
            ]
        )

        return self.workflow

    @property
    def input_workflow(self):
        return self._input_workflow

    @property
    def output_workflow(self):
        return self._output_workflow

    @input_workflow.setter
    def input_workflow(self, output_query, name: str = "input") -> Workflow:

        """Generic Input workflow
        :param name: The name of the workflow.
        :return: The input workflow.
        """

        from pydra.tasks.nipype1.utils import Nipype1Task

        self._input_workflow = Workflow(name=name, input_spec=["input_dir"])

        bids_data_grabber = nio.BIDSDataGrabber(output_query=output_query)

        self._input_workflow.add(
            Nipype1Task(
                name="bids_data_grabber",
                interface=bids_data_grabber,
                base_dir=self._input_workflow.lzin.input_dir,
            )
        )

        self._input_workflow.set_output(
            [
                ("t1w_files", self._input_workflow.bids_data_grabber.lzout.T1w),
            ]
        )

    @output_workflow.setter
    def output_workflow(self, name: str = "output") -> Workflow:
        """Generic output workflow.

        :param name: The name of the workflow.
        :return: The output workflow.
        """
        import pydra

        @pydra.mark.task
        @pydra.mark.annotate({"return": {"output_file": str}})
        def bids_writer_task(input_file, output_dir):
            """
            (Dummy) Task to write files to output_dir
            """
            import subprocess

            output_file = f"{output_dir}/{input_file}"
            print(f"{output_file}")

            return output_file

        self._output_workflow = Workflow(
            name=name, input_spec=["input_file", "output_dir"]
        )

        self._output_workflow.add(
            bids_writer_task(
                name="bids_writer",
                input_file=self._output_workflow.lzin.input_file,
                output_dir=self._output_workflow.lzin.output_dir,
            )
        )

        self._output_workflow.set_output(
            [("output_file", self._output_workflow.bids_writer.lzout.output_file)]
        )

    def run(self):
        with Submitter(plugin="cf") as submitter:
            submitter(self.workflow)

        results = self.workflow.result(return_inputs=True)

        print(results)
