from clinica.pipelines.interfaces import bids_reader_task, bids_writer_task

from pydra import Submitter, Workflow
from attrs import define, field


class Pipeline(Workflow):
    """
    Wrapper around Pydra Workflow to mutualize input and output tasks among pipelines
    input_dir: the BIDS input directory of the pipeline
    output_dir: the BIDS output directory
    query_bids: a dict containing the details of the data to be queried in `input_dir`
    """

    def __init__(self, dir_input, query_bids, dir_output, **kwargs):

        super().__init__(**kwargs)

        self.dir_input = dir_input
        self.query_bids = query_bids
        self.dir_output = dir_output

        self.add(bids_reader_task(self.dir_input, self.query_bids))
        self.add(bids_writer_task(self.dir_output))

    @property
    def input_task(self):
        """
        Read BIDS-compliant datasets
        TODO:replace with custom interface
        """
        return self.pydra_bids_data_reader

    @property
    def output_task(self):
        """
        Generate BIDS-derivative complient output
        TODO:replace with custom interface
        """
        return self.pydra_bids_data_writer

    @property
    def outputs(self):
        """
        Return pipeline outputs as dict
        """
        return self.output_task._interface._outputs().get()

    def run(self):
        with Submitter(plugin="cf") as submitter:
            submitter(self)
        results = self.result(return_inputs=True)

        return results
