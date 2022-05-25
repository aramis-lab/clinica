from pydra import Submitter, Workflow

from clinica.pipelines.interfaces import bids_reader_task, bids_writer_task


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
        self.dir_output = dir_output
        self.query_bids = query_bids

        self.add(bids_reader_task(self.dir_input, self.query_bids))
        self.list_outputs = []

        # optional_arg = kwargs.get('optional_arg')

    @property
    def input_task(self):
        return self.pydra_bids_data_reader

    @property
    def outputs(self):
        return self.list_outputs

    def set_output_task(self, inputs, output_dir):
        """
        Add a data writer for each file in the list
        TODO: replace with official pydra datasink if developped
        """

        for i, file in enumerate(inputs):

            self.list_outputs.append(file)

            self.add(
                bids_writer_task(
                    name="pydra_bids_data_writer" + str(i),
                    input_file=file,
                    output_dir=output_dir,
                )
            )
        return

    def get_nodes_names(self):
        """
        Return pipeline node names
        """
        return [x.name for x in self.nodes]

    def get_bids_files(self, modality):
        """
        Return input files of given modality
        """
        return self.input_task().get_output_field(modality)

    def run(self):
        """
        Submit the workflow
        """

        with Submitter(plugin="cf") as submitter:
            submitter(self)
        results = self.result(return_inputs=True)

        return results
