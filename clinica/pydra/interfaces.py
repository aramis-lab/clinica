from pathlib import PurePath

import nipype.interfaces.io as nio
import pydra
from pydra.tasks.nipype1.utils import Nipype1Task


@pydra.mark.task
@pydra.mark.annotate({"return": {"output_file": str}})
def bids_writer(output_dir, output_file):
    """
    (Dummy) Task to echo files generated by the core workflow
    """

    print("printing core vars in output_workflow", output_file)
    print("should go to location", output_dir)

    return "dummy_output_string"


def bids_reader(query_bids: dict, input_dir: PurePath):
    """
    :query_bids: input to BIDSDataGrabber (c.f https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.io.html#bidsdatagrabber)
    :input_dir: the BIDS input directory
    """
    bids_data_grabber = nio.BIDSDataGrabber(output_query=query_bids)
    bids_reader_task = Nipype1Task(
        name="bids_reader_task",
        interface=bids_data_grabber,
        base_dir=input_dir,
    )
    return bids_reader_task