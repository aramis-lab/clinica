from pydra.tasks.nipype1.utils import Nipype1Task
import nipype.interfaces.io as nio
from nipype import DataSink


def bids_reader_task(input_dir, query_bids):
    """
    Return pydra Nipype1Task for querying bids dataset
    """

    bids_data_grabber = nio.BIDSDataGrabber(
        base_dir=input_dir,
        output_query=query_bids,
    )

    task_bids_data_grabber = Nipype1Task(
        name="pydra_bids_data_reader", interface=bids_data_grabber
    )

    return task_bids_data_grabber


def bids_writer_task(output_dir, parametrization=False):
    """
    Return pydra Task for writting bids-compliant dataset
    """

    data_writer = DataSink(base_directory=output_dir, parametrization=parametrization)
    data_writer.inputs.structural = (
        "/Users/omar.elrifai/workspace/experimentations/pydra/IN/structural.nii.gz"
    )
    task_bids_data_writer = Nipype1Task(
        name="pydra_bids_data_writer", interface=data_writer
    )

    return task_bids_data_writer
