from pydra.tasks.nipype1.utils import Nipype1Task
import nipype.interfaces.io as nio
import pydra


def bids_reader_task(input_dir, query_bids):
    """
    Return pydra Nipype1Task for querying bids dataset
    """

    bids_data_grabber = nio.BIDSDataGrabber(
            output_query=query_bids,
    )
    
    task_bids_data_grabber = Nipype1Task(
        name="pydra_bids_data_reader", base_dir=input_dir, interface=bids_data_grabber
    )
    
    return task_bids_data_grabber




@pydra.mark.task

def bids_writer_task(input_file, output_dir):
    """
    Task to write files to output_dir
    """
    import subprocess
    subprocess.run(["cp", input_file, output_dir])

    