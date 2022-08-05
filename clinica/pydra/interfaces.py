# from pathlib import PurePath
from os import PathLike

import nipype.interfaces.io as nio
import pydra
from pydra.tasks.nipype1.utils import Nipype1Task


@pydra.mark.task
@pydra.mark.annotate({"return": {"output_file": str}})
def bids_writer(output_dir: PathLike, output_file: PathLike) -> str:
    """
    (toy) Interface to echo bids files

    Parameters
    ----------
    output_dir : PathLike
        output directory
    output_file : PathLike
        The output file to write
    Returns
    -------
    An BIDS reader based workflow
    """

    print("printing core vars in output_workflow", output_file)
    print("should go to location", output_dir)

    return output_file


def bids_reader(query_bids: dict, input_dir: PathLike):
    """
    Parameters
    ----------
    query_bids : dict
        Input to BIDSDataGrabber (c.f https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.io.html#bidsdatagrabber)
    input_dir :  PathLike
        the BIDS input directory

    Returns
    -------
        Nipype1Task
            The task used for reading files from BIDS
    """
    bids_data_grabber = nio.BIDSDataGrabber(output_query=query_bids)
    bids_reader_task = Nipype1Task(
        name="bids_reader_task",
        interface=bids_data_grabber,
        base_dir=input_dir,
    )
    return bids_reader_task
