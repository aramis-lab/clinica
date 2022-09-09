# from pathlib import PurePath
from os import PathLike

import nipype.interfaces.io as nio
import pydra
from nipype.interfaces.base import Directory, DynamicTraitedSpec, Str, isdefined, traits
from nipype.interfaces.io import IOBase, add_traits
from pydra.tasks.nipype1.utils import Nipype1Task


class CAPSDataGrabberInputSpec(DynamicTraitedSpec):
    base_dir = Directory(exists=True, desc="Path to CAPS Directory.", mandatory=True)
    output_query = traits.Dict(
        key_trait=Str, value_trait=traits.Dict, desc="Queries for outfield outputs"
    )
    raise_on_empty = traits.Bool(
        True,
        usedefault=True,
        desc="Generate exception if list is empty for a given field",
    )


class CAPSDataGrabber(IOBase):
    input_spec = CAPSDataGrabberInputSpec
    output_spec = DynamicTraitedSpec

    def __init__(self, **kwargs):
        super(CAPSDataGrabber, self).__init__(**kwargs)

        if not isdefined(self.inputs.output_query):
            self.inputs.output_query = {}
        # used for mandatory inputs check
        undefined_traits = {}
        self.inputs.trait_set(trait_change_notify=False, **undefined_traits)

    def _list_outputs(self):
        from clinica.utils.inputs import clinica_file_reader, clinica_group_reader
        from clinica.utils.participant import get_subject_session_list

        sessions, subjects = get_subject_session_list(
            self.inputs.base_dir,
            is_bids_dir=False,
        )
        query = {}
        for k, q in self.inputs.output_query.items():
            reader = q.pop("reader")
            if reader == "file":
                query[k] = clinica_file_reader(
                    subjects,
                    sessions,
                    self.inputs.base_dir,
                    q,
                )
            elif reader == "group":
                query[k] = clinica_group_reader(self.inputs.base_dir, q)
            else:
                raise ValueError(f"Unknown reader {reader}.")
        return query

    def _add_output_traits(self, base):
        return add_traits(base, list(self.inputs.output_query.keys()))


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


def bids_reader(query: dict, input_dir: PathLike):
    """
    Parameters
    ----------
    query : dict
        Input to BIDSDataGrabber (c.f https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.io.html#bidsdatagrabber)
    input_dir :  PathLike
        The BIDS input directory.

    Returns
    -------
    Nipype1Task
        The task used for reading files from BIDS.
    """
    bids_data_grabber = nio.BIDSDataGrabber(output_query=query)
    bids_reader_task = Nipype1Task(
        name="bids_reader_task",
        interface=bids_data_grabber,
        base_dir=input_dir,
    )
    return bids_reader_task


def caps_reader(query: dict, input_dir: PathLike):
    """
    Parameters
    ----------
    query : dict
        Input to CAPSDataGrabber.
    input_dir :  PathLike
        The CAPS input directory.

    Returns
    -------
    Nipype1Task
        The task used for reading files from CAPS.
    """
    caps_data_grabber = CAPSDataGrabber(output_query=query)
    caps_reader_task = Nipype1Task(
        name="caps_reader_task",
        interface=caps_data_grabber,
        base_dir=input_dir,
    )
    return caps_reader_task
