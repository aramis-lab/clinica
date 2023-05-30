# from pathlib import PurePath
import abc
from os import PathLike
from typing import Dict, List

import nipype.interfaces.io as nio
import pydra
from nipype.interfaces.base import Directory, DynamicTraitedSpec, Str, isdefined, traits
from nipype.interfaces.io import IOBase, add_traits
from pydra.tasks.nipype1.utils import Nipype1Task

from clinica.pydra.query import BIDSQuery, CAPSFileQuery, CAPSGroupQuery, CAPSQuery


class CAPSDataGrabberInputSpec(DynamicTraitedSpec):
    base_dir = Directory(exists=True, desc="Path to CAPS Directory.", mandatory=True)
    output_query = traits.Dict(key_trait=Str, desc="Queries for outfield outputs")
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
        self.sessions = None
        self.subjects = None

    def _list_outputs(self):
        from clinica.utils.participant import get_subject_session_list

        sessions, subjects = get_subject_session_list(
            self.inputs.base_dir,
            is_bids_dir=False,
        )
        self.sessions = sessions
        self.subjects = subjects

        output_query = {}
        for k, query in self.inputs.output_query.items():
            if isinstance(query, list):
                temp = [self._execute_single_query(q) for q in query]
                if len(temp) != len(self.subjects) and len(temp[0]) == len(
                    self.subjects
                ):
                    temp = [list(i) for i in zip(*temp)]  # Transpose
                output_query[k] = temp
            else:
                output_query[k] = self._execute_single_query(query)
        return output_query

    @abc.abstractmethod
    def _execute_single_query(self, query: dict) -> list:
        pass

    def _add_output_traits(self, base):
        return add_traits(base, list(self.inputs.output_query.keys()))


class CAPSFileDataGrabber(CAPSDataGrabber):
    def _execute_single_query(self, query: Dict) -> List[str]:
        from clinica.utils.inputs import clinica_file_reader

        return clinica_file_reader(
            self.subjects,
            self.sessions,
            self.inputs.base_dir,
            query,
        )[0]


class CAPSGroupDataGrabber(CAPSDataGrabber):
    def _execute_single_query(self, query: Dict) -> List[str]:
        from clinica.utils.inputs import clinica_group_reader

        return clinica_group_reader(self.inputs.base_dir, query)


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


def bids_reader(query: BIDSQuery, input_dir: PathLike):
    """
    Parameters
    ----------
    query : BIDSQuery
        BIDS dataset reader query.
    input_dir :  PathLike
        The BIDS input directory.

    Returns
    -------
    Nipype1Task
        The task used for reading files from BIDS.
    """
    from pydra.tasks.bids import read_bids_dataset

    return read_bids_dataset(
        output_queries=dict(query),
        dataset_path=input_dir,
        name="bids_reader_task",
    )


def caps_reader(query: CAPSQuery, input_dir: PathLike):
    """
    Parameters
    ----------
    query : CAPSQuery
        Input to CAPSDataGrabber.
    input_dir :  PathLike
        The CAPS input directory.

    Returns
    -------
    Nipype1Task
        The task used for reading files from CAPS.
    """
    if isinstance(query, CAPSFileQuery):
        grabber = CAPSFileDataGrabber
        name = "caps_file_reader_task"
    elif isinstance(query, CAPSGroupQuery):
        grabber = CAPSGroupDataGrabber
        name = "caps_group_reader_task"
    else:
        raise TypeError(
            f"caps_reader received an unexpected query type {type(query)}. "
            "Supported types are: CAPSFileQuery and CAPSGroupQuery."
        )
    caps_data_grabber = grabber(output_query=query.query)
    caps_reader_task = Nipype1Task(
        name=name,
        interface=caps_data_grabber,
        base_dir=input_dir,
        output_query=query.query,
    )
    return caps_reader_task
