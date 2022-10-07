# from pathlib import PurePath
from os import PathLike

import nipype.interfaces.io as nio
import pydra
from nipype.interfaces.base import Directory, DynamicTraitedSpec, Str, isdefined, traits
from nipype.interfaces.io import IOBase, add_traits
from pydra.tasks.nipype1.utils import Nipype1Task

from clinica.pydra.query import BIDSQuery, CAPSFileQuery, CAPSGroupQuery, CAPSQuery


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


class CAPSFileDataGrabber(IOBase):
    input_spec = CAPSDataGrabberInputSpec
    output_spec = DynamicTraitedSpec

    def __init__(self, **kwargs):
        super(CAPSFileDataGrabber, self).__init__(**kwargs)

        if not isdefined(self.inputs.output_query):
            self.inputs.output_query = {}
        # used for mandatory inputs check
        undefined_traits = {}
        self.inputs.trait_set(trait_change_notify=False, **undefined_traits)

    def _list_outputs(self):
        from clinica.utils.inputs import clinica_file_reader
        from clinica.utils.participant import get_subject_session_list

        sessions, subjects = get_subject_session_list(
            self.inputs.base_dir,
            is_bids_dir=False,
        )
        output_query = {}
        for k, query in self.inputs.output_query.items():
            if isinstance(query, list):
                temp = [
                    clinica_file_reader(
                        subjects,
                        sessions,
                        self.inputs.base_dir,
                        q,
                    )[0]
                    for q in query
                ]
                if len(temp) != len(subjects) and len(temp[0]) == len(subjects):
                    transpose = []
                    for x in zip(*temp):
                        transpose.append(x)
                    assert len(transpose) == len(subjects)
                    temp = transpose
                output_query[k] = temp
            else:
                output_query[k] = clinica_file_reader(
                    subjects,
                    sessions,
                    self.inputs.base_dir,
                    query,
                )[0]
        return output_query

    def _add_output_traits(self, base):
        return add_traits(base, list(self.inputs.output_query.keys()))


class CAPSGroupDataGrabber(IOBase):
    input_spec = CAPSDataGrabberInputSpec
    output_spec = DynamicTraitedSpec

    def __init__(self, **kwargs):
        super(CAPSGroupDataGrabber, self).__init__(**kwargs)

        if not isdefined(self.inputs.output_query):
            self.inputs.output_query = {}
        # used for mandatory inputs check
        undefined_traits = {}
        self.inputs.trait_set(trait_change_notify=False, **undefined_traits)

    def _list_outputs(self):
        from clinica.utils.inputs import clinica_group_reader
        from clinica.utils.participant import get_subject_session_list

        sessions, subjects = get_subject_session_list(
            self.inputs.base_dir,
            is_bids_dir=False,
        )
        output_query = {}
        for k, query in self.inputs.output_query.items():
            if isinstance(query, list):
                output_query[k] = [
                    clinica_group_reader(self.inputs.base_dir, sub_query)
                    for sub_query in query
                ]
            else:
                output_query[k] = clinica_group_reader(self.inputs.base_dir, query)
        return output_query

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


def bids_reader(query: BIDSQuery, input_dir: PathLike):
    """
    Parameters
    ----------
    query : BIDSQuery
        Input to BIDSDataGrabber (c.f https://nipype.readthedocs.io/en/latest/api/generated/nipype.interfaces.io.html#bidsdatagrabber)
    input_dir :  PathLike
        The BIDS input directory.

    Returns
    -------
    Nipype1Task
        The task used for reading files from BIDS.
    """
    bids_data_grabber = nio.BIDSDataGrabber(output_query=query.query)
    bids_reader_task = Nipype1Task(
        name="bids_reader_task",
        interface=bids_data_grabber,
        base_dir=input_dir,
    )
    return bids_reader_task


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
    elif isinstance(query, CAPSGroupQuery):
        grabber = CAPSGroupDataGrabber
    else:
        raise TypeError(
            f"caps_reader received an unexpected query type {type(query)}. Supported types are: CAPSFileQuery and CAPSGroupQuery."
        )
    caps_data_grabber = grabber(output_query=query.query)
    caps_reader_task = Nipype1Task(
        name="caps_reader_task",
        interface=caps_data_grabber,
        base_dir=input_dir,
    )
    return caps_reader_task
