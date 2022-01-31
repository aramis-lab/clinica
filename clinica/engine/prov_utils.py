from typing import Optional, List
from pathlib import Path

from .prov_model import *


def get_files_list(self, pipeline_fullname: str, dict_field="input_to") -> List[Path]:
    """
    params:
        pipeline_fullname: the current running pipeline name
        dict_field: variable to specify if fetching inputs or outputs to the pipeline

    return list of 'Path's to the files used in the pipeline
    """
    from clinica.utils.inputs import clinica_file_reader
    import clinica.utils.input_files as cif

    dict_field_options = ["input_to", "output_from"]
    if dict_field not in dict_field_options:
        raise (f"dict_field must be one of {dict_field_options}")

    # Retrieve all the data dict from the input_files module

    files_dicts = {
        k: v
        for k, v in vars(cif).items()
        if isinstance(v, dict)
        and dict_field in v.keys()
        and pipeline_fullname in v[dict_field]
    }
    # TODO: check if bids or caps as output

    ret_files = []
    for elem in files_dicts:
        ref_dir = (
            self.bids_directory if dict_field == "input_to" else self.caps_directory
        )
        current_file, _ = clinica_file_reader(
            self.subjects,
            self.sessions,
            ref_dir,
            files_dicts[elem],
            raise_exception=False,
        )
        if current_file:
            ret_files.extend([Path(x) for x in current_file])

    return ret_files


def is_entity_tracked(prov_context: dict, entity_id: str) -> bool:
    flag_exists = next(
        (True for item in prov_context["Entity"] if item["@id"] == entity_id),
        False,
    )
    return flag_exists


def is_agent_tracked(prov_context: dict, agent_id: str) -> bool:
    flag_exists = next(
        (True for item in prov_context["Agent"] if item["@id"] == agent_id),
        False,
    )
    return flag_exists


def is_activity_tracked(prov_context: dict, activity_id: str) -> bool:
    flag_exists = next(
        (True for item in prov_context["Activity"] if item["@id"] == activity_id),
        False,
    )
    return flag_exists


def get_entity_id(path_file: Path) -> str:
    id = Identifier(label=path_file.with_suffix("").name)
    return id


def get_activity_id(pipeline_name: str) -> Identifier:
    id = Identifier(label="clin:" + pipeline_name)
    return id


def get_agent_id() -> Identifier:
    id = Identifier(label="RRID:Clinica")
    return id


def get_last_activity(path_entity: Path) -> Optional[ProvActivity]:

    """
    return the last activity executed on the file
    """

    prov_record = read_prov_jsonld(get_path_prov(path_entity))
    if prov_record and prov_record.entries:
        last_activity = prov_record.entries[-1]["@id"]
        return last_activity
    return None


def get_path_prov(path_entity: Path) -> Path:
    """
    return: Path of the provenance file associated with an entity
    """

    while path_entity.suffix != "":
        path_entity = path_entity.with_suffix("")

    path_prov = path_entity.with_suffix(".jsonld")
    return path_prov


def read_prov_jsonld(path_prov: Path) -> Optional[ProvRecord]:
    """
    return: ProvRecord in a specific location stored in jsonld format
    """

    if path_prov.exists():
        elements, prov_record = deserialize_jsonld(path_prov)
        return prov_record

    return None


def deserialize_jsonld(path_prov) -> ProvRecord:
    """
    params:

    return list of ProvEntry objects from jsonld dictionary data
    """

    import rdflib

    g = rdflib.Graph(identifier="prov_graph_records")
    g.parse(path_prov, format="json-ld")

    elements = {}
    entries = []

    # fetch context:
    context = ProvContext([])
    for lbl, link in g.namespace_manager.namespaces():
        namespace = Namespace(lbl, link.n3())
        context._namespaces.append(namespace)

    for s, p, o in g:

        if str(p) == "http://www.w3.org/ns/prov#Activity":
            id = Identifier(label=g.namespace_manager.qname(o))
            elements[id.label] = ProvActivity(id)

        elif str(p) == "http://www.w3.org/ns/prov#Agent":
            id = Identifier(label=g.namespace_manager.qname(o))
            elements[id.label] = ProvAgent(id)

        elif str(p) == "http://www.w3.org/ns/prov#Entity":
            id = Identifier(label=g.namespace_manager.qname(o))
            elements[id.label] = ProvEntity(id)

    for s, p, o in g:
        if type(s) != rdflib.term.BNode:
            attr = g.namespace_manager.qname(p).split(":")[1]

            subj = elements[g.namespace_manager.qname(s)]
            subj.attributes[attr] = str(o)

            curr_entry = ProvEntry(
                subject=g.namespace_manager.qname(s), predicate=attr, object=o
            )

            entries.append(curr_entry)

    prov_rec = ProvRecord(context=context, entries=entries)

    return elements, prov_rec
