from pathlib import Path
from typing import List, Optional

from clinica.engine.prov_model import (
    Identifier,
    Namespace,
    ProvActivity,
    ProvAgent,
    ProvContext,
    ProvEntity,
    ProvRecord,
)


def mint_agent() -> ProvAgent:
    """
    return
        ProvAgent associated with running version of the software
    """
    from clinica import __name__, __version__
    from clinica.engine.prov_utils import generate_agent_id

    new_agent = ProvAgent(uid=generate_agent_id())

    new_agent.attributes["version"] = __version__
    new_agent.attributes["label"] = __name__

    return new_agent


def mint_activity(agent: Identifier, entities: List[ProvEntity]) -> ProvActivity:
    """
    return
        ProvActivity from related entities and associated agent
    """
    import sys

    from clinica.engine.prov_utils import generate_activity_id

    new_activity = ProvActivity(uid=generate_activity_id("testfullname"))

    new_activity.attributes["parameters"] = "testparameters"
    new_activity.attributes["label"] = "testfullname"
    new_activity.attributes["command"] = sys.argv[1:]
    new_activity.attributes["used"] = [str(x.uid) for x in entities]
    new_activity.attributes["wasAssociatedWith"] = str(agent.uid)

    return new_activity


def mint_entity(path_curr: Path) -> ProvEntity:
    """
    return an Entity object from the file in path_curr
    """

    from clinica.engine.prov_utils import generate_entity_id, get_last_activity

    new_entity = ProvEntity(uid=generate_entity_id(path_curr))
    new_entity.attributes["label"] = path_curr.name
    new_entity.attributes["path"] = str(path_curr)

    # TODO: implement function to return the latest associated activity
    new_entity.attributes["wasGeneratedBy"] = get_last_activity(path_curr)

    return new_entity


def generate_entity_id(path_file: Path) -> Identifier:
    id = Identifier(label=path_file.with_suffix("").name)
    return id


def generate_activity_id(pipeline_name: str) -> Identifier:
    id = Identifier(label="clin:" + pipeline_name)
    return id


def generate_agent_id() -> Identifier:
    id = Identifier(label="RRID:Clinica")
    return id


def get_last_activity(path_entity: Path) -> Optional[ProvActivity]:

    """
    return the last activity executed on the file
    """

    prov_record = read_prov_jsonld(get_path_prov(path_entity))
    if prov_record and prov_record.elements:
        # TODO: filter activities by date
        last_activity = [
            x for x in prov_record.elements if isinstance(x, ProvActivity)
        ][-1]
        return str(last_activity.uid)
    return None


def get_path_prov(path_entity: Path) -> Path:
    """
    return: Path of the provenance file associated with an entity
    """
    if path_entity.is_file():
        while path_entity.suffix != "":
            path_entity = path_entity.with_suffix("")
            path_prov = path_entity.with_suffix(".jsonld")
            return path_prov
    else:
        return None


def create_prov_file(prov_command, prov_path):
    """
    Create new provenance file based on command
    """
    import json

    with open(prov_path, "w") as fp:
        json.dump(prov_command.json(), fp, indent=4)

    return


def read_prov_jsonld(path_prov: Path) -> Optional[ProvRecord]:
    """
    return: ProvRecord in a specific location stored in jsonld format
    """

    if path_prov and path_prov.exists():
        prov_record = deserialize_jsonld(path_prov)
        return prov_record

    return None


def deserialize_jsonld(path_prov) -> ProvRecord:
    """
    params:

    return ProvRecord object from jsonld dictionary data
    """

    import rdflib

    g = rdflib.Graph(identifier="prov_graph_records")
    built_in_namepsaces = list(g.namespace_manager.namespaces())
    g.parse(path_prov, format="json-ld")
    json_namespaces = list(g.namespace_manager.namespaces())
    json_namespaces = list(set(json_namespaces) - set(built_in_namepsaces))

    elements = {}

    # fetch context:
    context = ProvContext([])

    for lbl, link in json_namespaces:
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

    prov_rec = ProvRecord(context=context, elements=list(elements.values()))

    return prov_rec


def clean_arguments(pipeline_args, file_func):
    import inspect

    argspec = inspect.getargspec(file_func)
    if not argspec.keywords:
        for key in pipeline_args.copy().keys():
            if key not in argspec.args:
                del pipeline_args[key]
    return pipeline_args


def validate_command(prov_history: ProvRecord, prov_current: ProvRecord) -> bool:
    """
    Check the command is valid on the data being run
    """
    flag = True

    for a in prov_history.elements:
        for b in prov_current.elements:
            # TODO: check that the record entries are compatible with the current entry
            flag = True
    return flag


def is_valid(command: dict) -> bool:
    valid_list = [
        {
            ("clin:clinica0.5.0", "clin:adni2Bids"): (
                "clin:clinica0.5.0",
                "clin:t1-linear",
            )
        }
    ]
    if command in valid_list:
        return True
    return False


def write_prov_file(
    list_prov_entries: ProvRecord, path_entity: Path, overwrite=False
) -> None:
    """
    Create provenance file with current pipeline information

    params:
    prov_entries: list of ProvEntry
    entity_path: path of the prov-associated element
    """

    prov_path = get_path_prov(path_entity)

    create_prov_file(list_prov_entries, prov_path)

    return
