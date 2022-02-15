from pathlib import Path
from typing import List, Optional

from clinica.utils.input_files import pet_linear_nii

from .prov_model import *


def get_files_list(
    self, pipeline_fullname: str, dict_field="input_to", pipeline_args={}
) -> List[Path]:
    """
    params:
        pipeline_fullname: the current running pipeline name
        dict_field: variable to specify if fetching inputs or outputs to the pipeline

    return list of 'Path's to the files used in the pipeline
    """
    import clinica.utils.input_files as cif
    from clinica.utils.input_files import pet_linear_nii
    from clinica.utils.inputs import clinica_file_reader

    funcs = {"pet-linear": pet_linear_nii}

    dict_field_options = ["input_to", "output_from"]
    if dict_field not in dict_field_options:
        raise (f"dict_field must be one of {dict_field_options}")

    # Retrieve all the data dict from the input_files module

    if pipeline_fullname in funcs and dict_field == "output_from":
        files_dicts = {
            "PET": funcs[pipeline_fullname](
                **clean_arguments(pipeline_args, funcs[pipeline_fullname])
            )
        }
    else:
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
            raise_exception=True,
        )
        if current_file:
            ret_files.extend([Path(x) for x in current_file])

    return ret_files


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

    while path_entity.suffix != "":
        path_entity = path_entity.with_suffix("")

    path_prov = path_entity.with_suffix(".jsonld")
    return path_prov


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

    if path_prov.exists():
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
    g.parse(path_prov, format="json-ld")

    elements = {}

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
