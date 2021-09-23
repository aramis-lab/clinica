import json
import functools
from os import read

from pathlib import Path
from typing import Optional


def provenance(func):
    from .provenance_utils import get_files_list

    @functools.wraps(func)
    def run_wrapper(self, **kwargs):
        ret = []
        pipeline_fullname = self.fullname
        in_files_paths = get_files_list(self, pipeline_fullname, dict_field="input_to")

        prov_context = get_context(files_paths=in_files_paths)
        prov_command = get_command(self, in_files_paths)

        if validate_command(prov_context, prov_command):
            ret = func(self)
        else:
            raise Exception(
                "The pipeline selected is incompatible with the input files provenance"
            )
        out_files_paths = get_files_list(
            self, pipeline_fullname, dict_field="output_from"
        )
        register_prov(prov_command, out_files_paths)

        return ret

    return run_wrapper


def register_prov(prov_command: dict, out_files: list) -> bool:

    # TODO: iterate over out_files and create a provenance file for each
    for file in out_files:
        write_prov_file(prov_command, file)
    print("Provenance registered succesfully")
    return True


def get_context(files_paths: str) -> dict:
    """
    Return a dictionary with the provenance info related to the files in the files_paths
    """
    from clinica.engine.provenance_utils import read_prov, get_associated_prov

    prov_data = {"Entity": [], "Agent": [], "Activity": []}
    for path in files_paths:
        prov_record = read_prov(get_associated_prov(path))
        if prov_record:
            prov_data = append_prov_dict(prov_data, prov_record)

    return prov_data


def get_command(self, input_files_paths: list) -> dict:
    """
    Read the user command and save information in a dict
    """
    import sys

    new_entities = []
    new_agent = get_agent()
    for path in input_files_paths:
        new_entities.append(get_entity(path))
    new_activity = get_activity(self, new_agent["@id"], new_entities)

    return {
        "Agent": [new_agent],
        "Activity": [new_activity],
        "Entity": new_entities,
    }


def write_prov_file(prov_command, files_paths):
    """
    Write the dictionary data to the file_path
    """
    from clinica.engine.provenance_utils import read_prov, get_associated_prov

    for file_path in files_paths:
        prov_path = get_associated_prov(file_path)

        if prov_path.exists():
            # append the pipeline provenance information to the old provenance file
            prov_main = read_prov(prov_path)
            prov_main = append_prov_dict(prov_main, prov_command)
        else:
            print("help")
            # create new provenance file with pipeline information
    return ""


def append_prov_dict(prov_main: dict, prov_new: dict) -> dict:
    """
    Append a specific prov data to the global prov dict
    """

    for k in prov_new.keys():
        for el in prov_new[k]:
            if prov_main[k] and el not in prov_main[k]:
                prov_main[k].append(el)
    return prov_main


def get_agent() -> dict:
    import clinica
    from .provenance_utils import get_agent_id

    agent_version = clinica.__version__
    agent_label = clinica.__name__
    agent_id = get_agent_id(agent_label + agent_version)

    new_agent = {"@id": agent_id, "label": agent_label, "version": agent_version}

    return new_agent


def get_activity(self, agent_id: str, entities: list) -> dict:
    """
    Add the current command to the list of activities
    """
    import sys
    from .provenance_utils import get_activity_id

    activity_parameters = self.parameters
    activity_label = self.fullname
    activity_id = get_activity_id(self.fullname)
    activity_command = (sys.argv[1:],)
    activity_agent = agent_id
    activity_used_files = [e["@id"] for e in entities]

    new_activity = {
        "@id": activity_id,
        "label": activity_label,
        "command": activity_command,
        "parameters": activity_parameters,
        "wasAssociatedWith": activity_agent,
        "used": activity_used_files,
    }

    return new_activity


def get_entity(img_path: str) -> dict:
    """
    Add the current file to the list of entities
    """
    from clinica.engine.provenance_utils import get_entity_id
    from clinica.engine.provenance_utils import get_last_activity
    from pathlib import Path

    entity_id = get_entity_id(img_path)
    entity_label = Path(img_path).name
    entity_path = img_path
    entity_source = get_last_activity(img_path)

    new_entity = {
        "@id": entity_id,
        "label": entity_label,
        "atLocation": entity_path,
        "wasGeneratedBy": entity_source,
    }

    return new_entity


def create_prov_file(command, path):
    """
    Create new provenance file based on command
    """
    # TODO: create a json-ld object next to the file and add it to the active prov object
    return


def validate_command(prov_context: dict, prov_command: dict) -> bool:
    """
    Check the command is valid on the data being run
    """
    flag = True
    new_activity_id = prov_command["Activity"][0]["@id"]
    new_agent_id = prov_command["Agent"][0]["@id"]

    for entity in prov_context["Entity"]:
        old_activity_id = entity["wasGeneratedBy"]
        if old_activity_id:
            ptr_activity = next(
                item
                for item in prov_context["Activity"]
                if item["@id"] == old_activity_id
            )
            old_agent_id = ptr_activity["wasAssociatedWith"]
            flag and is_valid(
                {(old_agent_id, old_activity_id): (new_agent_id, new_activity_id)}
            )
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
