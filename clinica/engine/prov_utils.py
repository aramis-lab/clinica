from typing import Union, Optional
from pathlib import Path


def get_files_list(self, pipeline_fullname: str, dict_field="input_to") -> list:
    """
    Calls clinica_file_reader with the appropriate extentions
    """
    from clinica.utils.inputs import clinica_file_reader
    import clinica.utils.input_files as cif

    dict_field_options = ["input_to", "output_from"]
    if dict_field not in dict_field_options:
        raise (f"dict_field must be one of {dict_field_options}")

    # retrieve all the data dictionaries from the input_files module
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
            ret_files.extend(current_file)

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


def get_entity_id(file_path: str) -> str:
    from pathlib import Path

    entity_id = Path(file_path).with_suffix("").name
    return entity_id


def get_activity_id(pipeline_name: str) -> str:
    return "clin:" + pipeline_name


def get_agent_id(agent_name: str) -> str:
    return "clin:" + agent_name


def get_last_activity(file_path: str) -> Optional[list]:

    """
    Return the last activity executed on the file
    """

    prov_record = read_prov(get_associated_prov(file_path))
    if prov_record and prov_record["Activity"]:
        last_activity = prov_record["Activity"][-1]["@id"]
        return last_activity
    return None


def get_associated_prov(file_path: str) -> Path:

    file_path = Path(file_path)
    while file_path.suffix != "":
        file_path = file_path.with_suffix("")

    associated_jsonld = file_path.with_suffix(".jsonld")
    return associated_jsonld


def read_prov(prov_path: Path) -> Optional[dict]:
    """
    Check if the given file is a valid provenance json-ld
    """
    import json

    # TODO: check that the provenance file associations and uses exists
    if prov_path.exists():
        with open(prov_path, "r") as fp:
            json_ld_data = json.load(fp)
            return json_ld_data
    return None
