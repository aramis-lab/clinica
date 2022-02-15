import functools
from os import read
from pathlib import Path
from typing import List

from clinica.engine.prov_utils import create_prov_file

from .prov_model import *


def provenance(func):
    from .prov_utils import get_files_list

    @functools.wraps(func)
    def run_wrapper(self, **kwargs):
        ret = []
        pipeline_args = self.parameters
        pipeline_fullname = self.fullname
        paths_input_files = get_files_list(
            self, pipeline_fullname, "input_to", pipeline_args
        )

        prov_history = get_history_record(paths_files=paths_input_files)
        prov_current = get_pipeline_record(self, paths_input_files)

        if validate_command(prov_history, prov_current):
            ret = func(self)
            print("The pipeline succesfully executed.")
        else:
            raise Exception(
                "The pipeline selected is incompatible with the input files provenance"
            )
        paths_out_files = get_files_list(
            self, pipeline_fullname, "output_from", pipeline_args
        )
        register_prov(prov_current, paths_out_files)

        return ret

    return run_wrapper


def register_prov(entries_current: ProvRecord, out_files: Path) -> None:

    # TODO: iterate over out_files and create a provenance file for each

    for file in out_files:
        write_prov_file(entries_current, file)
    print("Provenance registered succesfully")
    return True


def get_history_record(paths_files: List[Path]) -> ProvRecord:
    """
    return:
        a ProvRecord for the associated files in path_files
    """

    from .prov_utils import get_path_prov, read_prov_jsonld

    prov_record = ProvRecord({}, [])

    for path in paths_files:
        prov_record_tmp = read_prov_jsonld(get_path_prov(path))
        if prov_record_tmp:
            # TODO extend context as well
            prov_record.elements.extend(prov_record_tmp.elements)

    return prov_record


def get_pipeline_record(self, paths_inputs: List[Path]) -> ProvRecord:
    """
    params:
        paths_inputs: list of input entries paths
    return:
        ProvRecord associated with the launched pipeline
    """
    import sys

    elements = []
    new_agent = get_agent()
    elements.append(new_agent)
    new_entities = []

    for path in paths_inputs:
        entity_curr = get_entity(path)
        new_entities.append(entity_curr)
    elements.extend(new_entities)

    new_activity = get_activity(self, new_agent, new_entities)
    elements.append(new_activity)

    return ProvRecord(context={}, elements=elements)


def write_prov_file(
    list_prov_entries: ProvRecord, path_entity: Path, overwrite=False
) -> None:
    """
    Create provenance file with current pipeline information

    params:
    prov_entries: list of ProvEntry
    entity_path: path of the prov-associated element
    """

    from .prov_utils import get_path_prov

    prov_path = get_path_prov(path_entity)

    create_prov_file(list_prov_entries, prov_path)

    return


def extend_prov(prov_main: dict, prov_new: dict) -> dict:
    """
    Append a specific prov data to the global prov dict
    """

    for k in prov_new.keys():
        for el in prov_new[k]:
            if k in prov_main.keys() and el not in prov_main[k]:
                prov_main[k].append(el)
    return prov_main


def get_agent() -> ProvAgent:
    from clinica import __name__, __version__

    from .prov_utils import generate_agent_id

    new_agent = ProvAgent(uid=generate_agent_id())

    new_agent.attributes["version"] = __version__
    new_agent.attributes["label"] = __name__

    return new_agent


def get_activity(self, agent: Identifier, entities: List[ProvEntity]) -> ProvActivity:
    """
    return
        ProvActivity from related entities and associated agent
    """
    import sys

    from .prov_utils import generate_activity_id

    new_activity = ProvActivity(uid=generate_activity_id(self.fullname))

    new_activity.attributes["parameters"] = self.parameters
    new_activity.attributes["label"] = self.fullname
    new_activity.attributes["command"] = sys.argv[1:]
    new_activity.attributes["used"] = [str(x.uid) for x in entities]
    new_activity.attributes["wasAssociatedWith"] = str(agent.uid)

    return new_activity


def get_entity(path_curr: Path) -> ProvEntity:
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
