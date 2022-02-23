import functools

from os import read
from pathlib import Path
from typing import List


def provenance(func):
    @functools.wraps(func)
    def run_wrapper(self, **kwargs):
        ret = func(self)

        pipeline_args = self.parameters
        pipeline_fullname = self.fullname

        create_node_read(self)
        create_node_update(self, pipeline_args, pipeline_fullname)
        create_node_log(self)

        connect_nodes(self)

        return ret

    return run_wrapper


def connect_nodes(self):
    # fmt: off

    #self.output_node.outputs.get()[self.get_output_fields()[0]]

    self.connect(
        [
            (self.input_node, self.prov_input_node, [("t1w", "input_files")]),
            (self.input_node, self.prov_update_node, [("t1w", "input_files")]),
            (self.prov_input_node, self.prov_update_node, [("prov_in_record", "prov_in_record")]),
            (self.prov_update_node, self.prov_log_node,[("prov_upd_record", "prov_log_record")]),
            (self.output_node, self.prov_log_node, [(self.get_output_fields()[0], "out_file")]),        
        ]
    )
    return True
    # fmt: on


def create_node_read(self):
    import nipype.pipeline.engine as npe
    import nipype.interfaces.utility as nutil

    self.prov_input_node = npe.Node(
        nutil.Function(
            input_names=["input_files"],
            output_names=["prov_in_record"],
            function=read_prov,
        ),
        name="ReadProvRecord",
    )


def create_node_update(self, parameters, fullname):
    import nipype.pipeline.engine as npe
    import nipype.interfaces.utility as nutil

    self.prov_update_node = npe.Node(
        nutil.Function(
            input_names=["input_files", "prov_in_record", "parameters", "fullname"],
            output_names=["prov_upd_record"],
            function=update_prov,
        ),
        name="UpdateRecord",
    )

    return True


def create_node_log(self):
    import nipype.pipeline.engine as npe
    import nipype.interfaces.utility as nutil

    self.prov_log_node = npe.Node(
        nutil.Function(
            input_names=["prov_log_record", "out_file", "out_dir"],
            output_names=["output_record"],
            function=log_prov,
        ),
        name="LogProv",
    )

    self.prov_log_node.inputs.out_dir = self.caps_directory
    return


def read_prov(input_files):
    """
    return:
        a ProvRecord for the associated files in path_files
    """
    from clinica.engine.prov_utils import get_path_prov, read_prov_jsonld
    from clinica.engine.prov_model import ProvRecord
    from pathlib import Path

    prov_record = ProvRecord({}, [])
    if isinstance(input_files, list):
        paths_files = [Path(x) for x in input_files]
    elif isinstance(input_files, str):
        paths_files = [Path(input_files)]

    for path in paths_files:
        print("in read_prov, path for input:", path)
        prov_record_tmp = read_prov_jsonld(get_path_prov(path))
        if prov_record_tmp:
            # TODO extend context as well
            prov_record.elements.extend(prov_record_tmp.elements)

    return prov_record


def update_prov(input_files, prov_in_record):
    """
    params:
        input_files: list of input entries
    return:
        ProvRecord associated with the launched pipeline
    """
    from clinica.engine.prov_utils import (
        mint_activity,
        mint_agent,
        mint_entity,
        validate_command,
    )
    from pathlib import Path
    from clinica.engine.prov_model import ProvRecord

    elements = []
    new_agent = mint_agent()
    elements.append(new_agent)
    new_entities = []

    if isinstance(input_files, list):
        paths_files = [Path(x) for x in input_files]
    elif isinstance(input_files, str):
        paths_files = [Path(input_files)]

    for path in paths_files:
        entity_curr = mint_entity(path)
        new_entities.append(entity_curr)
    elements.extend(new_entities)

    new_activity = mint_activity(new_agent, new_entities)
    elements.append(new_activity)

    prov_current = ProvRecord(context={}, elements=elements)

    if not validate_command(prov_in_record, prov_current):
        raise ("Invalid commmand")
    return prov_current


def log_prov(prov_log_record, out_file, out_dir):
    from clinica.engine.prov_utils import write_prov_file
    from pathlib import Path

    out_file = out_file + "*"
    out_files_paths = []
    if isinstance(out_file, list):
        for x in out_file:
            out_files_paths.extend(list(Path(out_dir).rglob(x)))
    elif isinstance(out_file, str):
        out_files_paths = list(Path(out_dir).rglob(out_file))

    print("the file searched:", out_file)
    print("the folder searched:", out_dir)

    print("out_files_path:", out_files_paths)
    print("in log prov, prov_record", prov_log_record)
    for path_file in out_files_paths:
        write_prov_file(prov_log_record, path_file)
    print("Provenance registered succesfully")
    return True
