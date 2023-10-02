from os import PathLike
from pathlib import Path
from typing import Optional

from pydra.mark import annotate, task

from clinica.utils.dwi import generate_acq_file, generate_index_file

generate_acq_file_task = task(
    annotate({"return": {"out_file": Path}})(generate_acq_file)
)
generate_index_file_task = task(
    annotate({"return": {"out_file": Path}})(generate_index_file)
)
