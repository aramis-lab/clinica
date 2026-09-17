import subprocess
from pathlib import Path
from typing import Optional

from clinica.utils.exceptions import ClinicaSubprocessError
from clinica.utils.stream import cprint, log_and_raise

__all__ = ["run_command_as_subprocess"]


def run_command_as_subprocess(
    command_name: str,
    command: str,
    out_file: Optional[Path] = None,
    env: Optional[dict] = None,
):
    cprint(
        f"Running {command_name} with the following command:\n\n{command}", lvl="debug"
    )
    subprocess_ = subprocess.run(
        command,
        shell=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        env=env,
    )
    if (code := subprocess_.returncode) != 0:
        error_msg = (
            f"The subprocess '{command_name}' failed with a non-zero return code of {code}. "
            f"The following command was run:\n\n{command}"
        )
        log_and_raise(error_msg, ClinicaSubprocessError)
    if out_file and not out_file.exists():
        error_msg = (
            f"Something went wrong while trying to run '{command_name}' ."
            f"Expected output file was not generated. Command launched :\n\t " + command
        )
        log_and_raise(error_msg, ClinicaSubprocessError)
