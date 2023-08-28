from os import PathLike
from pathlib import Path
from typing import Optional

from pydra.mark import annotate, task


@task
@annotate({"return": {"out_file": Path}})
def generate_acq_file_task(
    dwi_filename: PathLike,
    fsl_phase_encoding_direction: str,
    total_readout_time: str,
    image_id: Optional[str] = None,
) -> Path:
    from clinica.utils.dwi import generate_acq_file

    return Path(
        generate_acq_file(
            dwi_filename=str(dwi_filename),
            fsl_phase_encoding_direction=fsl_phase_encoding_direction,
            total_readout_time=total_readout_time,
            image_id=image_id,
        )
    )


@task
@annotate({"return": {"out_file": Path}})
def generate_index_file_task(
    b_values_filename: PathLike, image_id: Optional[str] = None
) -> Path:
    from clinica.utils.dwi import generate_index_file

    return Path(
        generate_index_file(
            b_values_filename=str(b_values_filename),
            image_id=image_id,
        )
    )
