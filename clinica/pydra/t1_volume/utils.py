from os import PathLike
from typing import Union


def initialize_tissues_spm_segment(
    parameters: dict, tissue_probability_map: Union[PathLike, None] = None
) -> tuple:
    """Prepare data structure for SPM segment interface

    Parameters
    ----------
    tissue_probability_map: PathLike, optional
        Path to the nifti tissue probability map. Default=None.

    Returns
    -------
    tuple
        "tissues" data structure for SPMSegment
    """

    from clinica.pipelines.t1_volume_tissue_segmentation.t1_volume_tissue_segmentation_utils import (
        get_tissue_tuples,
    )
    from clinica.utils.spm import get_tpm

    parameters["tissue_probability_maps"] = tissue_probability_map or get_tpm()

    tissue_tuples = get_tissue_tuples(
        parameters["tissue_probability_maps"],
        parameters["tissue_classes"],
        parameters["dartel_tissues"],
        parameters["save_warped_unmodulated"],
        parameters["save_warped_modulated"],
    )
    return tissue_tuples
