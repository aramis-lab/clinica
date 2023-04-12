from typing import List

from pydra.mark import annotate, task


@task
@annotate({"return": {"nii_file": str}})
def peak_correction_task(
    t_map: str, t_threshold: float, output_name: str = None
) -> str:
    """Pydra task for correcting the peak.

    ..note::
        Please refer to the documentation of function
        `clinica.pipelines.statistics_volume_correction.statistics_volume_correction_utils.peak_correction`.
    """
    from clinica.pipelines.statistics_volume_correction.statistics_volume_correction_utils import (
        peak_correction,
    )

    return peak_correction(t_map, t_threshold, output_name)


@task
@annotate({"return": {"nii_file": str}})
def cluster_correction_task(
    t_map: str, t_thresh: float, c_thresh: int, output_name: str = None
) -> str:
    """Pydra task for cluster correction

    ..note::
        Please refer to the documentation of function
        `clinica.pipelines.statistics_volume_correction.statistics_volume_correction_utils.cluster_correction`.
    """
    from clinica.pipelines.statistics_volume_correction.statistics_volume_correction_utils import (
        cluster_correction,
    )

    return cluster_correction(t_map, t_thresh, c_thresh, output_name)


@task
@annotate({"return": {"figs": list}})
def produce_figures_task(
    nii_file: str,
    template: str,
    type_of_correction: str,
    t_thresh: str,
    c_thresh: int,
    n_cuts: int,
) -> List[str]:
    """Pydra task for figure production

    ..note::
        Please refer to the documentation of function
        `clinica.pipelines.statistics_volume_correction.statistics_volume_correction_utils.produce_figures`
    """
    from clinica.pipelines.statistics_volume_correction.statistics_volume_correction_utils import (
        produce_figures,
    )

    return produce_figures(
        nii_file, template, type_of_correction, t_thresh, c_thresh, n_cuts
    )


@task
@annotate({"return": {"bool": bool}})
def generate_output_task(t_map: str, figs: List[str], correction_name: str) -> None:
    """Pydra task for output production

    ..note::
        Please refer to the documentation of function
        `clinica.pipelines.statistics_volume_correction.statistics_volume_correction_utils.generate_output`
    """
    from clinica.pipelines.statistics_volume_correction.statistics_volume_correction_utils import (
        generate_output,
    )

    return generate_output(t_map, figs, correction_name)
