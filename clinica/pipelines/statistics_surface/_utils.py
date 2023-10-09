from pathlib import Path
from typing import Optional

__all__ = [
    "init_input_node",
    "run_clinica_surfstat",
    "save_to_caps",
    "get_t1_freesurfer_custom_file",
    "get_pet_surface_custom_file",
    "create_glm_info_dictionary",
    "build_design_matrix",
]


def get_t1_freesurfer_custom_file() -> str:
    import os

    return os.path.join(
        "%(subject)s",
        "%(session)s",
        "t1",
        "freesurfer_cross_sectional",
        "%(subject)s_%(session)s",
        "surf",
        "%(hemi)s.thickness.fwhm%(fwhm)s.fsaverage.mgh",
    )


def get_pet_surface_custom_file(acq_label: str, suvr_reference_region: str) -> str:
    import os

    return os.path.join(
        "%(subject)s",
        "%(session)s",
        "pet",
        "surface",
        f"%(subject)s_%(session)s_trc-{acq_label}_pet"
        f"_space-fsaverage_suvr-{suvr_reference_region}_pvc-iy_hemi-%(hemi)s_fwhm-%(fwhm)s_projection.mgh",
    )


def init_input_node(parameters: dict, base_dir, subjects_visits_tsv):
    """Initialize the pipeline.

    This function will:
        - Create `surfstat_results_dir` in `base_dir`/<group_id> for SurfStat;
        - Save pipeline parameters in JSON file;
        - Copy TSV file with covariates;
        - Print begin execution message.

    Parameters
    ----------
    parameters : dict
        The pipeline's parameters.

    base_dir : Path
        The path to the pipeline's base directory.
        This is a pathlib Path. No type hints because of Nipype.

    subjects_visits_tsv : Path
        The path to the subjects TSV file.
        This is a pathlib Path. No type hints because of Nipype.

    Returns
    -------
    group_label : str
        The group label.

    surfstat_results_dir : Path
        The folder which will contain the results for SurfStat.
    """
    import json
    import shutil

    from clinica.pipelines.statistics_surface._utils import create_glm_info_dictionary
    from clinica.utils.ux import print_begin_image

    group_id = "group-" + parameters["group_label"]
    surfstat_results_dir = base_dir / group_id
    surfstat_results_dir.mkdir(parents=True, exist_ok=True)

    # Save pipeline parameters in JSON file
    glm_dict = create_glm_info_dictionary(subjects_visits_tsv, parameters)
    with open(surfstat_results_dir / f"{group_id}_glm.json", "w") as json_file:
        json.dump(glm_dict, json_file, indent=4)

    # Copy TSV file with covariates
    tsv_filename = surfstat_results_dir / f"{group_id}_covariates.tsv"
    shutil.copyfile(subjects_visits_tsv, tsv_filename)

    # Print begin message
    list_keys = ["AnalysisType", "Covariates", "Contrast", "FWHM", "ClusterThreshold"]
    list_values = [
        parameters["glm_type"],
        parameters["covariates"],
        parameters["contrast"],
        str(parameters["full_width_at_half_maximum"]),
        str(parameters["cluster_threshold"]),
    ]
    group_id = "group-" + parameters["group_label"]
    print_begin_image(group_id, list_keys, list_values)

    return parameters["group_label"], surfstat_results_dir


def _get_string_format_from_tsv(tsv_file: Path) -> str:
    """Determine string format from TSV file.

    If the TSV file is like:

    participant_id  session_id  sex     group   age
    sub-CLNC0001    ses-M000     Female  CN      71.1
    sub-CLNC0002    ses-M000     Male    CN      81.3
    sub-CLNC0003    ses-M000     Male    CN      75.4

    The columns of the TSV file contains consecutively strings, strings,
    strings, strings and float. The string_format is therefore "%s %s %s %s %f".

    Parameters
    ----------
    tsv_file : Path
        The path to the TSV file.

    Returns
    -------
    str :
        The string formatting of the TSV file (e.g. "%s %s %s %s %f")
    """
    import pandas as pd

    demographics_df = pd.read_csv(tsv_file, sep="\t")

    return " ".join(
        [
            _convert_dtype_to_str_format(demographics_df[column].dtype)
            for column in demographics_df.columns
        ]
    )


def _convert_dtype_to_str_format(dtype) -> str:
    """Convert pandas dtypes (e.g. int64) to string format (e.g. %d)"""
    import numpy as np

    if dtype == np.int64:
        return "%d"
    if dtype == np.float64:
        return "%f"
    if dtype == np.object:
        return "%s"
    raise ValueError(f"Unknown dtype (given: {dtype})")


def build_design_matrix(contrast: str, covariates: Optional[str] = None) -> str:
    """Generate the design matrix for SurfStat based on the contrast and the optional list of covariates.

    Design matrix "1 + <contrast> + <covariate_1> + ... + <covariate_n>"

    Example
    -------
    >>> from clinica.pipelines.statistics_surface.statistics_surface_utils import _build_design_matrix
    >>> _build_design_matrix('group', 'age sex group')
    1 + group + age + sex
    >>> _build_design_matrix('group', 'age')
    1 + group + age
    >>> _build_design_matrix('group', None)
    1 + group
    """
    if covariates:
        # Convert string to list while handling case where several spaces are present
        list_covariates = list(covariates)
        try:
            list_covariates.remove("")
        except ValueError:
            pass
        if contrast in list_covariates:
            return "1 + " + " + ".join(covariate for covariate in list_covariates)
        return (
            "1 + "
            + contrast
            + " + "
            + " + ".join(covariate for covariate in list_covariates)
        )
    return "1 + " + contrast


def run_clinica_surfstat(
    caps_dir,  # Path. no type hint because of Nipype self-contained requirement
    output_dir,  # Path. no type hint because of Nipype self-contained requirement
    subjects_visits_tsv,  # Path. no type hint because of Nipype self-contained requirement
    pipeline_parameters: dict,
):
    """Call clinica_surfstat function.

    Parameters
    ----------
    caps_dir : Path
        The path to CAPS directory containing surface-based features.

    output_dir : Path
        The path to output directory that will contain outputs.

    subjects_visits_tsv : Path
        The path to TSV file containing the GLM information.

    pipeline_parameters : dict
        Parameters of StatisticsSurface pipeline.

    Returns
    -------
    output_dir : Path
        The path to the output directory.
    """
    from pathlib import Path

    from clinica.pipelines.statistics_surface._utils import build_design_matrix
    from clinica.pipelines.statistics_surface.surfstat import clinica_surfstat
    from clinica.utils.check_dependency import check_environment_variable

    freesurfer_home = Path(check_environment_variable("FREESURFER_HOME", "FreeSurfer"))

    clinica_surfstat(
        caps_dir / "subjects",
        output_dir,
        subjects_visits_tsv,
        build_design_matrix(
            pipeline_parameters["contrast"], pipeline_parameters["covariates"]
        ),
        pipeline_parameters["contrast"],
        pipeline_parameters["glm_type"],
        pipeline_parameters["group_label"],
        freesurfer_home,
        pipeline_parameters["measure_label"],
        surface_file=pipeline_parameters["custom_file"],
        fwhm=pipeline_parameters["full_width_at_half_maximum"],
        cluster_threshold=pipeline_parameters["cluster_threshold"],
    )
    return output_dir


def create_glm_info_dictionary(tsv_file: Path, pipeline_parameters: dict) -> dict:
    """Create dictionary containing the GLM information that will be stored in a JSON file."""
    glm_info = {
        # Clinica compulsory arguments
        "AnalysisType": pipeline_parameters["glm_type"],
        "DesignMatrix": build_design_matrix(
            pipeline_parameters["contrast"],
            pipeline_parameters["covariates"],
        ),
        "StringFormatTSV": _get_string_format_from_tsv(tsv_file),
        "Contrast": pipeline_parameters["contrast"],
        "GroupLabel": pipeline_parameters["group_label"],
        # Optional arguments
        "Covariates": pipeline_parameters["covariates"],
        "FWHM": pipeline_parameters["full_width_at_half_maximum"],
        # Optional arguments for custom pipeline
        "custom_file": pipeline_parameters["custom_file"],
        "measure_label": pipeline_parameters["measure_label"],
        # Advanced arguments (i.e. tricky parameters)
        "ThresholdUncorrectedPvalue": 0.001,
        "ThresholdCorrectedPvalue": 0.05,
        "ClusterThreshold": pipeline_parameters["cluster_threshold"],
    }
    # Optional arguments for inputs from pet-surface pipeline
    if (
        pipeline_parameters["acq_label"]
        and pipeline_parameters["suvr_reference_region"]
    ):
        glm_info["acq_label"] = pipeline_parameters["acq_label"]
        glm_info["suvr_reference_region"] = pipeline_parameters["suvr_reference_region"]

    return glm_info


def save_to_caps(
    source_dir,  # Path. no type hint because of Nipype self-contained requirement
    caps_dir,  # Path. no type hint because of Nipype self-contained requirement
    overwrite_caps: bool,
    group_label: str,
    glm_type: str,
) -> None:
    """Save `source_dir`/ to CAPS folder.

    This function copies outputs of `source_dir`/ to
    `caps_dir`/groups/<group_id>/<statistics>/surfstat_<glm_type>/

    The `source_dir`/ folder should contain the following elements:
        - group-<group_label>_<group_1_or_2>-lt-<group_1_or_2>_measure-<measure>_fwhm-<label>_suffix.ext
    or
        - group-<group_label>_correlation-<label>_contrast-{-|+}_measure-<measure>_fwhm-<label>_suffix.ext
    and
        - group-<group_label>_covariates.tsv
        - group-<group_label>_glm.json

    Raises
    ------
    NotImplementedError
        If overwrite_caps=True.
    """
    import shutil

    from clinica.pipelines.statistics_surface.surfstat.models import GLMModelType
    from clinica.utils.ux import print_end_image

    glm_type = GLMModelType(glm_type)
    if glm_type == GLMModelType.GROUP_COMPARISON:
        surfstat_folder = f"surfstat_{glm_type}"
    if glm_type == GLMModelType.CORRELATION:
        surfstat_folder = f"surfstat_{glm_type}_analysis"

    destination_dir = (
        caps_dir.expanduser()
        / "groups"
        / f"group-{group_label}"
        / "statistics"
        / surfstat_folder
    )

    if overwrite_caps:
        raise NotImplementedError("save_to_caps(overwrite_caps=True) not implemented")
    shutil.copytree(source_dir, destination_dir, symlinks=True)
    print_end_image(f"group-{group_label}")
