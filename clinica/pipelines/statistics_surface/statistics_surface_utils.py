def get_t1_freesurfer_custom_file():
    import os

    custom_file = os.path.join(
        "@subject",
        "@session",
        "t1",
        "freesurfer_cross_sectional",
        "@subject_@session",
        "surf",
        "@hemi.thickness.fwhm@fwhm.fsaverage.mgh",
    )
    return custom_file


def get_pet_surface_custom_file(acq_label, suvr_reference_region):
    import os

    custom_file = os.path.join(
        "@subject",
        "@session",
        "pet",
        "surface",
        f"@subject_@session_trc-{acq_label}_pet"
        f"_space-fsaverage_suvr-{suvr_reference_region}_pvc-iy_hemi-@hemi_fwhm-@fwhm_projection.mgh",
    )
    return custom_file


def init_input_node(parameters, base_dir, subjects_visits_tsv):
    """Initialize the pipeline.

    This function will:
        - Create `surfstat_results_dir` in `base_dir`/<group_id> for SurfStat;
        - Save pipeline parameters in JSON file;
        - Copy TSV file with covariates;
        - Print begin execution message.
    """
    import json
    import os
    import shutil

    from clinica.utils.ux import print_begin_image

    group_id = "group-" + parameters["group_label"]

    # Create surfstat_results_dir for SurfStat
    surfstat_results_dir = os.path.join(base_dir, group_id)
    os.makedirs(surfstat_results_dir, exist_ok=True)

    # Save pipeline parameters in JSON file
    glm_dict = create_glm_info_dictionary(subjects_visits_tsv, parameters)
    json_filename = os.path.join(surfstat_results_dir, group_id + "_glm.json")
    with open(json_filename, "w") as json_file:
        json.dump(glm_dict, json_file, indent=4)

    # Copy TSV file with covariates
    tsv_filename = os.path.join(surfstat_results_dir, group_id + "_covariates.tsv")
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


def get_string_format_from_tsv(tsv_file):
    """
    Determine string format from TSV file.

    If the TSV file is like:

    participant_id  session_id  sex     group   age
    sub-CLNC0001    ses-M00     Female  CN      71.1
    sub-CLNC0002    ses-M00     Male    CN      81.3
    sub-CLNC0003    ses-M00     Male    CN      75.4

    The columns of the TSV file contains consecutively strings, strings,
    strings, strings and float. The string_format is therefore "%s %s %s %s %f".

    Args:
        tsv_file: TSV file.

    Returns:
        String formatting of the TSV file (e.g. "%s %s %s %s %f")
    """
    import pandas as pd

    demographics_df = pd.read_csv(tsv_file, sep="\t")

    def dtype_to_str_format(dtype):
        """Convert pandas dtypes (e.g. int64) to string format (e.g. %d)"""
        import numpy as np

        if dtype == np.int64:
            str_format = "%d"
        elif dtype == np.float64:
            str_format = "%f"
        elif dtype == np.object:
            str_format = "%s"
        else:
            raise ValueError("Unknown dtype (given: %s)" % dtype)
        return str_format

    list_str_format = []
    for column in demographics_df.columns:
        list_str_format.append(dtype_to_str_format(demographics_df[column].dtype))

    return " ".join(list_str_format)


def covariates_to_design_matrix(contrast, covariates=None):
    """
    Generate design matrix for SurfStat based on the contrast and the optional list of covariates.

    Design matrix "1 + <contrast> + <covariate_1> + ... + <covariate_n>"

    Example:
        >>> from clinica.pipelines.statistics_surface.statistics_surface_utils import covariates_to_design_matrix
        >>> covariates_to_design_matrix('group', 'age sex group')
        1 + group + age + sex
        >>> covariates_to_design_matrix('group', 'age')
        1 + group + age
        >>> covariates_to_design_matrix('group', None)
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
            design_matrix = "1 + " + " + ".join(
                covariate for covariate in list_covariates
            )
        else:
            design_matrix = (
                "1 + "
                + contrast
                + " + "
                + " + ".join(covariate for covariate in list_covariates)
            )
    else:
        design_matrix = "1 + " + contrast

    return design_matrix


def run_matlab(caps_dir, output_dir, subjects_visits_tsv, pipeline_parameters):
    """
    Wrap the call of SurfStat using clinicasurfstat.m Matlab script.

    Args:
        caps_dir (str): CAPS directory containing surface-based features
        output_dir (str): Output directory that will contain outputs of clinicasurfstat.m
        subjects_visits_tsv (str): TSV file containing the GLM information
        pipeline_parameters (dict): parameters of StatisticsSurface pipeline
    """
    import os

    from clinica.pipelines.statistics_surface.clinica_surfstat import clinica_surfstat
    from clinica.utils.check_dependency import check_environment_variable

    freesurfer_home = check_environment_variable("FREESURFER_HOME", "FreeSurfer")

    clinica_surfstat(
        os.path.join(caps_dir, "subjects"),
        output_dir,
        subjects_visits_tsv,
        covariates_to_design_matrix(
            pipeline_parameters["contrast"], pipeline_parameters["covariates"]
        ),
        pipeline_parameters["contrast"],
        pipeline_parameters["glm_type"],
        pipeline_parameters["group_label"],
        freesurfer_home,
        pipeline_parameters["custom_file"],
        pipeline_parameters["measure_label"],
        {
            "sizeoffwhm": pipeline_parameters["full_width_at_half_maximum"],
            "thresholduncorrectedpvalue": 0.001,
            "thresholdcorrectedpvalue": 0.05,
            "clusterthreshold": pipeline_parameters["cluster_threshold"],
        },
    )
    return output_dir


def create_glm_info_dictionary(tsv_file, pipeline_parameters):
    """Create dictionary containing the GLM information that will be stored in a JSON file."""
    out_dict = {
        # Clinica compulsory arguments
        "AnalysisType": pipeline_parameters["glm_type"],
        "DesignMatrix": covariates_to_design_matrix(
            pipeline_parameters["contrast"], pipeline_parameters["covariates"]
        ),
        "StringFormatTSV": get_string_format_from_tsv(tsv_file),
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
        out_dict["acq_label"] = pipeline_parameters["acq_label"]
        out_dict["suvr_reference_region"] = pipeline_parameters["suvr_reference_region"]

    return out_dict


def save_to_caps(source_dir, caps_dir, overwrite_caps, pipeline_parameters):
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

    Raise:
        NotImplementedError: If overwrite_caps=True.
    """
    import os
    import shutil

    from clinica.utils.ux import print_end_image

    group_id = "group-" + pipeline_parameters["group_label"]

    if pipeline_parameters["glm_type"] == "group_comparison":
        surfstat_folder = "surfstat_" + pipeline_parameters["glm_type"]
    elif pipeline_parameters["glm_type"] == "correlation":
        surfstat_folder = "surfstat_" + pipeline_parameters["glm_type"] + "_analysis"
    else:
        raise NotImplementedError(
            "The other GLM situations have not been implemented in this pipeline."
        )

    destination_dir = os.path.join(
        os.path.expanduser(caps_dir), "groups", group_id, "statistics", surfstat_folder
    )

    if overwrite_caps:
        raise NotImplementedError("save_to_caps(overwrite_caps=True) not implemented")
    shutil.copytree(source_dir, destination_dir, symlinks=True)
    print_end_image(group_id)
