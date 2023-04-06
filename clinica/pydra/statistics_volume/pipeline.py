import pydra
from nipype.algorithms.misc import Gunzip
from pydra import Workflow
from pydra.tasks.nipype1.utils import Nipype1Task

from clinica.pydra.engine import clinica_io


def _check_pipeline_parameters(parameters: dict) -> dict:
    """Check the parameters passed to the pipeline.

    Parameters
    ----------
    parameters : dict
        Dictionary of parameters to analyze.

    Returns
    -------
    dict :
        Cleaned dictionary of parameters.
    """
    from clinica.utils.exceptions import ClinicaException

    # PET-Volume pipeline
    if parameters["orig_input_data_volume"] == "pet-volume":
        for param in ("acq_label", "suvr_reference_region"):
            if not parameters[param]:
                raise ClinicaException(
                    f"You selected pet-volume pipeline without setting --{param} flag. "
                    "Clinica will now exit."
                )

    # Custom pipeline
    if parameters["orig_input_data_volume"] == "custom-pipeline":
        if not all([parameters["custom_file"], parameters["measure_label"]]):
            raise ClinicaException(
                "You must set --measure_label and --custom_file flags."
            )
        if not parameters["custom_file"]:
            raise ClinicaException(
                "Custom pipeline was selected but no 'custom_file' was specified."
            )
    return parameters


def _build_query(parameters: dict) -> dict:
    input_name = parameters["orig_input_data_volume"].replace("-", "_")
    query = {
        "group_label": parameters["group_label_dartel"],
        "fwhm": parameters["full_width_at_half_maximum"],
    }
    if input_name == "pet_volume":
        query.update(
            {
                "acq_label": parameters["acq_label"],
                "suvr_reference_region": parameters["suvr_reference_region"],
                "use_brainmasked_image": True,
                "use_pvc_data": parameters["use_pvc_data"],
            }
        )
    elif input_name == "t1_volume":
        query.update({"tissue_number": 1, "modulation": True})

    elif input_name == "custom_pipeline":
        query = {
            "pattern": parameters["custom_file"],
            "description": "custom file provided by user",
        }
    return query


@clinica_io
def build_core_workflow(name: str = "core", parameters={}) -> Workflow:
    """Build the core workflow for the Statistics Volume pipeline.

    Parameters
    ----------
    name : str, optional
        The name of the workflow. Default="core".

    parameters : dict, optional
        Dictionary of parameters to be used
        within the workflow.
        Default={}.

    Returns
    -------
    wf : Workflow
        The core workflow.
    """
    from os.path import dirname, join
    from typing import Any

    import clinica.pydra.statistics_volume.task as utils
    from clinica.utils.exceptions import ClinicaException
    from clinica.utils.filemanip import get_parent
    from clinica.utils.spm import spm_standalone_is_available, use_spm_standalone

    if spm_standalone_is_available():
        use_spm_standalone()

    parameters = _check_pipeline_parameters(parameters)

    input_name = parameters["orig_input_data_volume"].replace("-", "_")
    query = _build_query(parameters)

    input_spec = pydra.specs.SpecInfo(
        name="Input",
        fields=[
            ("_graph_checksums", Any),
            (input_name, dict, query, {"mandatory": True}),
        ],
        bases=(pydra.specs.BaseSpec,),
    )
    wf = Workflow(name, input_spec=input_spec)
    wf.add(
        Nipype1Task(
            name=f"unzip_nii",
            interface=Gunzip(),
            in_file=getattr(wf.lzin, input_name),
        )
        .split("in_file")
        .combine("in_file")
    )
    wf.add(
        utils.get_group_1_and_2_task(
            name="get_groups",
            tsv=parameters["tsv_file"],
            contrast=parameters["contrast"],
        )
    )
    # 1. Model creation
    # We use overwrite option to be sure this node is always run so that it can delete the output dir if it
    # already exists (this may cause error in output files otherwise)
    wf.add(
        utils.write_matlab_model_task(
            name="model_creation",
            tsv=parameters["tsv_file"],
            contrast=parameters["contrast"],
            template_file=join(
                get_parent(__file__, n=3),
                "pipelines",
                "statistics_volume",
                "template_model_creation.m",
            ),
            file_list=wf.unzip_nii.lzout.out_file,
            idx_group1=wf.get_groups.lzout.idx_group1,
            idx_group2=wf.get_groups.lzout.idx_group2,
        )
    )
    wf.add(
        utils.run_m_script_task(
            name="run_spm_model_creation",
            m_file=wf.model_creation.lzout.current_model,
        )
    )
    # 2. Model estimation
    wf.add(
        utils.clean_template_file_task(
            name="model_estimation",
            mat_file=wf.run_spm_model_creation.lzout.output_mat_file,
            template_file=join(
                get_parent(__file__, n=3),
                "pipelines",
                "statistics_volume",
                "template_model_estimation.m",
            ),
        )
    )
    wf.add(
        utils.run_m_script_task(
            name="run_spm_model_estimation",
            m_file=wf.model_estimation.lzout.current_model_estimation,
        )
    )
    # 3. Contrast
    wf.add(
        utils.clean_spm_contrast_file_task(
            name="model_contrast",
            mat_file=wf.run_spm_model_estimation.lzout.output_mat_file,
            template_file=join(
                get_parent(__file__, n=3),
                "pipelines",
                "statistics_volume",
                "template_model_contrast.m",
            ),
            covariates=wf.model_creation.lzout.covariates,
            class_names=wf.get_groups.lzout.class_names,
        )
    )
    wf.add(
        utils.run_m_script_task(
            name="run_spm_model_contrast",
            m_file=wf.model_contrast.lzout.current_model_estimation,
        )
    )
    # 4. Results
    wf.add(
        utils.clean_spm_result_file_task(
            name="model_result_no_correction",
            mat_file=wf.run_spm_model_contrast.lzout.output_mat_file,
            template_file=join(
                get_parent(__file__, n=3),
                "pipelines",
                "statistics_volume",
                "template_model_results.m",
            ),
            method="none",
            threshold=parameters["cluster_threshold"],
        )
    )
    wf.add(
        utils.run_m_script_task(
            name="run_spm_model_result_no_correction",
            m_file=wf.model_result_no_correction.lzout.current_model_result,
        )
    )
    wf.add(
        utils.copy_and_rename_spm_output_files_task(
            name="read_output_node",
            spm_mat=wf.run_spm_model_result_no_correction.lzout.output_mat_file,
            class_names=wf.get_groups.lzout.class_names,
            covariates=wf.model_creation.lzout.covariates,
            group_label=parameters["group_label"],
            fwhm=parameters["full_width_at_half_maximum"],
            measure=parameters["measure_label"],
        )
    )
    wf.set_output(
        [
            ("spmT_0001", wf.read_output_node.lzout.spmT_0001),
            ("spmT_0002", wf.read_output_node.lzout.spmT_0002),
            ("spm_figures", wf.read_output_node.lzout.spm_figures),
            ("variance_of_error", wf.read_output_node.lzout.variance_of_error),
            ("resels_per_voxels", wf.read_output_node.lzout.resels_per_voxels),
            ("mask", wf.read_output_node.lzout.mask),
            ("regression_coeff", wf.read_output_node.lzout.regression_coeff),
            ("contrasts", wf.read_output_node.lzout.contrasts),
        ]
    )
    return wf
