import pydra
from pydra import Workflow

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
    if parameters["orig_input_data_ml"] == "pet-volume":
        for param in ("acq_label", "suvr_reference_region"):
            if not parameters[param]:
                raise ClinicaException(
                    f"You selected pet-volume pipeline without setting --{param} flag. "
                    "Clinica will now exit."
                )
    return parameters


def _build_query(parameters: dict) -> dict:
    if parameters["orig_input_data_ml"] == "t1-volume":
        return {
            "modulation": True,
            "tissue_number": 1,
            "group_label": parameters["group_label"],
        }
    if parameters["orig_input_data_ml"] == "pet-volume":
        return {
            "acq_label": parameters["acq_label"],
            "suvr_reference_region": parameters["suvr_reference_region"],
            "use_brainmasked_image": False,
            "use_pvc_data": parameters["use_pvc_data"],
            "fwhm": 0,
            "group_label": parameters["group_label"],
        }


@clinica_io
def build_core_workflow(name: str = "core", parameters={}) -> Workflow:
    """Build the core workflow for the machine learning spatial svm pipeline.

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
        The core
    """
    from typing import Any

    import clinica.pydra.machine_learning_spatial_svm.tasks as utils

    parameters = _check_pipeline_parameters(parameters)
    query = _build_query(parameters)
    input_name = parameters["orig_input_data_ml"].replace("-", "_")
    input_spec = pydra.specs.SpecInfo(
        name="Input",
        fields=[
            ("_graph_checksums", Any),
            (
                input_name,
                dict,
                query,
                {"mandatory": True},
            ),
            (
                "dartel_template",
                dict,
                {"group_label": parameters["group_label"]},
                {"mandatory": True},
            ),
        ],
        bases=(pydra.specs.BaseSpec,),
    )
    wf = Workflow(name, input_spec=input_spec)

    wf.split((input_name))

    wf.add(
        utils.obtain_g_fisher_tensor_task(
            name="obtain_g_fisher_tensor_task",
            dartel_input=wf.lzin.dartel_template,
            FWHM=parameters["fwhm"],
        )
    )
    wf.add(
        utils.obtain_time_step_estimation_task(
            name="obtain_time_step_estimation_task",
            dartel_input=wf.lzin.dartel_template,
            g=wf.obtain_g_fisher_tensor_task.lzout.g_fisher_tensor,
            FWHM=parameters["fwhm"],
        )
    )
    wf.add(
        utils.heat_solver_equation_task(
            name="heat_solver_equation_task",
            input_image=getattr(wf.lzin, input_name),
            g=wf.obtain_g_fisher_tensor_task.lzout.g_fisher_tensor,
            FWHM=parameters["fwhm"],
            t_step=wf.obtain_time_step_estimation_task.lzout.t_step,
            dartel_input=wf.lzin.dartel_template,
        )
    )
    wf.set_output(
        [
            (
                "fisher_tensor_path",
                wf.obtain_g_fisher_tensor_task.lzout.output_fisher_tensor,
            ),
            ("json_file", wf.obtain_time_step_estimation_task.lzout.output_data),
            ("regularized_image", wf.heat_solver_equation_task.lzout.regularized),
        ]
    )
    return wf
