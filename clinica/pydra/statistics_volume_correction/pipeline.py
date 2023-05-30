import pydra
from pydra import Workflow

from clinica.pydra.engine import clinica_io


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
    from os.path import abspath, dirname, exists, join, pardir
    from typing import Any

    import numpy as np

    import clinica.pydra.statistics_volume_correction.task as utils
    from clinica.pydra.tasks import download_mni_template_2009a
    from clinica.utils.spm import spm_standalone_is_available, use_spm_standalone

    if spm_standalone_is_available():
        use_spm_standalone()

    query = {"pattern": parameters["t_map"] + "*", "description": "statistics t map"}

    input_spec = pydra.specs.SpecInfo(
        name="Input",
        fields=[
            ("_graph_checksums", Any),
            ("t_map", dict, query, {"mandatory": True}),
        ],
        bases=(pydra.specs.BaseSpec,),
    )
    wf = Workflow(name, input_spec=input_spec)

    for threshold in ("FWE", "FDR"):
        wf.add(
            utils.peak_correction_task(
                name=f"{threshold}_peak_correction_task",
                t_map=wf.lzin.t_map,
                t_threshold=parameters[f"{threshold}p"],
            )
        )
    for threshold in ("FWE", "FDR"):
        wf.add(
            utils.cluster_correction_task(
                name=f"{threshold}_cluster_correction_task",
                t_map=wf.lzin.t_map,
                t_thresh=parameters["height_threshold"],
                c_thresh=parameters[f"{threshold}c"],
            )
        )

    wf.add(download_mni_template_2009a(name="download_mni_template"))

    for threshold in ("FWE", "FDR"):
        for kind in ("peak", "cluster"):
            t_thresh_key = f"{threshold}p" if kind == "peak" else "height_threshold"
            c_thresh = parameters[f"{threshold}c"] if kind == "cluster" else np.nan
            wf.add(
                utils.produce_figures_task(
                    name=f"produce_figure_{threshold}_{kind}_correction",
                    nii_file=getattr(
                        wf, f"{threshold}_{kind}_correction_task"
                    ).lzout.nii_file,
                    template=wf.download_mni_template.lzout.mni_template_file,
                    type_of_correction=threshold,
                    t_thresh=parameters[t_thresh_key],
                    c_thresh=c_thresh,
                    n_cuts=parameters["n_cuts"],
                )
            )
            wf.add(
                utils.generate_output_task(
                    name=f"save_figure_{kind}_correction_{threshold}",
                    t_map=wf.lzin.t_map,
                    figs=getattr(
                        wf, f"produce_figure_{threshold}_{kind}_correction"
                    ).lzout.figs,
                    correction_name=f"{threshold}{kind[0]}",
                )
            )
    wf.set_output([("figs", wf.produce_figure_FDR_peak_correction.lzout.figs)])
    return wf
