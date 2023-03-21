import pydra
from nipype.algorithms.misc import Gunzip
from pydra import Workflow
from pydra.tasks.nipype1.utils import Nipype1Task

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
    from clinica.utils.exceptions import ClinicaException
    from clinica.utils.inputs import RemoteFileStructure, fetch_file
    from clinica.utils.spm import spm_standalone_is_available, use_spm_standalone
    from clinica.utils.stream import cprint

    if spm_standalone_is_available():
        use_spm_standalone()

    # parameters = _check_pipeline_parameters(parameters)

    input_name = "t_map"
    query = {"pattern": parameters["t_map"] + "*", "description": "statistics t map"}

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
        utils.peak_correction_task(
            name="FWE_peak_correction_task",
            t_map=wf.lzin.t_map,
            t_threshold=parameters["FWEp"],
        )
    )
    wf.add(
        utils.peak_correction_task(
            name="FDR_peak_correction_task",
            t_map=wf.lzin.t_map,
            t_threshold=parameters["FDRp"],
        )
    )
    wf.add(
        utils.cluster_correction_task(
            name="FWE_cluster_correction_task",
            t_map=wf.lzin.t_map,
            t_thresh=parameters["height_threshold"],
            c_thresh=parameters["FWEc"],
        )
    )
    wf.add(
        utils.cluster_correction_task(
            name="FDR_cluster_correction_task",
            t_map=wf.lzin.t_map,
            t_thresh=parameters["height_threshold"],
            c_thresh=parameters["FDRc"],
        )
    )
    root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
    path_to_mask = join(root, "resources", "masks")
    url_aramis = "https://aramislab.paris.inria.fr/files/data/img_t1_linear/"
    FILE1 = RemoteFileStructure(
        filename="mni_icbm152_t1_tal_nlin_sym_09a.nii.gz",
        url=url_aramis,
        checksum="3b244ee7e287319d36a25263744c468ef0ab2fe5a94b15a2138844db73b49adf",
    )
    if not (exists(join(path_to_mask, FILE1.filename))):
        try:
            fetch_file(FILE1, path_to_mask)
        except IOError as err:
            cprint(
                msg=f"Unable to download required template (mni_icbm152) for processing: {err}",
                lvl="error",
            )

    wf.add(
        utils.produce_figures_task(
            name="produce_figure_FWE_peak_correction",
            nii_file=wf.FWE_peak_correction_task.nii_file,
            template=join(path_to_mask, FILE1.filename),
            type_of_correction="FWE",
            t_thresh=parameters["FWEp"],
            c_thresh=np.nan,
            n_cuts=parameters["n_cuts"],
        )
    )
    wf.add(
        utils.produce_figures_task(
            name="produce_figure_FDR_peak_correction",
            nii_file=wf.FDR_peak_correction_task.nii_file,
            template=join(path_to_mask, FILE1.filename),
            type_of_correction="FDR",
            t_thresh=parameters["FDRp"],
            c_thresh=np.nan,
            n_cuts=parameters["n_cuts"],
        )
    )

    wf.add(
        utils.produce_figures_task(
            name="produce_figure_FWE_cluster_correction",
            nii_file=wf.FWE_cluster_correction_task.nii_file,
            template=join(path_to_mask, FILE1.filename),
            type_of_correction="FWE",
            t_thresh=parameters["height_threshold"],
            c_thresh=parameters["FWEc"],
            n_cuts=parameters["n_cuts"],
        )
    )
    wf.add(
        utils.produce_figures_task(
            name="produce_figure_FDR_cluster_correction",
            nii_file=wf.FDR_cluster_correction_task.nii_file,
            template=join(path_to_mask, FILE1.filename),
            type_of_correction="FDR",
            t_thresh=parameters["height_threshold"],
            c_thresh=parameters["FDRc"],
            n_cuts=parameters["n_cuts"],
        )
    )
    wf.add(
        utils.generate_output_task(
            name="save_figure_peak_correction_FWE",
            t_map=wf.lzin.t_map,
            figs=wf.produce_figure_FWE_peak_correction.lzout.figs,
            name="FWEp",
        )
    )
    wf.add(
        utils.generate_output_task(
            name="save_figure_peak_correction_FDR",
            t_map=wf.lzin.t_map,
            figs=wf.produce_figure_FDR_peak_correction.lzout.figs,
            name="FDRp",
        )
    )
    wf.add(
        utils.generate_output_task(
            name="save_figure_peak_correction_FWE",
            t_map=wf.lzin.t_map,
            figs=wf.produce_figure_FWE_cluster_correction.lzout.figs,
            name="FWEc",
        )
    )
    wf.add(
        utils.generate_output_task(
            name="save_figure_peak_correction_FWE",
            t_map=wf.lzin.t_map,
            figs=wf.produce_figure_FDR_cluster_correction.lzout.figs,
            name="FDRc",
        )
    )
    return wf
