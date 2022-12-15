import typing as ty

from nipype.interfaces.spm import Smooth
from pydra.engine import Workflow
from pydra.tasks.nipype1.utils import Nipype1Task


def build_smoothing_workflow(
    name: str,
    fwhm: ty.List[ty.List[float]],
) -> Workflow:
    """Build and parametrize a smoothing workflow.

    Parameters
    ----------
    name : str
        The name of the Workflow.

    fwhm : List[List[float]]
        The smoothing kernel(s) to use for smoothing.
        If multiple kernels are provided, the workflow
        will output files for each specified kernel.
        Each file name has the kernel size added.

    Returns
    -------
    wf : Workflow
        The resulting smoothing workflow.
    """
    wf = Workflow(
        name,
        input_spec=["input_file"],
    )
    outputs = []
    for fwhm_value in fwhm:
        if not isinstance(fwhm_value, list) or len(fwhm_value) != 3:
            raise ValueError(
                f"FWHM value is not valid. Expected a length 3 list and got: {fwhm_value}."
            )
        fwhm_name = f"{'-'.join([str(int(f)) for f in fwhm_value])}mm"
        task = Nipype1Task(
            name=f"smoothing_node_{fwhm_name}",
            interface=Smooth(),
        )
        task.inputs.fwhm = fwhm_value
        task.inputs.out_prefix = f"fwhm-{fwhm_name}_"
        task.inputs.in_files = wf.lzin.input_file
        wf.add(task)
        outputs.append(
            (
                "smoothed_files",
                getattr(wf, f"smoothing_node_{fwhm_name}").lzout.smoothed_files,
            )
        )
    wf.set_output(outputs)
    return wf
