import typing as ty

from nipype.interfaces.spm import Smooth
from pydra.engine import Workflow
from pydra.tasks.nipype1.utils import Nipype1Task


def build_smoothing_workflow(
    name: str,
    fwhm: ty.Iterable[ty.Tuple[float, float, float]],
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
    Workflow
        The resulting smoothing workflow.
    """
    wf = Workflow(
        name,
        input_spec=["input_file"],
    )
    outputs = []
    for x, y, z in fwhm:
        fwhm_name = f"{'-'.join([str(int(f)) for f in (x, y, z)])}mm"
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
