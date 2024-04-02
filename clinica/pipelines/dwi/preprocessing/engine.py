"""This module contains the definition of the class DWIPreprocessingPipeline.

All DWI preprocessing pipelines should inherit from this class.
"""

from clinica.pipelines.engine import Pipeline


class DWIPreprocessingPipeline(Pipeline):
    """Class defining common properties of DWI preprocessing pipelines."""

    def _check_pipeline_parameters(self) -> None:
        """Check pipeline parameters."""
        from clinica.utils.stream import cprint

        self.parameters.setdefault("low_bval", 5)
        low_bval = self.parameters["low_bval"]
        if low_bval < 0:
            raise ValueError(
                f"The low_bval is negative ({low_bval}): it should be zero or close to zero."
            )
        if self.parameters["low_bval"] > 100:
            cprint(
                f"The low_bval parameter is {low_bval}: it should be close to zero.",
                lvl="warning",
            )
        self.parameters.setdefault("use_cuda", False)
        self.parameters.setdefault("initrand", False)

    def _check_custom_dependencies(self) -> None:
        """Check dependencies that can not be listed in the `info.json` file."""
        from clinica.utils.check_dependency import is_binary_present
        from clinica.utils.exceptions import ClinicaMissingDependencyError

        if self.parameters["use_cuda"]:
            if not is_binary_present("eddy_cuda"):
                raise ClinicaMissingDependencyError(
                    "[Error] FSL eddy with CUDA was set but Clinica could "
                    "not find eddy_cuda in your PATH environment. Check that "
                    "https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/eddy/UsersGuide#The_eddy_executables "
                    "is correctly set."
                )
