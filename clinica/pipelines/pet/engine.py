from clinica.pipelines.engine import Pipeline
from clinica.utils.input_files import QueryPattern
from clinica.utils.pet import ReconstructionMethod, Tracer
from clinica.utils.stream import log_and_raise


class PETPipeline(Pipeline):
    def _check_pipeline_parameters(self) -> None:
        """Check pipeline parameters."""
        from clinica.utils.exceptions import ClinicaPipelineConfigurationError

        if "acq_label" not in self.parameters:
            log_and_raise(
                "Missing compulsory 'acq_label' key in pipeline parameter.",
                ClinicaPipelineConfigurationError,
            )
        self.parameters["acq_label"] = Tracer(self.parameters["acq_label"])
        if "reconstruction_method" in self.parameters:
            if self.parameters["reconstruction_method"]:
                self.parameters["reconstruction_method"] = ReconstructionMethod(
                    self.parameters["reconstruction_method"]
                )
        else:
            self.parameters["reconstruction_method"] = None

    def _get_pet_scans_query(self) -> QueryPattern:
        """Return the query to retrieve PET scans."""
        from clinica.utils.input_files import get_pet_nifti

        return get_pet_nifti(
            self.parameters["acq_label"], self.parameters["reconstruction_method"]
        )
