from clinica.pipelines.engine import Pipeline


class PETPipeline(Pipeline):
    def _check_pipeline_parameters(self) -> None:
        """Check pipeline parameters."""
        if "acq_label" not in self.parameters.keys():
            raise KeyError("Missing compulsory acq_label key in pipeline parameter.")
        self.parameters.setdefault("reconstruction_method", None)

    def _get_pet_scans_query(self) -> dict:
        """Return the query to retrieve PET scans."""
        from clinica.utils.input_files import bids_pet_nii
        from clinica.utils.pet import ReconstructionMethod, Tracer

        pet_tracer = None
        if self.parameters["acq_label"] is not None:
            pet_tracer = Tracer(self.parameters["acq_label"])

        reconstruction_method = None
        if self.parameters["reconstruction_method"] is not None:
            reconstruction_method = ReconstructionMethod(
                self.parameters["reconstruction_method"]
            )

        return bids_pet_nii(pet_tracer, reconstruction_method)
