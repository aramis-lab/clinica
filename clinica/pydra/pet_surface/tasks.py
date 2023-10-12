from pydra import Workflow


def gtm_segmentation() -> Workflow:
    """Parcellate the brain using the anatomical T1w image."""


def pet_mri_coregistration() -> Workflow:
    """Register the PET image to the anatomical T1w image."""


def intensity_normalization() -> Workflow:
    """Normalize PET image intensities using the SUVR. """
    ...


def partial_volume_correction() -> Workflow:
    """Perform partial volume correction of the normalized PET image."""
    ...


def pet_signal_projection() -> Workflow:
    """Project the PET signal onto the subject's cortical surface."""
    ...


def template_registration() -> Workflow:
    """Project the PET signal onto fsaverage with multiple smoothing levels."""
    ...


def atlas_statistics() -> Workflow:
    """Compute statistical measures for the Desikan and Destrieux atlases."""
    ...
