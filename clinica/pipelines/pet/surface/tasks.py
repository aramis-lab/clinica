def remove_nan_from_image_task(image_path: str) -> str:
    from pathlib import Path

    from clinica.pipelines.pet.surface.utils import remove_nan_from_image

    return str(remove_nan_from_image(Path(image_path)))


def perform_gtmseg_task(
    caps_dir: str, subject_id: str, session_id: str, is_longitudinal: bool
) -> str:
    from pathlib import Path

    from clinica.pipelines.pet.surface.utils import perform_gtmseg

    return str(perform_gtmseg(Path(caps_dir), subject_id, session_id, is_longitudinal))


def make_label_conversion_task(gtmseg_file: str, csv_file: str) -> list:
    from pathlib import Path

    from clinica.pipelines.pet.surface.utils import make_label_conversion

    return [str(p) for p in make_label_conversion(Path(gtmseg_file), Path(csv_file))]


def run_apply_inverse_deformation_field_task(
    target_image: str,
    deformation_field: str,
    image: str,
    matscript_folder: str,
) -> str:
    from pathlib import Path

    from clinica.pipelines.pet.surface.utils import run_apply_inverse_deformation_field

    return str(
        run_apply_inverse_deformation_field(
            Path(target_image),
            Path(deformation_field),
            Path(image),
            Path(matscript_folder),
        )
    )


def run_apply_inverse_deformation_field_spm_standalone_task(
    target_image: str,
    deformation_field: str,
    image: str,
    matscript_folder: str,
) -> str:
    from pathlib import Path

    from clinica.pipelines.pet.surface.utils import (
        run_apply_inverse_deformation_field_spm_standalone,
    )

    return str(
        run_apply_inverse_deformation_field_spm_standalone(
            Path(target_image),
            Path(deformation_field),
            Path(image),
            Path(matscript_folder),
        )
    )


def normalize_suvr_task(pet_image: str, mask: str, output_dir=None) -> str:
    from pathlib import Path

    from clinica.pipelines.pet.surface.utils import normalize_suvr

    if output_dir:
        output_dir = Path(output_dir)
    return str(normalize_suvr(Path(pet_image), Path(mask), output_dir))


def reformat_surfname_task(
    hemisphere: str, left_surface: str, right_surface: str
) -> str:
    from pathlib import Path

    from clinica.pipelines.pet.surface.utils import reformat_surfname
    from clinica.utils.image import HemiSphere

    return str(
        reformat_surfname(
            HemiSphere(hemisphere),
            Path(left_surface),
            Path(right_surface),
        )
    )


def run_mris_expand_task(surface: str, output_dir=None):
    from pathlib import Path

    from clinica.pipelines.pet.surface.utils import run_mris_expand

    if output_dir:
        output_dir = Path(output_dir)
    return [str(p) for p in run_mris_expand(Path(surface), output_dir)]


def run_mri_surf2surf_task(
    surface: str,
    registration: str,
    gtmsegfile: str,
    subject_id: str,
    session_id: str,
    caps_dir: str,
    is_longitudinal: bool,
    output_dir=None,
) -> str:
    from pathlib import Path

    from clinica.pipelines.pet.surface.utils import run_mri_surf2surf

    if output_dir:
        output_dir = Path(output_dir)

    return str(
        run_mri_surf2surf(
            Path(surface),
            Path(registration),
            Path(gtmsegfile),
            subject_id,
            session_id,
            Path(caps_dir),
            is_longitudinal,
            output_dir,
        )
    )


def run_mri_vol2surf_task(
    pet_volume: str,
    surface: str,
    subject_id: str,
    session_id: str,
    caps_dir: str,
    gtmsegfile: str,
    is_longitudinal: bool,
    output_dir=None,
) -> str:
    from pathlib import Path

    from clinica.pipelines.pet.surface.utils import run_mri_vol2surf

    if output_dir:
        output_dir = Path(output_dir)

    return str(
        run_mri_vol2surf(
            Path(pet_volume),
            Path(surface),
            subject_id,
            session_id,
            Path(caps_dir),
            Path(gtmsegfile),
            is_longitudinal,
            output_dir,
        )
    )


def compute_weighted_mean_surface_task(surfaces, output_dir=None) -> str:
    from pathlib import Path

    from clinica.pipelines.pet.surface.utils import compute_weighted_mean_surface

    if output_dir:
        output_dir = Path(output_dir)

    return str(
        compute_weighted_mean_surface(
            [Path(surface) for surface in surfaces],
            output_dir,
        )
    )


def project_onto_fsaverage_task(
    projection: str,
    subject_id: str,
    session_id: str,
    caps_dir: str,
    fwhm: int,
    is_longitudinal: bool,
    output_dir=None,
) -> str:
    from pathlib import Path

    from clinica.pipelines.pet.surface.utils import project_onto_fsaverage

    if output_dir:
        output_dir = Path(output_dir)

    return str(
        project_onto_fsaverage(
            Path(projection),
            subject_id,
            session_id,
            Path(caps_dir),
            fwhm,
            is_longitudinal,
            output_dir,
        )
    )


def get_mid_surface_task(surfaces) -> str:
    from pathlib import Path

    from clinica.pipelines.pet.surface.utils import get_mid_surface

    return str(get_mid_surface([Path(surface) for surface in surfaces]))


def compute_average_pet_signal_based_on_annotations_task(
    pet_projections: tuple,
    atlas_files: dict,
    output_dir=None,
) -> tuple:
    from pathlib import Path

    from clinica.pipelines.pet.surface.utils import (
        compute_average_pet_signal_based_on_annotations,
    )
    from clinica.utils.stream import log_and_raise

    if output_dir:
        output_dir = Path(output_dir)
    if len(atlas_files) != 2:
        msg = (
            "The compute_average_pet_signal_based_on_annotations task requires two atlases "
            "for the argument 'atlas_files', one for desikan, one for destrieux. "
            f"The following {len(atlas_files)} were received:\n"
            + "\n".join([str(_) for _ in atlas_files])
        )
        log_and_raise(msg, ValueError)

    tsv_files = compute_average_pet_signal_based_on_annotations(
        (Path(pet_projections[0]), Path(pet_projections[1])),
        atlas_files,
        output_dir,
    )
    if len(tsv_files) != 2:
        msg = (
            "The compute_average_pet_signal_based_on_annotations task expects to return a tuple of 2 filenames. "
            f"Instead, {len(tsv_files)} files were computed:\n"
            + "\n".join([str(_) for _ in tsv_files])
        )
        log_and_raise(msg, ValueError)
    return str(tsv_files[0]), str(tsv_files[1])
