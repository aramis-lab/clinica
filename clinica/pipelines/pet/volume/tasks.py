def apply_binary_mask_task(image: str, binary_mask: str) -> str:
    from pathlib import Path

    from clinica.pipelines.pet.volume.utils import apply_binary_mask

    return str(apply_binary_mask(Path(image), Path(binary_mask)))


def compute_atlas_statistics_task(image: str, atlas_names: list) -> list:
    from pathlib import Path

    from clinica.pipelines.pet.volume.utils import compute_atlas_statistics

    return [str(p) for p in compute_atlas_statistics(Path(image), atlas_names)]


def create_binary_mask_task(
    tissues: list,
    threshold: float = 0.3,
) -> str:
    from pathlib import Path

    from clinica.pipelines.pet.volume.utils import create_binary_mask

    return str(create_binary_mask([Path(tissue) for tissue in tissues], threshold))


def create_pvc_mask_task(tissues: list) -> str:
    from pathlib import Path

    from clinica.pipelines.pet.volume.utils import create_pvc_mask

    return str(create_pvc_mask([Path(tissue) for tissue in tissues]))


def normalize_to_reference_task(pet_image: str, region_mask: str) -> str:
    from pathlib import Path

    from clinica.pipelines.pet.volume.utils import normalize_to_reference

    return str(normalize_to_reference(Path(pet_image), Path(region_mask)))
