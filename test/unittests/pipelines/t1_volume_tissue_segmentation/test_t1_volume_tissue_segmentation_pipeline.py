import shutil

from packaging.version import Version

from clinica.utils.bids import Visit
from clinica.utils.testing_utils import build_bids_directory, build_caps_directory


def test_t1_volume_tissue_segmentation_info_loading(tmp_path):
    from clinica.pipelines.t1_volume_tissue_segmentation.t1_volume_tissue_segmentation_pipeline import (
        T1VolumeTissueSegmentation,
    )

    caps = build_caps_directory(
        tmp_path / "caps",
        {
            "pipelines": {"t1-volume-tissue-segmentation": {}},
            "subjects": {"sub-01": ["ses-M000"]},
        },
    )
    pipeline = T1VolumeTissueSegmentation(caps_directory=str(caps))

    assert pipeline.info == {
        "id": "aramislab/t1-volume-tissue-segmentation",
        "author": "Jorge Samper-Gonzalez",
        "credits": ["Alexandre Routier"],
        "space_caps": "35M",
        "space_wd": "147M",
        "version": "0.1.0",
        "dependencies": [{"type": "software", "name": "spm", "version": ">=12"}],
    }


def test_t1_volume_tissue_segmentation_dependencies(tmp_path, mocker):
    from packaging.specifiers import SpecifierSet

    from clinica.pipelines.t1_volume_tissue_segmentation.t1_volume_tissue_segmentation_pipeline import (
        T1VolumeTissueSegmentation,
    )
    from clinica.utils.check_dependency import SoftwareDependency, ThirdPartySoftware

    mocker.patch(
        "clinica.utils.check_dependency._get_spm_version",
        return_value=Version("12.7219"),
    )
    caps = build_caps_directory(
        tmp_path / "caps",
        {
            "pipelines": {"t1-volume-tissue-segmentation": {}},
            "subjects": {"sub-01": ["ses-M000"]},
        },
    )
    pipeline = T1VolumeTissueSegmentation(caps_directory=str(caps))

    assert pipeline.dependencies == [
        SoftwareDependency(
            ThirdPartySoftware.SPM, SpecifierSet(">=12"), Version("12.7219")
        ),
    ]


def test_t1_volume_tissue_segmentation_get_processed_visits_empty(tmp_path, mocker):
    from clinica.pipelines.t1_volume_tissue_segmentation.t1_volume_tissue_segmentation_pipeline import (
        T1VolumeTissueSegmentation,
    )

    mocker.patch(
        "clinica.utils.check_dependency._get_spm_version",
        return_value=Version("12.7219"),
    )
    bids = build_bids_directory(
        tmp_path / "bids", {"sub-01": ["ses-M000", "ses-M006"], "sub-02": ["ses-M000"]}
    )
    caps = build_caps_directory(tmp_path / "caps", {})

    pipeline = T1VolumeTissueSegmentation(
        bids_directory=str(bids),
        caps_directory=str(caps),
        parameters={
            "tissue_classes": (1, 2, 3),
            "dartel_tissues": (1, 2, 3),
        },
    )
    assert pipeline.get_processed_visits() == []


def test_t1_volume_tissue_segmentation_get_processed_visits(tmp_path, mocker):
    """Test the get_processed_visits for T1VolumeTissueSegmentation.

    We build the following CAPS dataset:

    caps2
    ├── dataset_description.json
    └── subjects
        └── sub-01
            └── ses-M000
                └── t1
                    └── spm
                        └── segmentation
                            ├── dartel_input
                            │         ├── sub-01_ses-M000_T1w_segm-csf_dartelinput.nii.gz
                            │         ├── sub-01_ses-M000_T1w_segm-graymatter_dartelinput.nii.gz
                            │         └── sub-01_ses-M000_T1w_segm-whitematter_dartelinput.nii.gz
                            ├── native_space
                            │         ├── sub-01_ses-M000_T1w_segm-csf_probability.nii.gz
                            │         ├── sub-01_ses-M000_T1w_segm-graymatter_probability.nii.gz
                            │         └── sub-01_ses-M000_T1w_segm-whitematter_probability.nii.gz
                            └── normalized_space
                                ├── sub-01_ses-M000_T1w_segm-csf_space-Ixi549Space_modulated-off_probability.nii.gz
                                ├── sub-01_ses-M000_T1w_segm-graymatter_space-Ixi549Space_modulated-off_probability.nii.gz
                                └── sub-01_ses-M000_T1w_segm-whitematter_space-Ixi549Space_modulated-off_probability.nii.gz
    ... same for other subjects and sessions.
    Verify that removing files for a given visit will remove it from the "processed" visits.
    """
    from clinica.pipelines.t1_volume_tissue_segmentation.t1_volume_tissue_segmentation_pipeline import (
        T1VolumeTissueSegmentation,
    )

    mocker.patch(
        "clinica.utils.check_dependency._get_spm_version",
        return_value=Version("12.7219"),
    )
    bids = build_bids_directory(
        tmp_path / "bids",
        {"sub-01": ["ses-M000", "ses-M006"], "sub-02": ["ses-M000", "ses-M012"]},
    )
    caps = build_caps_directory(
        tmp_path / "caps",
        {
            "pipelines": {"t1_volume_tissue_segmentation": {}},
            "subjects": {
                "sub-01": ["ses-M006"],
                "sub-02": ["ses-M000", "ses-M012"],
            },
        },
    )
    pipeline = T1VolumeTissueSegmentation(
        bids_directory=str(bids),
        caps_directory=str(caps),
        parameters={
            "tissue_classes": (1, 2, 3),
            "dartel_tissues": (1, 2, 3),
        },
    )
    assert pipeline.get_processed_visits() == [
        Visit("sub-01", "ses-M006"),
        Visit("sub-02", "ses-M000"),
        Visit("sub-02", "ses-M012"),
    ]
    # Delete the folder "dartel_input" altogether for subject 02 session M000 (but keep the other folders)
    shutil.rmtree(
        tmp_path
        / "caps"
        / "subjects"
        / "sub-02"
        / "ses-M000"
        / "t1"
        / "spm"
        / "segmentation"
        / "dartel_input"
    )
    # Check that subject 02 session M000 is not considered a processed visit anymore
    assert pipeline.get_processed_visits() == [
        Visit("sub-01", "ses-M006"),
        Visit("sub-02", "ses-M012"),
    ]
    # Delete a single file in the "native_space" folder for subject 01 session M006 (keep other files and folders)
    (
        tmp_path
        / "caps"
        / "subjects"
        / "sub-01"
        / "ses-M006"
        / "t1"
        / "spm"
        / "segmentation"
        / "native_space"
        / "sub-01_ses-M006_T1w_segm-graymatter_probability.nii.gz"
    ).unlink()
    # Check that subject 01 session M006 is not considered a processed visit anymore
    assert pipeline.get_processed_visits() == [Visit("sub-02", "ses-M012")]
