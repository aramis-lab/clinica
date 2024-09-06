from clinica.utils.testing_utils import build_caps_directory


def test_t1_volume_tissue_segmentation_info_loading(tmp_path):
    from clinica.pipelines.t1_volume_tissue_segmentation.t1_volume_tissue_segmentation_pipeline import (
        T1VolumeTissueSegmentation,
    )

    caps = build_caps_directory(
        tmp_path / "caps",
        {
            "pipelines": ["t1-volume-tissue-segmentation"],
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
    from packaging.version import Version

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
            "pipelines": ["t1-volume-tissue-segmentation"],
            "subjects": {"sub-01": ["ses-M000"]},
        },
    )
    pipeline = T1VolumeTissueSegmentation(caps_directory=str(caps))

    assert pipeline.dependencies == [
        SoftwareDependency(
            ThirdPartySoftware.SPM, SpecifierSet(">=12"), Version("12.7219")
        ),
    ]
