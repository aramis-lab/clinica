from clinica.utils.testing_utils import build_caps_directory


def test_dwi_preprocessing_using_fmap_info_loading(tmp_path):
    from clinica.pipelines.dwi.preprocessing.fmap.pipeline import (
        DwiPreprocessingUsingPhaseDiffFMap,
    )

    caps = build_caps_directory(
        tmp_path / "caps",
        {
            "pipelines": ["dwi-preprocessing-using-phasediff-fmap"],
            "subjects": {"sub-01": ["ses-M000"]},
        },
    )
    pipeline = DwiPreprocessingUsingPhaseDiffFMap(caps_directory=str(caps))

    assert pipeline.info == {
        "id": "aramislab/dwi-preprocessing-using-phasediff-fmap",
        "author": ["Alexandre Routier"],
        "credits": ["Nipype"],
        "version": "0.2.0",
        "space_caps": "250M",
        "space_wd": "6G",
        "dependencies": [
            {"type": "software", "name": "fsl", "version": ">=6.0.1"},
            {"type": "software", "name": "ants", "version": ""},
            {"type": "software", "name": "mrtrix", "version": ""},
        ],
    }


def test_dwi_preprocessing_using_fmap_dependencies(tmp_path, mocker):
    from packaging.specifiers import SpecifierSet
    from packaging.version import Version

    from clinica.pipelines.dwi.preprocessing.fmap.pipeline import (
        DwiPreprocessingUsingPhaseDiffFMap,
    )
    from clinica.utils.check_dependency import SoftwareDependency, ThirdPartySoftware

    mocker.patch(
        "clinica.utils.check_dependency._get_ants_version",
        return_value=Version("2.3.2"),
    )
    mocker.patch(
        "clinica.utils.check_dependency._get_fsl_version",
        return_value=Version("6.0.5.2"),
    )
    mocker.patch(
        "clinica.utils.check_dependency._get_mrtrix_version",
        return_value=Version("3.0.3"),
    )
    caps = build_caps_directory(
        tmp_path / "caps",
        {
            "pipelines": ["dwi-preprocessing-using-phasediff-fmap"],
            "subjects": {"sub-01": ["ses-M000"]},
        },
    )
    pipeline = DwiPreprocessingUsingPhaseDiffFMap(caps_directory=str(caps))

    assert pipeline.dependencies == [
        SoftwareDependency(
            ThirdPartySoftware.FSL, SpecifierSet(">=6.0.1"), Version("6.0.5.2")
        ),
        SoftwareDependency(
            ThirdPartySoftware.ANTS, SpecifierSet(">=0.0.0"), Version("2.3.2")
        ),
        SoftwareDependency(
            ThirdPartySoftware.MRTRIX, SpecifierSet(">=0.0.0"), Version("3.0.3")
        ),
    ]
