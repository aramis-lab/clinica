from clinica.utils.testing_utils import build_caps_directory


def test_dwi_dti_info_loading(tmp_path):
    from clinica.pipelines.dwi.dti.pipeline import DwiDti

    caps = build_caps_directory(
        tmp_path / "caps",
        {"pipelines": {"dwi-dti": {}}, "subjects": {"sub-01": ["ses-M000"]}},
    )
    pipeline = DwiDti(caps_directory=str(caps))

    assert pipeline.info == {
        "id": "aramislab/dwi-dti",
        "author": ["Alexandre Routier", "Thomas Jacquemont"],
        "version": "0.1.1",
        "space_caps": "400M",
        "space_wd": "700M",
        "dependencies": [
            {"type": "software", "name": "fsl", "version": ">=5.0.9"},
            {"type": "software", "name": "ants", "version": ""},
            {"type": "software", "name": "mrtrix", "version": ""},
        ],
    }


def test_dwi_dti_dependencies(tmp_path, mocker):
    from packaging.specifiers import SpecifierSet
    from packaging.version import Version

    from clinica.pipelines.dwi.dti.pipeline import DwiDti
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
        {"pipelines": {"dwi-dti": {}}, "subjects": {"sub-01": ["ses-M000"]}},
    )
    pipeline = DwiDti(caps_directory=str(caps))

    assert pipeline.dependencies == [
        SoftwareDependency(
            ThirdPartySoftware.FSL, SpecifierSet(">=5.0.9"), Version("6.0.5.2")
        ),
        SoftwareDependency(
            ThirdPartySoftware.ANTS, SpecifierSet(">=0.0.0"), Version("2.3.2")
        ),
        SoftwareDependency(
            ThirdPartySoftware.MRTRIX, SpecifierSet(">=0.0.0"), Version("3.0.3")
        ),
    ]
