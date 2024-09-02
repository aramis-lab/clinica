from clinica.utils.testing_utils import build_caps_directory


def test_dwi_connectome_info_loading(tmp_path):
    from clinica.pipelines.dwi.connectome.pipeline import DwiConnectome

    caps = build_caps_directory(
        tmp_path / "caps",
        {"pipelines": ["compute-atlas"], "subjects": {"sub-01": ["ses-M000"]}},
    )
    pipeline = DwiConnectome(caps_directory=str(caps))

    assert pipeline.info == {
        "id": "aramislab/dwi-connectome",
        "author": ["Jeremy Guillon", "Alexandre Routier", "Thomas Jacquemont"],
        "version": "0.1.0",
        "space_caps": "1G",
        "space_wd": "2G",
        "dependencies": [
            {"type": "software", "name": "freesurfer", "version": ""},
            {"type": "software", "name": "fsl", "version": ">=5.0.9"},
            {"type": "software", "name": "mrtrix", "version": ""},
        ],
    }


def test_dwi_connectome_dependencies(tmp_path, mocker):
    from packaging.specifiers import SpecifierSet
    from packaging.version import Version

    from clinica.pipelines.dwi.connectome.pipeline import DwiConnectome
    from clinica.utils.check_dependency import SoftwareDependency, ThirdPartySoftware

    mocker.patch(
        "clinica.utils.check_dependency._get_freesurfer_version",
        return_value=Version("7.2.1"),
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
        {"pipelines": ["dwi-connectome"], "subjects": {"sub-01": ["ses-M000"]}},
    )
    pipeline = DwiConnectome(caps_directory=str(caps))

    assert pipeline.dependencies == [
        SoftwareDependency(
            ThirdPartySoftware.FREESURFER, SpecifierSet(">=0.0.0"), Version("7.2.1")
        ),
        SoftwareDependency(
            ThirdPartySoftware.FSL, SpecifierSet(">=5.0.9"), Version("6.0.5.2")
        ),
        SoftwareDependency(
            ThirdPartySoftware.MRTRIX, SpecifierSet(">=0.0.0"), Version("3.0.3")
        ),
    ]
