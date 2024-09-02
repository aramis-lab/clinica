from clinica.utils.testing_utils import build_caps_directory


def test_t1_freesurfer_atlas_info_loading(tmp_path):
    from clinica.pipelines.anatomical.freesurfer.atlas.pipeline import T1FreeSurferAtlas

    caps = build_caps_directory(
        tmp_path / "caps",
        {"pipelines": ["compute-atlas"], "subjects": {"sub-01": ["ses-M000"]}},
    )
    pipeline = T1FreeSurferAtlas(caps_directory=str(caps))

    assert pipeline.info == {
        "id": "aramislab/compute-atlas",
        "author": ["Matthieu Joulot"],
        "version": "0.1",
        "space_caps": "1G",
        "space_wd": "1G",
        "dependencies": [
            {"type": "software", "name": "freesurfer", "version": ">=6.0.0"}
        ],
    }


def test_t1_freesurfer_atlas_dependencies(tmp_path, mocker):
    from packaging.specifiers import SpecifierSet
    from packaging.version import Version

    from clinica.pipelines.anatomical.freesurfer.atlas.pipeline import T1FreeSurferAtlas
    from clinica.utils.check_dependency import SoftwareDependency, ThirdPartySoftware

    mocker.patch(
        "clinica.utils.check_dependency._get_freesurfer_version",
        return_value=Version("7.2.1"),
    )
    caps = build_caps_directory(
        tmp_path / "caps",
        {"pipelines": ["compute-atlas"], "subjects": {"sub-01": ["ses-M000"]}},
    )
    pipeline = T1FreeSurferAtlas(caps_directory=str(caps))

    assert pipeline.dependencies == [
        SoftwareDependency(
            ThirdPartySoftware.FREESURFER, SpecifierSet(">=6.0.0"), Version("7.2.1")
        )
    ]
