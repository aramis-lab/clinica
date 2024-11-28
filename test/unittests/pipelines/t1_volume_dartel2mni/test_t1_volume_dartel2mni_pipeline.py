from clinica.utils.testing_utils import build_caps_directory


def test_t1_volume_dartel2mni_info_loading(tmp_path):
    from clinica.pipelines.t1_volume_dartel2mni.t1_volume_dartel2mni_pipeline import (
        T1VolumeDartel2MNI,
    )

    caps = build_caps_directory(
        tmp_path / "caps",
        {
            "pipelines": {"t1-volume-dartel2mni": {}},
            "subjects": {"sub-01": ["ses-M000"]},
        },
    )
    pipeline = T1VolumeDartel2MNI(caps_directory=str(caps), group_label="test")

    assert pipeline.info == {
        "id": "aramislab/t1-volume-dartel2mni",
        "author": "Jorge Samper-Gonzalez; Simona Bottani",
        "version": "0.1.0",
        "space_caps": "30M",
        "space_wd": "150M",
        "dependencies": [{"type": "software", "name": "spm", "version": ">=12"}],
    }


def test_t1_volume_dartel2mni_dependencies(tmp_path, mocker):
    from packaging.specifiers import SpecifierSet
    from packaging.version import Version

    from clinica.pipelines.t1_volume_dartel2mni.t1_volume_dartel2mni_pipeline import (
        T1VolumeDartel2MNI,
    )
    from clinica.utils.check_dependency import SoftwareDependency, ThirdPartySoftware

    mocker.patch(
        "clinica.utils.check_dependency._get_spm_version",
        return_value=Version("12.7219"),
    )
    caps = build_caps_directory(
        tmp_path / "caps",
        {
            "pipelines": {"t1-volume-dartel2mni": {}},
            "subjects": {"sub-01": ["ses-M000"]},
        },
    )
    pipeline = T1VolumeDartel2MNI(caps_directory=str(caps), group_label="test")

    assert pipeline.dependencies == [
        SoftwareDependency(
            ThirdPartySoftware.SPM, SpecifierSet(">=12"), Version("12.7219")
        ),
    ]
