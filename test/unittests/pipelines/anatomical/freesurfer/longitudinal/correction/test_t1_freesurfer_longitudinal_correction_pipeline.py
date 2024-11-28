from clinica.utils.testing_utils import build_caps_directory


def test_t1_freesurfer_longitudinal_correction_info_loading(tmp_path):
    from clinica.pipelines.anatomical.freesurfer.longitudinal.correction.pipeline import (
        T1FreeSurferLongitudinalCorrection,
    )

    caps = build_caps_directory(
        tmp_path / "caps",
        {"pipelines": {"compute-atlas": {}}, "subjects": {"sub-01": ["ses-M000"]}},
    )
    pipeline = T1FreeSurferLongitudinalCorrection(caps_directory=str(caps))

    assert pipeline.info == {
        "id": "aramislab/t1-freesurfer-longitudinal",
        "author": ["Alexis Guyot", "Alexandre Routier"],
        "version": "0.1",
        "space_caps": "300M",
        "space_wd": "300M",
        "dependencies": [
            {"type": "software", "name": "freesurfer", "version": ">=6.0"}
        ],
    }


def test_t1_freesurfer_longitudinal_correction_dependencies(tmp_path, mocker):
    from packaging.specifiers import SpecifierSet
    from packaging.version import Version

    from clinica.pipelines.anatomical.freesurfer.longitudinal.correction.pipeline import (
        T1FreeSurferLongitudinalCorrection,
    )
    from clinica.utils.check_dependency import SoftwareDependency, ThirdPartySoftware

    mocker.patch(
        "clinica.utils.check_dependency._get_freesurfer_version",
        return_value=Version("7.2.1"),
    )
    caps = build_caps_directory(
        tmp_path / "caps",
        {
            "pipelines": {"t1-freesurfer-longitudinal-correction": {}},
            "subjects": {"sub-01": ["ses-M000"]},
        },
    )
    pipeline = T1FreeSurferLongitudinalCorrection(caps_directory=str(caps))

    assert pipeline.dependencies == [
        SoftwareDependency(
            ThirdPartySoftware.FREESURFER, SpecifierSet(">=6.0.0"), Version("7.2.1")
        )
    ]
