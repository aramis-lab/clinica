from unittest import mock

from clinica.utils.testing_utils import build_caps_directory


def filter_qc_mock() -> tuple[list[str], list[str]]:
    return [], []


def test_dwi_preprocessing_using_t1_info_loading(tmp_path):
    from clinica.pipelines.dwi.preprocessing.t1.pipeline import DwiPreprocessingUsingT1

    caps = build_caps_directory(
        tmp_path / "caps",
        {
            "pipelines": ["dwi-preprocessing-using-t1"],
            "subjects": {"sub-01": ["ses-M000"]},
        },
    )
    with mock.patch(
        "clinica.pipelines.dwi.preprocessing.t1.pipeline.DwiPreprocessingUsingT1.filter_qc",
        wraps=filter_qc_mock,
    ) as qc_filter_mock:
        pipeline = DwiPreprocessingUsingT1(caps_directory=str(caps))
        qc_filter_mock.assert_called_once()

    assert pipeline.info == {
        "id": "aramislab/dwi-preprocessing-using-t1",
        "author": ["Junhao Wen", "Thomas Jacquemont", "Alexandre Routier"],
        "credits": ["Nipype"],
        "space_caps": "150M",
        "space_wd": "3G",
        "version": "0.2.0",
        "dependencies": [
            {"type": "software", "name": "fsl", "version": ">=6.0.1"},
            {"type": "software", "name": "convert3d", "version": ""},
            {"type": "software", "name": "ants", "version": ""},
            {"type": "software", "name": "mrtrix", "version": ""},
        ],
    }


def test_dwi_preprocessing_using_t1_dependencies(tmp_path, mocker):
    from packaging.specifiers import SpecifierSet
    from packaging.version import Version

    from clinica.pipelines.dwi.preprocessing.t1.pipeline import DwiPreprocessingUsingT1
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
    mocker.patch(
        "clinica.utils.check_dependency._get_convert3d_version",
        return_value=Version("1.0.0"),
    )
    caps = build_caps_directory(
        tmp_path / "caps",
        {
            "pipelines": ["dwi-preprocessing-using-t1"],
            "subjects": {"sub-01": ["ses-M000"]},
        },
    )
    with mock.patch(
        "clinica.pipelines.dwi.preprocessing.t1.pipeline.DwiPreprocessingUsingT1.filter_qc",
        wraps=filter_qc_mock,
    ) as qc_filter_mock:
        pipeline = DwiPreprocessingUsingT1(caps_directory=str(caps))
        qc_filter_mock.assert_called_once()

    assert pipeline.dependencies == [
        SoftwareDependency(
            ThirdPartySoftware.FSL, SpecifierSet(">=6.0.1"), Version("6.0.5.2")
        ),
        SoftwareDependency(
            ThirdPartySoftware.CONVERT3D, SpecifierSet(">=0.0.0"), Version("1.0.0")
        ),
        SoftwareDependency(
            ThirdPartySoftware.ANTS, SpecifierSet(">=0.0.0"), Version("2.3.2")
        ),
        SoftwareDependency(
            ThirdPartySoftware.MRTRIX, SpecifierSet(">=0.0.0"), Version("3.0.3")
        ),
    ]
