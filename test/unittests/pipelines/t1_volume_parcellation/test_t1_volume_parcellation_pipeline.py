from clinica.utils.testing_utils import build_caps_directory


def test_t1_volume_parcellation_info_loading(tmp_path):
    from clinica.pipelines.t1_volume_parcellation.t1_volume_parcellation_pipeline import (
        T1VolumeParcellation,
    )

    caps = build_caps_directory(
        tmp_path / "caps",
        {"pipelines": ["t1-volume-parcellation"], "subjects": {"sub-01": ["ses-M000"]}},
    )
    pipeline = T1VolumeParcellation(caps_directory=str(caps))

    assert pipeline.info == {
        "id": "aramislab/t1-volume-parcellation",
        "author": "Simona Bottani",
        "version": "0.1.0",
        "space_caps": "760K",
        "space_wd": "520K",
        "dependencies": [],
    }


def test_t1_volume_parcellation_dependencies(tmp_path):
    from clinica.pipelines.t1_volume_parcellation.t1_volume_parcellation_pipeline import (
        T1VolumeParcellation,
    )

    caps = build_caps_directory(
        tmp_path / "caps",
        {"pipelines": ["t1-volume-parcellation"], "subjects": {"sub-01": ["ses-M000"]}},
    )
    pipeline = T1VolumeParcellation(caps_directory=str(caps))

    assert pipeline.dependencies == []
