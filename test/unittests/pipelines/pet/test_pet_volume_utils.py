import pytest


def test_compute_atlas_statistics(tmp_path, mocker):
    from clinica.pipelines.pet.volume.utils import compute_atlas_statistics

    atlases = ["foo", "bar", "foobar"]

    mocker.patch("clinica.utils.statistics.statistics_on_atlas")
    result = compute_atlas_statistics(tmp_path, atlases)
    assert len(result) == len(atlases)
