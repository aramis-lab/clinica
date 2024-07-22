import operator
from pathlib import Path

import pandas as pd
import pytest

truth = operator.truth
not_truth = operator.not_


def test_get_atlas_name(tmp_path):
    from clinica.iotools.utils.pipeline_handling import _get_atlas_name

    with pytest.raises(
        ValueError,
        match="Not supported pipeline foo",
    ):
        _get_atlas_name(tmp_path, "foo")
    atlas_path = (
        tmp_path
        / "t1"
        / "spm"
        / "dartel"
        / "group-foo"
        / "atlas_statistics"
        / "sub-01_ses-M000_T1w_space-AAL2_map-graymatter_statistics.tsv"
    )
    assert _get_atlas_name(atlas_path, "t1-volume") == "AAL2"
    with pytest.raises(
        ValueError,
        match="Unable to infer the atlas name",
    ):
        _get_atlas_name(atlas_path, "dwi_dti")


def test_get_mod_path_errors(tmp_path):
    from clinica.iotools.utils.pipeline_handling import _get_mod_path

    with pytest.raises(
        ValueError,
        match="Not supported pipeline foo",
    ):
        _get_mod_path(tmp_path, "foo")


@pytest.mark.parametrize(
    "pipeline,expected_path",
    [
        ("dwi_dti", ["dwi", "dti_based_processing", "atlas_statistics"]),
        ("t1-freesurfer", ["t1", "freesurfer_cross_sectional", "regional_measures"]),
        ("t1-volume", ["t1", "spm", "dartel"]),
        ("pet-volume", ["pet", "preprocessing"]),
    ],
)
def test_get_mod_path(tmp_path, pipeline, expected_path):
    from functools import reduce

    from clinica.iotools.utils.pipeline_handling import _get_mod_path

    assert _get_mod_path(tmp_path, pipeline) == reduce(
        lambda x, y: x / y, [tmp_path] + expected_path
    )


def test_get_mod_longitudinal(tmp_path):
    from clinica.iotools.utils.pipeline_handling import _get_mod_path

    assert _get_mod_path(tmp_path, "t1_freesurfer_longitudinal") is None
    (tmp_path / "t1").mkdir()
    (tmp_path / "t1" / "longfoo").mkdir()
    assert (
        _get_mod_path(tmp_path, "t1_freesurfer_longitudinal")
        == tmp_path / "t1" / "longfoo" / "freesurfer_longitudinal" / "regional_measures"
    )


def test_get_label_list(tmp_path):
    from clinica.iotools.utils.pipeline_handling import _get_label_list

    with pytest.raises(
        ValueError,
        match="Not supported pipeline bar",
    ):
        _get_label_list(tmp_path, "foo", "bar", "baz")


@pytest.mark.parametrize(
    "pipeline", ["foo", "t1-freesurfer_longitudinal", "t1-volume", "pet-volume"]
)
def test_skip_atlas_default(tmp_path, pipeline):
    from clinica.iotools.utils.pipeline_handling import _skip_atlas

    assert not _skip_atlas(tmp_path, pipeline)


@pytest.mark.parametrize("txt", ["-wm_", "-ba_"])
def test_skip_atlas_longitudinal(tmp_path, txt):
    from clinica.iotools.utils.pipeline_handling import _skip_atlas

    (tmp_path / f"foo{txt}bar.tsv").touch()
    assert _skip_atlas(tmp_path / f"foo{txt}bar.tsv", "t1-freesurfer_longitudinal")


@pytest.fixture
def expected_operator_pvc(atlas_path, pvc_restriction):
    if atlas_path == "stats.tsv":
        if pvc_restriction:
            return truth
        return not_truth
    if pvc_restriction:
        return not_truth
    return truth


@pytest.mark.parametrize("pipeline", ["t1-volume", "pet-volume"])
@pytest.mark.parametrize("atlas_path", ["stats.tsv", "pvc-rbv_stats.tsv"])
@pytest.mark.parametrize("pvc_restriction", [True, False])
def test_skip_atlas_volume_pvc(
    tmp_path, pipeline, atlas_path, pvc_restriction, expected_operator_pvc
):
    from clinica.iotools.utils.pipeline_handling import _skip_atlas

    assert expected_operator_pvc(
        _skip_atlas(tmp_path / atlas_path, pipeline, pvc_restriction=pvc_restriction)
    )


@pytest.mark.parametrize("pipeline", ["t1-volume", "pet-volume"])
@pytest.mark.parametrize("tracers", [["fdg"], ["trc1", "fdg", "trc2"]])
def test_skip_atlas_volume_tracers(tmp_path, pipeline, tracers):
    from clinica.iotools.utils.pipeline_handling import _skip_atlas

    assert _skip_atlas(
        tmp_path / "pvc-rbv_stats.tsv", pipeline, tracers_selection=tracers
    )


@pytest.fixture
def expected_operator_volume(atlas_path, pvc_restriction):
    if pvc_restriction:
        if atlas_path == "trc-fdg_stats.tsv":
            return truth
        return not_truth
    if pvc_restriction is False and atlas_path == "trc-fdg_pvc-rbv_stats.tsv":
        return truth
    return not_truth


@pytest.mark.parametrize("pipeline", ["t1-volume", "pet-volume"])
@pytest.mark.parametrize(
    "atlas_path", ["trc-fdg_pvc-rbv_stats.tsv", "trc-fdg_stats.tsv"]
)
@pytest.mark.parametrize("tracers", [["fdg"], ["trc2", "fdg"]])
@pytest.mark.parametrize("pvc_restriction", [None, True, False])
def test_skip_atlas_volume(
    tmp_path, pipeline, atlas_path, tracers, pvc_restriction, expected_operator_volume
):
    from clinica.iotools.utils.pipeline_handling import _skip_atlas

    assert expected_operator_volume(
        _skip_atlas(
            tmp_path / atlas_path,
            pipeline,
            pvc_restriction=pvc_restriction,
            tracers_selection=tracers,
        )
    )


def test_extract_metrics_from_pipeline_errors(tmp_path):
    from clinica.iotools.utils.pipeline_handling import _extract_metrics_from_pipeline

    with pytest.raises(
        KeyError,
        match="Fields `participant_id` and `session_id` are required.",
    ):
        _extract_metrics_from_pipeline(tmp_path, pd.DataFrame(), ["metrics"], "foo")
    df = pd.DataFrame([["bar", "bar"]], columns=["participant_id", "session_id"])
    (tmp_path / "groups").mkdir()
    (tmp_path / "groups" / "UnitTest").mkdir()
    with pytest.raises(
        ValueError,
        match="Not supported pipeline foo",
    ):
        _extract_metrics_from_pipeline(tmp_path, df, ["metrics"], "foo")


def test_extract_metrics_from_pipeline(tmp_path):
    from clinica.iotools.utils.pipeline_handling import _extract_metrics_from_pipeline

    df = pd.DataFrame([["bar", "bar"]], columns=["participant_id", "session_id"])
    assert _extract_metrics_from_pipeline(tmp_path, df, ["metrics"], "foo") == (
        df,
        None,
    )
    (tmp_path / "groups").mkdir()
    (tmp_path / "groups" / "UnitTest").mkdir()
    x, y = _extract_metrics_from_pipeline(tmp_path, df, ["metrics"], "t1-volume")
    assert isinstance(x, pd.DataFrame)
    assert isinstance(y, pd.DataFrame)
    assert len(x) == 1
    assert len(y) == 0
    assert set(x.columns) == {"participant_id", "session_id"}
    assert set(y.columns) == {
        "atlas_id",
        "pvc",
        "last_column_name",
        "tracer",
        "pipeline_name",
        "first_column_name",
        "group_id",
        "regions_number",
    }


def test_t1_freesurfer_pipeline_nothing_found(tmp_path):
    from pandas.testing import assert_frame_equal

    from clinica.iotools.utils.pipeline_handling import t1_freesurfer_pipeline

    caps = tmp_path / "caps"
    caps.mkdir()
    merged_df = pd.DataFrame(
        {
            "participant_id": ["sub-01"],
            "session_id": ["ses-M000"],
            "age": [85],
        }
    )
    merged_df.set_index(
        ["participant_id", "session_id"], inplace=True, verify_integrity=True
    )
    merged_df_with_t1_freesurfer_metrics, summary = t1_freesurfer_pipeline(
        caps, merged_df
    )
    merged_df_with_t1_freesurfer_metrics.set_index(
        ["participant_id", "session_id"], inplace=True, verify_integrity=True
    )

    assert_frame_equal(merged_df_with_t1_freesurfer_metrics, merged_df)
    assert summary.empty


def test_t1_freesurfer_longitudinal_pipeline_nothing_found(tmp_path):
    from pandas.testing import assert_frame_equal

    from clinica.iotools.utils.pipeline_handling import (
        t1_freesurfer_longitudinal_pipeline,
    )

    caps = tmp_path / "caps"
    caps.mkdir()
    merged_df = pd.DataFrame(
        {
            "participant_id": ["sub-01"],
            "session_id": ["ses-M000"],
            "age": [85],
        }
    )
    merged_df.set_index(
        ["participant_id", "session_id"], inplace=True, verify_integrity=True
    )
    merged_df_with_t1_freesurfer_metrics, summary = t1_freesurfer_longitudinal_pipeline(
        caps, merged_df
    )
    merged_df_with_t1_freesurfer_metrics.set_index(
        ["participant_id", "session_id"], inplace=True, verify_integrity=True
    )

    assert_frame_equal(merged_df_with_t1_freesurfer_metrics, merged_df)
    assert summary.empty


def write_fake_statistics(tsv_file: Path):
    fake = pd.DataFrame(
        {
            "label_name": ["foo", "bar", "baz"],
            "label_value": [1.2, 2.3, 3.4],
        }
    )
    fake.to_csv(tsv_file, sep="\t")


def test_t1_freesurfer_pipeline(tmp_path):
    from clinica.iotools.utils.pipeline_handling import t1_freesurfer_pipeline

    caps = tmp_path / "caps"
    regional_measures_folder = (
        caps
        / "subjects"
        / "sub-01"
        / "ses-M000"
        / "t1"
        / "freesurfer_cross_sectional"
        / "regional_measures"
    )
    regional_measures_folder.mkdir(parents=True)
    for file in (
        "sub-01_ses-M000_parcellation-desikan_thickness.tsv",
        "sub-01_ses-M000_parcellation-desikan_area.tsv",
        "sub-01_ses-M000_parcellation-destrieux_thickness.tsv",
        "sub-01_ses-M000_segmentationVolumes.tsv",
    ):
        write_fake_statistics(regional_measures_folder / file)

    assert len([f for f in regional_measures_folder.iterdir()]) == 4
    merged_df = pd.DataFrame(
        {
            "participant_id": ["sub-01"],
            "session_id": ["ses-M000"],
            "age": [85],
        }
    )
    merged_df.set_index(
        ["participant_id", "session_id"], inplace=True, verify_integrity=True
    )
    merged_df_with_t1_freesurfer_metrics, summary = t1_freesurfer_pipeline(
        caps, merged_df
    )
    merged_df_with_t1_freesurfer_metrics.set_index(
        ["participant_id", "session_id"], inplace=True, verify_integrity=True
    )

    assert len(merged_df_with_t1_freesurfer_metrics) == 1
    assert set(merged_df_with_t1_freesurfer_metrics.columns) == {
        "age",
        "t1-freesurfer_atlas-desikan_ROI-bar_thickness",
        "t1-freesurfer_atlas-desikan_ROI-baz_thickness",
        "t1-freesurfer_atlas-desikan_ROI-foo_thickness",
        "t1-freesurfer_atlas-destrieux_ROI-bar_thickness",
        "t1-freesurfer_atlas-destrieux_ROI-baz_thickness",
        "t1-freesurfer_atlas-destrieux_ROI-foo_thickness",
        "t1-freesurfer_segmentation-volumes_ROI-bar_volume",
        "t1-freesurfer_segmentation-volumes_ROI-baz_volume",
        "t1-freesurfer_segmentation-volumes_ROI-foo_volume",
    }
    assert summary.pipeline_name.to_list() == ["t1-freesurfer"] * 3
    assert summary.atlas_id.to_list() == ["desikan", "destrieux", "volumes"]
    assert summary.regions_number.to_list() == [3, 3, 3]
    assert summary.first_column_name.to_list() == [
        "t1-freesurfer_atlas-desikan_ROI-foo_thickness",
        "t1-freesurfer_atlas-destrieux_ROI-foo_thickness",
        "t1-freesurfer_segmentation-volumes_ROI-foo_volume",
    ]
    assert summary.last_column_name.to_list() == [
        "t1-freesurfer_atlas-desikan_ROI-baz_thickness",
        "t1-freesurfer_atlas-destrieux_ROI-baz_thickness",
        "t1-freesurfer_segmentation-volumes_ROI-baz_volume",
    ]
