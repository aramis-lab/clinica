import operator
from pathlib import Path

import pandas as pd
import pytest

from clinica.iotools.utils.pipeline_handling import PipelineNameForMetricExtraction

truth = operator.truth
not_truth = operator.not_


@pytest.fixture
def atlas_path(tmp_path) -> Path:
    filename = "sub-01_ses-M000_T1w_space-AAL2_map-graymatter_statistics.tsv"
    return (
        tmp_path / "t1" / "spm" / "dartel" / "group-foo" / "atlas_statistics" / filename
    )


def test_get_atlas_name(atlas_path):
    from clinica.iotools.utils.pipeline_handling import _get_atlas_name

    assert (
        _get_atlas_name(atlas_path, PipelineNameForMetricExtraction.T1_VOLUME) == "AAL2"
    )


def test_get_atlas_name_error(atlas_path):
    from clinica.iotools.utils.pipeline_handling import _get_atlas_name

    with pytest.raises(
        ValueError,
        match="Unable to infer the atlas name",
    ):
        _get_atlas_name(atlas_path, PipelineNameForMetricExtraction.DWI_DTI)


@pytest.mark.parametrize(
    "pipeline,expected_path",
    [
        (
            PipelineNameForMetricExtraction.DWI_DTI,
            ["dwi", "dti_based_processing", "atlas_statistics"],
        ),
        (
            PipelineNameForMetricExtraction.T1_FREESURFER,
            ["t1", "freesurfer_cross_sectional", "regional_measures"],
        ),
        (PipelineNameForMetricExtraction.T1_VOLUME, ["t1", "spm", "dartel"]),
        (PipelineNameForMetricExtraction.PET_VOLUME, ["pet", "preprocessing"]),
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

    assert (
        _get_mod_path(tmp_path, PipelineNameForMetricExtraction.T1_FREESURFER_LONGI)
        is None
    )
    (tmp_path / "t1").mkdir()
    (tmp_path / "t1" / "longfoo").mkdir()
    assert (
        _get_mod_path(tmp_path, PipelineNameForMetricExtraction.T1_FREESURFER_LONGI)
        == tmp_path / "t1" / "longfoo" / "freesurfer_longitudinal" / "regional_measures"
    )


@pytest.mark.parametrize(
    "pipeline",
    [
        PipelineNameForMetricExtraction.T1_FREESURFER_LONGI,
        PipelineNameForMetricExtraction.T1_VOLUME,
        PipelineNameForMetricExtraction.PET_VOLUME,
    ],
)
def test_skip_atlas_default(tmp_path, pipeline):
    from clinica.iotools.utils.pipeline_handling import _skip_atlas

    assert not _skip_atlas(tmp_path, pipeline)


@pytest.mark.parametrize("txt", ["-wm_", "-ba_"])
def test_skip_atlas_longitudinal(tmp_path, txt):
    from clinica.iotools.utils.pipeline_handling import _skip_atlas

    (tmp_path / f"foo{txt}bar.tsv").touch()
    assert _skip_atlas(
        tmp_path / f"foo{txt}bar.tsv",
        PipelineNameForMetricExtraction.T1_FREESURFER_LONGI,
    )


@pytest.fixture
def expected_operator_pvc(atlas_path, pvc_restriction):
    if atlas_path == "stats.tsv":
        if pvc_restriction:
            return truth
        return not_truth
    if pvc_restriction:
        return not_truth
    return truth


@pytest.mark.parametrize(
    "pipeline",
    [
        PipelineNameForMetricExtraction.T1_VOLUME,
        PipelineNameForMetricExtraction.PET_VOLUME,
    ],
)
@pytest.mark.parametrize("atlas_path", ["stats.tsv", "pvc-rbv_stats.tsv"])
@pytest.mark.parametrize("pvc_restriction", [True, False])
def test_skip_atlas_volume_pvc(
    tmp_path, pipeline, atlas_path, pvc_restriction, expected_operator_pvc
):
    from clinica.iotools.utils.pipeline_handling import _skip_atlas

    assert expected_operator_pvc(
        _skip_atlas(tmp_path / atlas_path, pipeline, pvc_restriction=pvc_restriction)
    )


@pytest.mark.parametrize(
    "pipeline",
    [
        PipelineNameForMetricExtraction.T1_VOLUME,
        PipelineNameForMetricExtraction.PET_VOLUME,
    ],
)
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


@pytest.mark.parametrize(
    "pipeline",
    [
        PipelineNameForMetricExtraction.T1_VOLUME,
        PipelineNameForMetricExtraction.PET_VOLUME,
    ],
)
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
        _extract_metrics_from_pipeline(
            tmp_path,
            pd.DataFrame(),
            ["metrics"],
            PipelineNameForMetricExtraction.DWI_DTI,
        )


def test_extract_metrics_from_pipeline(tmp_path):
    from clinica.iotools.utils.pipeline_handling import _extract_metrics_from_pipeline

    df = pd.DataFrame([["bar", "bar"]], columns=["participant_id", "session_id"])
    assert _extract_metrics_from_pipeline(
        tmp_path, df, ["metrics"], PipelineNameForMetricExtraction.T1_VOLUME
    ) == (
        df,
        None,
    )
    (tmp_path / "groups").mkdir()
    (tmp_path / "groups" / "UnitTest").mkdir()
    x, y = _extract_metrics_from_pipeline(
        tmp_path, df, ["metrics"], PipelineNameForMetricExtraction.T1_VOLUME
    )
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
