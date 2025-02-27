import operator
from pathlib import Path

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from clinica.iotools.pipeline_handling import PipelineNameForMetricExtraction

truth = operator.truth
not_truth = operator.not_


@pytest.fixture
def atlas_path(tmp_path, pipeline: PipelineNameForMetricExtraction) -> Path:
    if pipeline in (
        PipelineNameForMetricExtraction.T1_VOLUME,
        PipelineNameForMetricExtraction.PET_VOLUME,
    ):
        filename = "sub-01_ses-M000_T1w_space-AAL2_map-graymatter_statistics.tsv"
        return (
            tmp_path
            / "t1"
            / "spm"
            / "dartel"
            / "group-foo"
            / "atlas_statistics"
            / filename
        )
    if pipeline == PipelineNameForMetricExtraction.DWI_DTI:
        filename = "sub-01_ses-M000_dwi_space-foobar.tsv"
        return tmp_path / "dwi" / filename
    if pipeline in (
        PipelineNameForMetricExtraction.T1_FREESURFER_LONGI,
        PipelineNameForMetricExtraction.T1_FREESURFER,
    ):
        filename = "sub-01_ses-M000_parcellation-foobaz.tsv"
        return tmp_path / "freesurfer" / filename


@pytest.mark.parametrize(
    "pipeline,name",
    [
        (PipelineNameForMetricExtraction.T1_VOLUME, "AAL2"),
        (PipelineNameForMetricExtraction.PET_VOLUME, "AAL2"),
        (PipelineNameForMetricExtraction.DWI_DTI, "foobar"),
        (PipelineNameForMetricExtraction.T1_FREESURFER_LONGI, "foobaz"),
        (PipelineNameForMetricExtraction.T1_FREESURFER, "foobaz"),
    ],
)
def test_get_atlas_name(atlas_path, pipeline, name):
    from clinica.iotools.pipeline_handling import _get_atlas_name

    assert _get_atlas_name(atlas_path, pipeline) == f"atlas-{name}"


def _get_wrong_atlas_path(
    folder: Path, pipeline: PipelineNameForMetricExtraction
) -> Path:
    if pipeline in (
        PipelineNameForMetricExtraction.T1_VOLUME,
        PipelineNameForMetricExtraction.PET_VOLUME,
    ):
        filename = (
            "sub-01_ses-M000_T1w_spaaaacccceee-AAL2_map-graymatter_statistics.tsv"
        )
        return (
            folder
            / "t1"
            / "spm"
            / "dartel"
            / "group-foo"
            / "atlas_statistics"
            / filename
        )
    if pipeline == PipelineNameForMetricExtraction.DWI_DTI:
        filename = "sub-01_ses-M000_dwiii_space-foobar.tsv"
        return folder / "dwi" / filename


@pytest.mark.parametrize(
    "pipeline",
    [
        PipelineNameForMetricExtraction.T1_VOLUME,
        PipelineNameForMetricExtraction.PET_VOLUME,
        PipelineNameForMetricExtraction.DWI_DTI,
    ],
)
def test_get_atlas_name_error(tmp_path, pipeline):
    from clinica.iotools.pipeline_handling import _get_atlas_name

    with pytest.raises(
        ValueError,
        match="Unable to infer the atlas name",
    ):
        _get_atlas_name(
            _get_wrong_atlas_path(tmp_path, pipeline),
            PipelineNameForMetricExtraction.DWI_DTI,
        )


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
def test_get_modality_path(tmp_path, pipeline, expected_path):
    from functools import reduce

    from clinica.iotools.pipeline_handling import _get_modality_path

    assert _get_modality_path(tmp_path, pipeline) == reduce(
        lambda x, y: x / y, [tmp_path] + expected_path
    )


def test_get_modality_path_longitudinal(tmp_path):
    from clinica.iotools.pipeline_handling import _get_modality_path

    assert (
        _get_modality_path(
            tmp_path, PipelineNameForMetricExtraction.T1_FREESURFER_LONGI
        )
        is None
    )
    (tmp_path / "t1").mkdir()
    (tmp_path / "t1" / "longfoo").mkdir()

    assert (
        _get_modality_path(
            tmp_path, PipelineNameForMetricExtraction.T1_FREESURFER_LONGI
        )
        == tmp_path / "t1" / "longfoo" / "freesurfer_longitudinal" / "regional_measures"
    )


@pytest.mark.parametrize("pipeline", PipelineNameForMetricExtraction)
def test_skip_atlas_non_existing_file(tmp_path, pipeline):
    from clinica.iotools.pipeline_handling import _skip_atlas

    assert _skip_atlas(tmp_path / "foo.tsv", pipeline)


@pytest.mark.parametrize("txt", ["-wm_", "-ba_"])
def test_skip_atlas_longitudinal(tmp_path, txt: str):
    from clinica.iotools.pipeline_handling import _skip_atlas

    (tmp_path / f"foo{txt}bar.tsv").touch()

    assert _skip_atlas(
        tmp_path / f"foo{txt}bar.tsv",
        PipelineNameForMetricExtraction.T1_FREESURFER_LONGI,
    )


@pytest.fixture
def expected_operator_pvc(atlas_path, pvc_restriction):
    if atlas_path == "sub-01_ses-M000_space-foobar_stats.tsv":
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
@pytest.mark.parametrize(
    "atlas_path",
    ["sub-01_ses-M000_space-foobar_stats.tsv", "pvc-rbv_space-foo_stats.tsv"],
)
@pytest.mark.parametrize("pvc_restriction", [True, False])
def test_skip_atlas_volume_pvc(
    tmp_path, pipeline, atlas_path: str, pvc_restriction: bool, expected_operator_pvc
):
    from clinica.iotools.pipeline_handling import _skip_atlas

    (tmp_path / atlas_path).touch()
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
def test_skip_atlas_volume_tracers(
    tmp_path, pipeline: PipelineNameForMetricExtraction, tracers: list[str]
):
    from clinica.iotools.pipeline_handling import _skip_atlas

    assert _skip_atlas(
        tmp_path / "pvc-rbv_stats.tsv", pipeline, tracers_selection=tracers
    )


@pytest.fixture
def expected_operator_volume(atlas_path: str, pvc_restriction: bool):
    if pvc_restriction:
        if atlas_path == "sub-01_ses-M000_space-bar_trc-fdg_stats.tsv":
            return truth
        return not_truth
    if (
        pvc_restriction is False
        and atlas_path == "sub-01_ses-M000_space-foo_trc-fdg_pvc-rbv_stats.tsv"
    ):
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
    "atlas_path",
    [
        "sub-01_ses-M000_space-foo_trc-fdg_pvc-rbv_stats.tsv",
        "sub-01_ses-M000_space-bar_trc-fdg_stats.tsv",
    ],
)
@pytest.mark.parametrize("tracers", [["fdg"], ["trc2", "fdg"]])
@pytest.mark.parametrize("pvc_restriction", [None, True, False])
def test_skip_atlas_volume(
    tmp_path,
    pipeline: PipelineNameForMetricExtraction,
    atlas_path: str,
    tracers: list[str],
    pvc_restriction: bool,
    expected_operator_volume,
):
    from clinica.iotools.pipeline_handling import _skip_atlas

    (tmp_path / atlas_path).touch()

    assert expected_operator_volume(
        _skip_atlas(
            tmp_path / atlas_path,
            pipeline,
            pvc_restriction=pvc_restriction,
            tracers_selection=tracers,
        )
    )


def test_extract_metrics_from_pipeline_errors(tmp_path):
    from clinica.iotools.pipeline_handling import _extract_metrics_from_pipeline

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


def test_extract_metrics_from_pipeline_empty(tmp_path):
    from clinica.iotools.pipeline_handling import _extract_metrics_from_pipeline

    (tmp_path / "groups").mkdir()
    df = pd.DataFrame([["bar", "bar"]], columns=["participant_id", "session_id"])
    x, y = _extract_metrics_from_pipeline(
        tmp_path, df, ["metrics"], PipelineNameForMetricExtraction.T1_VOLUME
    )

    assert_frame_equal(
        x, pd.DataFrame([["bar", "bar"]], columns=["participant_id", "session_id"])
    )
    assert y.empty


def test_extract_metrics_from_pipeline(tmp_path):
    from clinica.iotools.pipeline_handling import _extract_metrics_from_pipeline

    (tmp_path / "groups" / "UnitTest").mkdir(parents=True)
    df = pd.DataFrame([["bar", "bar"]], columns=["participant_id", "session_id"])
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


def test_extract_metrics_from_t1_freesurfer_nothing_found(tmp_path):
    from pandas.testing import assert_frame_equal

    from clinica.iotools.pipeline_handling import pipeline_metric_extractor_factory

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
    merged_df_with_t1_freesurfer_metrics, summary = pipeline_metric_extractor_factory(
        PipelineNameForMetricExtraction.T1_FREESURFER
    )(caps, merged_df)
    merged_df_with_t1_freesurfer_metrics.set_index(
        ["participant_id", "session_id"], inplace=True, verify_integrity=True
    )

    assert_frame_equal(merged_df_with_t1_freesurfer_metrics, merged_df)
    assert summary.empty


def test_extract_metrics_from_t1_freesurfer_longitudinal_nothing_found(tmp_path):
    from pandas.testing import assert_frame_equal

    from clinica.iotools.pipeline_handling import pipeline_metric_extractor_factory

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
    merged_df_with_t1_freesurfer_metrics, summary = pipeline_metric_extractor_factory(
        PipelineNameForMetricExtraction.T1_FREESURFER_LONGI
    )(caps, merged_df)
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


def test_extract_metrics_from_t1_freesurfer(tmp_path):
    from clinica.iotools.pipeline_handling import pipeline_metric_extractor_factory

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
    merged_df_with_t1_freesurfer_metrics, summary = pipeline_metric_extractor_factory(
        PipelineNameForMetricExtraction.T1_FREESURFER
    )(caps, merged_df)
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
