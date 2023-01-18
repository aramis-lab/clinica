import pytest


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
    assert _get_atlas_name(atlas_path, "t1_volume") == "AAL2"
    with pytest.raises(
        ValueError,
        match="Unable to infer the atlas name",
    ):
        _get_atlas_name(atlas_path, "dwi_dti")


def test_get_mod_path(tmp_path):
    from clinica.iotools.utils.pipeline_handling import _get_mod_path

    with pytest.raises(
        ValueError,
        match="Not supported pipeline foo",
    ):
        _get_mod_path(tmp_path, "foo")
    assert (
        _get_mod_path(tmp_path, "dwi_dti")
        == tmp_path / "dwi" / "dti_based_processing" / "atlas_statistics"
    )
    assert (
        _get_mod_path(tmp_path, "t1_freesurfer")
        == tmp_path / "t1" / "freesurfer_cross_sectional" / "regional_measures"
    )
    assert _get_mod_path(tmp_path, "t1_volume") == tmp_path / "t1" / "spm" / "dartel"
    assert _get_mod_path(tmp_path, "pet_volume") == tmp_path / "pet" / "preprocessing"
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


def test_skip_atlas(tmp_path):
    from clinica.iotools.utils.pipeline_handling import _skip_atlas

    for pipeline in ("foo", "t1-freesurfer_longitudinal", "t1_volume", "pet_volume"):
        assert not _skip_atlas(tmp_path, pipeline)
    for txt in ("-wm_", "-ba_"):
        (tmp_path / f"foo{txt}bar.tsv").touch()
        assert _skip_atlas(tmp_path / f"foo{txt}bar.tsv", "t1-freesurfer_longitudinal")
    atlas_path = tmp_path / "stats.tsv"
    for pipeline in ("t1_volume", "pet_volume"):
        assert _skip_atlas(atlas_path, pipeline, pvc_restriction=True)
        assert not _skip_atlas(atlas_path, pipeline, pvc_restriction=False)
    atlas_path = tmp_path / "pvc-rbv_stats.tsv"
    for pipeline in ("t1_volume", "pet_volume"):
        assert not _skip_atlas(atlas_path, pipeline, pvc_restriction=True)
        assert _skip_atlas(atlas_path, pipeline, pvc_restriction=False)
    for pipeline in ("t1_volume", "pet_volume"):
        assert _skip_atlas(atlas_path, pipeline, tracers_selection=["fdg"])
        assert _skip_atlas(
            atlas_path, pipeline, tracers_selection=["trc1", "fdg", "trc2"]
        )
    atlas_path = tmp_path / "trc-fdg_pvc-rbv_stats.tsv"
    for pipeline in ("t1_volume", "pet_volume"):
        assert not _skip_atlas(atlas_path, pipeline, tracers_selection=["fdg"])
        assert not _skip_atlas(atlas_path, pipeline, tracers_selection=["trc2", "fdg"])
        assert not _skip_atlas(
            atlas_path,
            pipeline,
            pvc_restriction=True,
            tracers_selection=["trc2", "fdg"],
        )
        assert _skip_atlas(
            atlas_path,
            pipeline,
            pvc_restriction=False,
            tracers_selection=["trc2", "fdg"],
        )
    atlas_path = tmp_path / "trc-fdg_stats.tsv"
    for pipeline in ("t1_volume", "pet_volume"):
        assert not _skip_atlas(atlas_path, pipeline, tracers_selection=["fdg"])
        assert not _skip_atlas(atlas_path, pipeline, tracers_selection=["trc2", "fdg"])
        assert _skip_atlas(
            atlas_path,
            pipeline,
            pvc_restriction=True,
            tracers_selection=["trc2", "fdg"],
        )
        assert not _skip_atlas(
            atlas_path,
            pipeline,
            pvc_restriction=False,
            tracers_selection=["trc2", "fdg"],
        )
