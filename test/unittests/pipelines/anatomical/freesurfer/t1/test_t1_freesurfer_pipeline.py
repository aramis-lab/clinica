from pathlib import Path

import pytest

from clinica.utils.testing_utils import build_bids_directory


def test_t1_freesurfer_info_loading(tmp_path):
    from clinica.pipelines.anatomical.freesurfer.t1.pipeline import T1FreeSurfer

    bids = build_bids_directory(tmp_path / "bids", {"sub-01": ["ses-M00"]})
    pipeline = T1FreeSurfer(bids_directory=str(bids))

    assert pipeline.info == {
        "id": "aramislab/t1-freesurfer",
        "author": ["Junhao Wen", "Alexandre Routier"],
        "version": "0.1",
        "space_caps": "1G",
        "space_wd": "1G",
        "dependencies": [
            {"type": "software", "name": "freesurfer", "version": ">=6.0.0"}
        ],
    }


def test_t1_freesurfer_dependencies(tmp_path, mocker):
    from packaging.specifiers import SpecifierSet
    from packaging.version import Version

    from clinica.pipelines.anatomical.freesurfer.t1.pipeline import T1FreeSurfer
    from clinica.utils.check_dependency import SoftwareDependency, ThirdPartySoftware

    mocker.patch(
        "clinica.utils.check_dependency._get_freesurfer_version",
        return_value=Version("7.2.1"),
    )
    bids = build_bids_directory(tmp_path / "bids", {"sub-01": ["ses-M00"]})
    pipeline = T1FreeSurfer(bids_directory=str(bids))

    assert pipeline.dependencies == [
        SoftwareDependency(
            ThirdPartySoftware.FREESURFER, SpecifierSet(">=6.0.0"), Version("7.2.1")
        )
    ]


def _build_caps_freesurfer(tmp_path: Path, caps_name: str) -> Path:
    caps = tmp_path / caps_name
    output = (
        caps
        / "subjects"
        / "sub-01"
        / "ses-M000"
        / "t1"
        / "freesurfer_cross_sectional"
        / "sub-01_ses-M000"
        / "mri"
        / "aparc.a2009s+aseg.mgz"
    )
    output.parent.mkdir(parents=True)
    output.touch()
    (caps / "dataset_description.json").write_text(
        '{"DatasetType": "derivative"}', encoding="utf-8"
    )
    return caps


def test_get_processed_images_returns_completed_image(tmp_path, monkeypatch):
    from clinica.pipelines.anatomical.freesurfer.t1.pipeline import T1FreeSurfer

    monkeypatch.chdir(tmp_path)
    # Relative paths avoid passing a Windows drive letter to the globber.
    caps = _build_caps_freesurfer(tmp_path, "caps").relative_to(tmp_path)

    assert T1FreeSurfer.get_processed_images(
        caps, ["sub-01", "sub-02"], ["ses-M000", "ses-M000"]
    ) == ["sub-01_ses-M000"]


def test_get_processed_images_rejects_conflicting_ids(tmp_path, monkeypatch):
    from clinica.pipelines.anatomical.freesurfer.t1.pipeline import T1FreeSurfer

    monkeypatch.chdir(tmp_path)
    caps = _build_caps_freesurfer(tmp_path, "sub-other_ses-M999").relative_to(tmp_path)

    with pytest.raises(ValueError, match="more than one distinct"):
        T1FreeSurfer.get_processed_images(caps, ["sub-01"], ["ses-M000"])
