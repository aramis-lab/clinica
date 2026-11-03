from pathlib import Path

import pytest

from clinica.pipelines.anatomical.freesurfer.t1.pipeline import T1FreeSurfer


@pytest.mark.parametrize("caps_name", ["caps", "sub-other_ses-M999"])
def test_get_processed_images_uses_caps_subject_and_session(
    tmp_path, monkeypatch, caps_name
):
    monkeypatch.chdir(tmp_path)
    caps = Path(caps_name)
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

    assert T1FreeSurfer.get_processed_images(
        caps, ["sub-01", "sub-02"], ["ses-M000", "ses-M000"]
    ) == ["sub-01_ses-M000"]
