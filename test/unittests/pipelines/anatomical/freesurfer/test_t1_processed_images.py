from pathlib import Path

import pytest

from clinica.pipelines.anatomical.freesurfer.t1.pipeline import T1FreeSurfer


@pytest.mark.parametrize("caps_name", ["caps", "sub-other_ses-M999"])
def test_get_processed_images_reports_conflicting_ids(tmp_path, monkeypatch, caps_name):
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

    subjects = ["sub-01", "sub-02"]
    sessions = ["ses-M000", "ses-M000"]
    if caps_name == "caps":
        assert T1FreeSurfer.get_processed_images(caps, subjects, sessions) == [
            "sub-01_ses-M000"
        ]
    else:
        with pytest.raises(ValueError, match="more than one distinct"):
            T1FreeSurfer.get_processed_images(caps, subjects, sessions)
