import re

import pytest


def test_anat_linear_pipeline_no_input_error(tmp_path):
    from clinica.pipelines.t1_linear.anat_linear_pipeline import AnatLinear

    with pytest.raises(
        RuntimeError,
        match=(
            "The AnatLinear pipeline does not contain BIDS "
            "nor CAPS directory at the initialization."
        ),
    ):
        AnatLinear()


def test_anat_linear_pipeline_single_bids_input_error(tmp_path):
    from clinica.pipelines.t1_linear.anat_linear_pipeline import AnatLinear
    from clinica.utils.exceptions import ClinicaBIDSError

    with pytest.raises(
        ClinicaBIDSError,
        match=re.escape(
            f"The raw directory ({tmp_path}) you provided "
            "is missing a dataset_description.json file."
        ),
    ):
        AnatLinear(bids_directory=str(tmp_path))


def test_anat_linear_pipeline_single_caps_input_error(tmp_path):
    from clinica.pipelines.t1_linear.anat_linear_pipeline import AnatLinear
    from clinica.utils.exceptions import ClinicaCAPSError

    with pytest.raises(
        ClinicaCAPSError,
        match=re.escape(
            f"The derivative directory ({tmp_path}) you provided "
            "is missing a dataset_description.json file."
        ),
    ):
        AnatLinear(caps_directory=str(tmp_path))


def test_anat_linear_pipeline_write_caps_dataset_description(tmp_path):
    from clinica.pipelines.t1_linear.anat_linear_pipeline import AnatLinear
    from clinica.utils.caps import CAPSDatasetDescription, DatasetType

    bids = tmp_path / "bids"
    caps = tmp_path / "caps"
    bids.mkdir()
    caps.mkdir()
    (bids / "dataset_description.json").touch()
    (bids / "sub-01").mkdir()

    AnatLinear(bids_directory=str(bids), caps_directory=str(caps))

    files = [f for f in caps.iterdir()]
    assert len(files) == 1

    desc = CAPSDatasetDescription.from_file(caps / "dataset_description.json")

    assert desc.bids_version == "1.7.0"
    assert desc.caps_version == "1.0.0"
    assert desc.dataset_type == DatasetType.DERIVATIVE
    assert desc.processing[0].name == "AnatLinear"
