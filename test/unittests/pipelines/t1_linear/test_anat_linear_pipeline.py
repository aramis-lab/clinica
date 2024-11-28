import re

import pytest
from packaging.version import Version

from clinica.utils.bids import Visit
from clinica.utils.testing_utils import build_bids_directory, build_caps_directory


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

    bids = build_bids_directory(tmp_path / "bids", {"sub-01": ["ses-M00"]})
    caps = build_caps_directory(tmp_path / "caps", {})
    AnatLinear(bids_directory=str(bids), caps_directory=str(caps))

    files = [f for f in caps.iterdir()]
    assert len(files) == 1

    desc = CAPSDatasetDescription.from_file(caps / "dataset_description.json")

    assert desc.bids_version == Version("1.7.0")
    assert desc.caps_version == Version("1.0.0")
    assert desc.dataset_type == DatasetType.DERIVATIVE
    assert desc.processing[0].name == "AnatLinear"


def test_anat_linear_info_loading(tmp_path):
    from clinica.pipelines.t1_linear.anat_linear_pipeline import AnatLinear

    bids = build_bids_directory(tmp_path / "bids", {"sub-01": ["ses-M00"]})
    caps = build_caps_directory(tmp_path / "caps", {})
    pipeline = AnatLinear(bids_directory=str(bids), caps_directory=str(caps))

    assert pipeline.info == {
        "author": "Mauricio Diaz",
        "dependencies": [
            {
                "name": "ants",
                "type": "software",
                "version": ">=2.2.0",
            },
        ],
        "id": "aramislab/t1-linear",
        "space_caps": "45M",
        "space_wd": "45M",
        "version": "0.1.0",
    }


def test_anat_linear_dependencies(tmp_path, mocker):
    from packaging.specifiers import SpecifierSet

    from clinica.pipelines.t1_linear.anat_linear_pipeline import AnatLinear
    from clinica.utils.check_dependency import SoftwareDependency, ThirdPartySoftware

    mocker.patch(
        "clinica.utils.check_dependency._get_ants_version",
        return_value=Version("2.2.1"),
    )
    bids = build_bids_directory(tmp_path / "bids", {"sub-01": ["ses-M00"]})
    caps = build_caps_directory(tmp_path / "caps", {})
    pipeline = AnatLinear(bids_directory=str(bids), caps_directory=str(caps))

    assert pipeline.dependencies == [
        SoftwareDependency(
            ThirdPartySoftware.ANTS, SpecifierSet(">=2.2.0"), Version("2.2.1")
        )
    ]

    # When using AntsPy, the ANTs dependency is not considered
    pipeline = AnatLinear(
        bids_directory=str(bids), caps_directory=str(caps), use_antspy=True
    )

    assert pipeline.dependencies == []


@pytest.mark.parametrize(
    "config,expected",
    [
        ({}, []),
        (
            {
                "pipelines": {"t1_linear": {"uncropped_image": False}},
                "subjects": {"sub-01": ["ses-M006"]},
            },
            [Visit("sub-01", "ses-M006")],
        ),
        (
            {
                "pipelines": {"t1_linear": {"uncropped_image": False}},
                "subjects": {"sub-01": ["ses-M000", "ses-M006"]},
            },
            [Visit("sub-01", "ses-M000"), Visit("sub-01", "ses-M006")],
        ),
        (
            {
                "pipelines": {"t1_linear": {"uncropped_image": False}},
                "subjects": {
                    "sub-01": ["ses-M000", "ses-M006"],
                    "sub-02": ["ses-M000"],
                },
            },
            [
                Visit("sub-01", "ses-M000"),
                Visit("sub-01", "ses-M006"),
                Visit("sub-02", "ses-M000"),
            ],
        ),
    ],
)
def test_anat_linear_get_processed_visits_uncropped_images(
    tmp_path, mocker, config, expected
):
    from clinica.pipelines.t1_linear.anat_linear_pipeline import AnatLinear

    mocker.patch(
        "clinica.utils.check_dependency._get_ants_version",
        return_value=Version("2.2.1"),
    )
    bids = build_bids_directory(
        tmp_path / "bids", {"sub-01": ["ses-M000", "ses-M006"], "sub-02": ["ses-M000"]}
    )
    caps = build_caps_directory(tmp_path / "caps", config)
    pipeline = AnatLinear(
        bids_directory=str(bids),
        caps_directory=str(caps),
        parameters={"uncropped_image": False},
    )
    pipeline2 = AnatLinear(
        bids_directory=str(bids),
        caps_directory=str(caps),
        parameters={"uncropped_image": True},
    )

    assert pipeline.get_processed_images() == expected
    # pipeline2 always find an empty list of processed images because we want un-cropped images
    # and the CAPS folder only contains cropped images
    assert pipeline2.get_processed_images() == []


@pytest.mark.parametrize(
    "config,expected",
    [
        ({}, []),
        (
            {
                "pipelines": {"t1_linear": {"uncropped_image": True}},
                "subjects": {"sub-01": ["ses-M006"]},
            },
            [Visit("sub-01", "ses-M006")],
        ),
        (
            {
                "pipelines": {"t1_linear": {"uncropped_image": True}},
                "subjects": {"sub-01": ["ses-M000", "ses-M006"]},
            },
            [Visit("sub-01", "ses-M000"), Visit("sub-01", "ses-M006")],
        ),
        (
            {
                "pipelines": {"t1_linear": {"uncropped_image": True}},
                "subjects": {
                    "sub-01": ["ses-M000", "ses-M006"],
                    "sub-02": ["ses-M000"],
                },
            },
            [
                Visit("sub-01", "ses-M000"),
                Visit("sub-01", "ses-M006"),
                Visit("sub-02", "ses-M000"),
            ],
        ),
    ],
)
def test_anat_linear_get_processed_visits_cropped_images(
    tmp_path, mocker, config, expected
):
    from clinica.pipelines.t1_linear.anat_linear_pipeline import AnatLinear

    mocker.patch(
        "clinica.utils.check_dependency._get_ants_version",
        return_value=Version("2.2.1"),
    )
    bids = build_bids_directory(
        tmp_path / "bids", {"sub-01": ["ses-M000", "ses-M006"], "sub-02": ["ses-M000"]}
    )
    caps = build_caps_directory(tmp_path / "caps", config)
    pipeline = AnatLinear(
        bids_directory=str(bids),
        caps_directory=str(caps),
        parameters={"uncropped_image": False},
    )
    pipeline2 = AnatLinear(
        bids_directory=str(bids),
        caps_directory=str(caps),
        parameters={"uncropped_image": True},
    )

    assert pipeline2.get_processed_images() == expected
    # pipeline always find an empty list of processed images because we want cropped images
    # and the CAPS folder only contains un-cropped images
    assert pipeline.get_processed_images() == []
