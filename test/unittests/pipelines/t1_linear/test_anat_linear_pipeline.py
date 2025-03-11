import re

import pytest
from packaging.version import Version

from clinica.dataset import Visit
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
    from clinica.utils.exceptions import ClinicaDatasetError

    with pytest.raises(
        ClinicaDatasetError,
        match=re.escape(
            f"The directory ({tmp_path}) you provided is missing a dataset_description.json file."
        ),
    ):
        AnatLinear(bids_directory=str(tmp_path))


def test_anat_linear_pipeline_single_caps_input_error(tmp_path):
    from clinica.pipelines.t1_linear.anat_linear_pipeline import AnatLinear
    from clinica.utils.exceptions import ClinicaDatasetError

    with pytest.raises(
        ClinicaDatasetError,
        match=re.escape(
            f"The directory ({tmp_path}) you provided is missing a dataset_description.json file."
        ),
    ):
        AnatLinear(caps_directory=str(tmp_path))


def test_anat_linear_pipeline_write_caps_dataset_description(tmp_path):
    from clinica.dataset import CAPSDatasetDescription, DatasetType
    from clinica.pipelines.t1_linear.anat_linear_pipeline import AnatLinear

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

    assert pipeline.get_processed_visits() == expected
    # pipeline2 always find an empty list of processed images because we want un-cropped images
    # and the CAPS folder only contains cropped images
    assert pipeline2.get_processed_visits() == []


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

    assert pipeline2.get_processed_visits() == expected
    # pipeline always find an empty list of processed images because we want cropped images
    # and the CAPS folder only contains un-cropped images
    assert pipeline.get_processed_visits() == []


@pytest.mark.parametrize(
    "config,remaining_subjects,remaining_sessions",
    [
        ({}, ["sub-01", "sub-01", "sub-02"], ["ses-M000", "ses-M006", "ses-M000"]),
        (
            {
                "pipelines": {"t1_linear": {"uncropped_image": False}},
                "subjects": {"sub-01": ["ses-M006"]},
            },
            ["sub-01", "sub-02"],
            ["ses-M000", "ses-M000"],
        ),
        (
            {
                "pipelines": {"t1_linear": {"uncropped_image": False}},
                "subjects": {"sub-01": ["ses-M000", "ses-M006"]},
            },
            ["sub-02"],
            ["ses-M000"],
        ),
        (
            {
                "pipelines": {"t1_linear": {"uncropped_image": False}},
                "subjects": {
                    "sub-01": ["ses-M000", "ses-M006"],
                    "sub-02": ["ses-M000"],
                },
            },
            [],
            [],
        ),
    ],
)
def test_determine_subject_and_session_to_process(
    tmp_path, mocker, config, remaining_subjects, remaining_sessions
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
    pipeline.determine_subject_and_session_to_process()
    pipeline2 = AnatLinear(
        bids_directory=str(bids),
        caps_directory=str(caps),
        parameters={"uncropped_image": True},
    )
    pipeline2.determine_subject_and_session_to_process()

    assert pipeline.subjects == remaining_subjects
    assert pipeline.sessions == remaining_sessions
    assert pipeline2.subjects == ["sub-01", "sub-01", "sub-02"]
    assert pipeline2.sessions == ["ses-M000", "ses-M006", "ses-M000"]


@pytest.mark.parametrize(
    "configuration,expected_message",
    [
        (
            {
                "pipelines": {"t1_linear": {"uncropped_image": False}},
                "subjects": {
                    "sub-01": ["ses-M000", "ses-M006"],
                    "sub-02": ["ses-M000"],
                },
            },
            (
                "Clinica found already processed images for 3 visit(s):"
                "\n- sub-01 ses-M000\n- sub-01 ses-M006\n- sub-02 ses-M000"
                "\nThose visits will be ignored by Clinica."
            ),
        ),
        (
            {
                "pipelines": {"t1_linear": {"uncropped_image": False}},
                "subjects": {"sub-01": ["ses-M000", "ses-M006"]},
            },
            (
                "Clinica found already processed images for 2 visit(s):"
                "\n- sub-01 ses-M000\n- sub-01 ses-M006"
                "\nThose visits will be ignored by Clinica."
            ),
        ),
    ],
)
def test_determine_subject_and_session_to_process_warning(
    tmp_path, mocker, configuration, expected_message
):
    from clinica.pipelines.t1_linear.anat_linear_pipeline import AnatLinear

    mocker.patch(
        "clinica.utils.check_dependency._get_ants_version",
        return_value=Version("2.2.1"),
    )
    bids = build_bids_directory(
        tmp_path / "bids", {"sub-01": ["ses-M000", "ses-M006"], "sub-02": ["ses-M000"]}
    )
    caps = build_caps_directory(tmp_path / "caps", configuration)
    pipeline = AnatLinear(
        bids_directory=str(bids),
        caps_directory=str(caps),
        parameters={"uncropped_image": False},
    )
    with pytest.warns(
        UserWarning,
        match=re.escape(f"In the provided CAPS folder {caps}, {expected_message}"),
    ):
        pipeline.determine_subject_and_session_to_process()
