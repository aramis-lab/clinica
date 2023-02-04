import pytest


@pytest.mark.parametrize(
    "input_list,expected",
    [
        (
            ["ses-M000", "ses-M006", "ses-M012", "ses-M024", "ses-M048", "ses-M003"],
            ["ses-M000", "ses-M003", "ses-M006", "ses-M012", "ses-M024", "ses-M048"],
        ),
        (
            ["ses-M00", "ses-M06", "ses-M12", "ses-M24", "ses-M48", "ses-M03"],
            ["ses-M00", "ses-M03", "ses-M06", "ses-M12", "ses-M24", "ses-M48"],
        ),
        (
            ["ses-M0", "ses-M6", "ses-M12", "ses-M24", "ses-M48", "ses-M3"],
            ["ses-M0", "ses-M3", "ses-M6", "ses-M12", "ses-M24", "ses-M48"],
        ),
    ],
)
def test_sort_session_list(input_list, expected):
    """Test function `sort_session_list`."""
    from clinica.iotools.converter_utils import sort_session_list

    assert sort_session_list(input_list) == expected


@pytest.mark.parametrize(
    "input,expected",
    [
        ("bl", "ses-M000"),
        ("m0", "ses-M000"),
        ("m3", "ses-M003"),
        ("m03", "ses-M003"),
        ("m003", "ses-M003"),
        ("m0003", "ses-M003"),
        ("m00", "ses-M000"),
        ("m0000000", "ses-M000"),
    ],
)
def test_viscode_to_session(input, expected):
    """Test function `viscode_to_session`."""

    from clinica.iotools.converter_utils import viscode_to_session

    assert viscode_to_session(input) == expected


@pytest.mark.parametrize(
    "sessions,modalities,expected",
    [
        ([], None, {}),
        (
            ["ses-M000", "ses-M006"],
            None,
            {
                k: {
                    "session": 0,
                    "dwi": 0,
                    "func": 0,
                    "fieldmap": 0,
                    "flair": 0,
                    "t1w": 0,
                }
                for k in ["ses-M000", "ses-M006"]
            },
        ),
        ([], ["fmri"], {}),
        (
            ["ses-M000", "ses-M006"],
            ["fmri", "dwi"],
            {k: {"session": 0, "dwi": 0, "fmri": 0} for k in ["ses-M000", "ses-M006"]},
        ),
    ],
)
def test_missing_mods_tracker_instantiation(sessions, modalities, expected):
    from clinica.iotools.converter_utils import MissingModsTracker

    tracker = MissingModsTracker(sessions, modalities)
    assert tracker.missing == expected
    assert tracker.ses == sessions
    assert tracker.get_missing_list() == expected


def test_missing_mods_tracker_add_errors():
    from clinica.iotools.converter_utils import MissingModsTracker

    tracker = MissingModsTracker(["ses-M000", "ses-M006"], ["fmri", "dwi"])

    with pytest.raises(
        ValueError,
        match="Session foo was not provided to the MissingModsTracker constructor.",
    ):
        tracker.add_missing_mod("foo", "fmri")
        tracker.increase_missing_ses("foo")

    with pytest.raises(
        ValueError,
        match="Modality foo is not tracked by this instance of MissingModsTracker.",
    ):
        tracker.add_missing_mod("ses-M006", "foo")


def test_missing_mods_tracker_add():
    from clinica.iotools.converter_utils import MissingModsTracker

    tracker = MissingModsTracker(["ses-M000", "ses-M006"], ["fmri", "dwi"])
    tracker.increase_missing_ses("ses-M000")
    assert tracker.missing["ses-M000"]["session"] == 1
    assert tracker.missing["ses-M006"]["session"] == 0
    tracker.add_missing_mod("ses-M006", "dwi")
    assert tracker.missing["ses-M006"]["dwi"] == 1
    assert tracker.missing["ses-M000"]["dwi"] == 0
    tracker.add_missing_mod("ses-M000", "dwi")
    assert tracker.missing["ses-M006"]["dwi"] == 1
    assert tracker.missing["ses-M000"]["dwi"] == 1


def test_compute_statistics():
    from clinica.iotools.converter_utils import MissingModsTracker, compute_statistics

    tracker = MissingModsTracker(["ses-M000", "ses-M006"], ["fmri", "dwi"])
    tracker.increase_missing_ses("ses-M000")
    tracker.add_missing_mod("ses-M006", "dwi")
    assert compute_statistics(2, ["ses-M000", "ses-M006"], tracker) == (
        "**********************************************\n"
        "Number of subjects converted: 2\n"
        "Sessions available: ['ses-M000', 'ses-M006']\n\n"
        "Number of sessions ses-M000 found: 1 (50.0%)\n\n"
        "Number of sessions ses-M006 found: 2 (100.0%)\n\n"
        "**********************************************\n\n"
        "Number of missing modalities for each session:\n\n"
        "ses-M000\nfmri: 0 (0.0%) \ndwi: 0 (0.0%) \n\n"
        "ses-M006\nfmri: 0 (0.0%) \ndwi: 1 (50.0%) \n"
    )


def test_compute_table():
    from clinica.iotools.converter_utils import compute_table

    mods = {
        "foo": {
            "missing": 42,
            "n/a": 16,
            "bar": 23,
        },
        "bar": {
            "missing": 0,
            "baz": 7,
        },
    }
    assert compute_table(mods) == (
        "\tbar\t| baz\t| missing\t| n/a\n------------------------------------------------\n"
        "foo\t23\t| 0\t| 42\t| 16\nbar\t0\t| 7\t| 0\t| 0\n"
    )


def test_compute_longitudinal_analysis(tmp_path):
    import pandas as pd

    from clinica.iotools.converter_utils import compute_longitudinal_analysis

    sessions = ["ses-M000", "ses-M006"]
    for subject in ("sub-01", "sub-03"):
        (tmp_path / "bids" / subject).mkdir(parents=True)
    for ses in sessions:
        df = pd.DataFrame(
            [["sub-01", 1, 0, 1], ["sub-03", 0, 1, 1]],
            columns=["participant_id", "foo", "bar", "baz"],
        )
        df.to_csv(tmp_path / f"missing_mods_{ses}.tsv", sep="\t", index=False)
    df = pd.DataFrame(
        [
            ["ses-M000", 18, "foooo", "CN", "ADNI"],
            ["ses-M006", 19, "foooo", "AD", "ADNI"],
        ],
        columns=["session_id", "age", "foobarbaz", "diagnosis", "study"],
    )
    df.to_csv(
        tmp_path / "bids" / "sub-01" / "sub-01_sessions.tsv", sep="\t", index=False
    )
    df = pd.DataFrame(
        [["ses-M000", 66, "bar", 42, "ADNI"], ["ses-M012", 67, "baz", 69, "ADNI"]],
        columns=["session_id", "age", "bar", "diagnosis", "study"],
    )
    df.to_csv(
        tmp_path / "bids" / "sub-03" / "sub-03_sessions.tsv", sep="\t", index=False
    )
    assert compute_longitudinal_analysis(
        tmp_path / "bids", tmp_path, sessions, "missing_mods_"
    ) == (
        "**********************************************\n\n"
        "Number of present diagnoses and modalities "
        "for each session:\nses-M000\n\tCN\t| n/a\n"
        "--------------------------------\nfoo\t1\t| 0\nbar\t0\t| "
        "1\nbaz\t1\t| 1\n\n\nses-M006\n\tAD\t| missing\n"
        "--------------------------------\nfoo\t1\t| "
        "0\nbar\t0\t| 1\nbaz\t1\t| 1\n\n\n"
    )
