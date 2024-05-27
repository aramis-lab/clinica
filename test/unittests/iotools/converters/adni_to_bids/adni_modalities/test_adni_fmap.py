import json

import pytest

from clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fmap import (
    BIDSFMAPCase,
)


@pytest.mark.parametrize(
    "case, input_value, expected",
    [
        (BIDSFMAPCase.ONE_PHASE_TWO_MAGNITUDES, "e2", "magnitude2"),
        (BIDSFMAPCase.ONE_PHASE_TWO_MAGNITUDES, "e2_ph", "phasediff"),
        (BIDSFMAPCase.TWO_PHASES_TWO_MAGNITUDES, "e1", "magnitude1"),
        (BIDSFMAPCase.TWO_PHASES_TWO_MAGNITUDES, "e2_ph", "phase2"),
    ],
)
def test_phase_magnitude_renamer_success(case, input_value, expected):
    from clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fmap import (
        phase_magnitude_renamer,
    )

    assert phase_magnitude_renamer(input_value, case) == expected


@pytest.mark.parametrize(
    "case, input_value",
    [
        (BIDSFMAPCase.ONE_PHASE_TWO_MAGNITUDES, "_e2"),
        (BIDSFMAPCase.ONE_PHASE_TWO_MAGNITUDES, "e7"),
        (BIDSFMAPCase.ONE_PHASE_TWO_MAGNITUDES, "ee1"),
        (BIDSFMAPCase.TWO_PHASES_TWO_MAGNITUDES, "e2ph"),
        (BIDSFMAPCase.TWO_PHASES_TWO_MAGNITUDES, " e2"),
        (BIDSFMAPCase.TWO_PHASES_TWO_MAGNITUDES, "ea"),
        (BIDSFMAPCase.TWO_PHASES_TWO_MAGNITUDES, "21"),
    ],
)
def test_phase_magnitude_renamer_value_error(case, input_value):
    from clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fmap import (
        phase_magnitude_renamer,
    )

    with pytest.raises(
        ValueError, match=f"Extension {input_value} not taken in charge."
    ):
        phase_magnitude_renamer(input_value, case)


@pytest.mark.parametrize(
    "case, input_value",
    [
        (BIDSFMAPCase.DIRECT_FIELDMAPS, "e2"),
        (BIDSFMAPCase.NOT_SUPPORTED, "e1"),
        (BIDSFMAPCase.EMPTY_FOLDER, "toto"),
    ],
)
def test_phase_magnitude_renamer_implemented_error(case, input_value):
    from clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fmap import (
        phase_magnitude_renamer,
    )

    with pytest.raises(
        NotImplementedError,
        match=f"No renaming should be performed for case {case.value}",
    ):
        phase_magnitude_renamer(input_value, case)


def test_rename_files_success_case1(tmp_path):
    from clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fmap import (
        rename_files,
    )

    fmap = tmp_path / "fmap"
    fmap.mkdir()
    filenames = [
        "sub01-sesM0_fmap_e1.nii.gz",
        "sub01-sesM0_fmap_e1.json",
        "sub01-sesM0_fmap_e2.nii.gz",
        "sub01-sesM0_fmap_e2.json",
        "sub01-sesM0_fmap_e2_ph.nii.gz",
        "sub01-sesM0_fmap_e2_ph.json",
        ".toto.txt",
    ]
    for f in filenames:
        (fmap / f).touch()

    expected_result = {
        fmap / "sub01-sesM0_magnitude1.nii.gz",
        fmap / "sub01-sesM0_magnitude1.json",
        fmap / "sub01-sesM0_magnitude2.nii.gz",
        fmap / "sub01-sesM0_magnitude2.json",
        fmap / "sub01-sesM0_phasediff.nii.gz",
        fmap / "sub01-sesM0_phasediff.json",
        fmap / ".toto.txt",
    }

    rename_files(fmap, BIDSFMAPCase.ONE_PHASE_TWO_MAGNITUDES)

    assert set(fmap.iterdir()) == expected_result


def test_rename_files_success_case2(tmp_path):
    from clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fmap import (
        rename_files,
    )

    fmap = tmp_path / "fmap"
    fmap.mkdir()
    filenames = [
        "sub01_sesM0_fmap_e1.nii.gz",
        "sub01_sesM0_fmap_e1.json",
        "sub01_sesM0_fmap_e2.nii.gz",
        "sub01_sesM0_fmap_e2.json",
        "sub01_sesM0_fmap_e2_ph.nii.gz",
        "sub01_sesM0_fmap_e2_ph.json",
        "sub01_sesM0_fmap_e1_ph.nii.gz",
        "sub01_sesM0_fmap_e1_ph.json",
        ".toto.txt",
    ]
    for f in filenames:
        (fmap / f).touch()

    expected_result = {
        fmap / "sub01_sesM0_magnitude1.nii.gz",
        fmap / "sub01_sesM0_magnitude1.json",
        fmap / "sub01_sesM0_magnitude2.nii.gz",
        fmap / "sub01_sesM0_magnitude2.json",
        fmap / "sub01_sesM0_phase2.nii.gz",
        fmap / "sub01_sesM0_phase2.json",
        fmap / "sub01_sesM0_phase1.nii.gz",
        fmap / "sub01_sesM0_phase1.json",
        fmap / ".toto.txt",
    }

    rename_files(fmap, BIDSFMAPCase.TWO_PHASES_TWO_MAGNITUDES)

    assert set(fmap.iterdir()) == expected_result


@pytest.mark.parametrize(
    "invalid_filename",
    [
        "sub01_sesM0_fmape1.nii.gz",
        "sub01_sesM0_magnitude.nii.gz",
        "sub01_sesM0fmap_e1.nii.gz",
        "toto.txt",
    ],
)
def test_rename_files_error(tmp_path, invalid_filename):
    from clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fmap import (
        rename_files,
    )

    fmap = tmp_path / "fmap"
    fmap.mkdir()
    (fmap / invalid_filename).touch()

    with pytest.raises(ValueError, match=f"Invalid file {invalid_filename} was found."):
        rename_files(fmap, BIDSFMAPCase.TWO_PHASES_TWO_MAGNITUDES)


def test_get_json_file_matching_pattern_success(tmp_path):
    from clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fmap import (
        get_json_file_matching_pattern,
    )

    fmap = tmp_path / "fmap"
    fmap.mkdir()
    for f in [
        "sub01_sesM0_fmap_e2_ph.json",
        "sub01_sesM0_fmap_e2_ph.nii.gz",
        "sub01_sesM0_fmap_e1.json",
        ".toto.txt",
    ]:
        (fmap / f).touch()

    assert (
        get_json_file_matching_pattern(fmap, pattern="ph.json")
        == fmap / "sub01_sesM0_fmap_e2_ph.json"
    )


@pytest.mark.parametrize(
    "filenames, pattern",
    [
        (("sub01_sesM0_fmap_e2.json",), "ph.json"),
        (("sub01_sesM0_fmap_e2_ph.json", "sub01_sesM0_fmap_e1_ph.json"), "ph.json"),
        (("sub01_sesM0_fmap_e2_ph.json"), "fmap.json"),
    ],
)
def test_get_json_file_matching_pattern_error(tmp_path, filenames, pattern):
    from clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fmap import (
        get_json_file_matching_pattern,
    )

    fmap = tmp_path / "fmap"
    fmap.mkdir()
    for f in filenames:
        (fmap / f).touch()

    with pytest.raises(ValueError):
        get_json_file_matching_pattern(fmap, pattern)


@pytest.mark.parametrize(
    "keys",
    [
        ("Echo1",),
        ("Echo1", "Echo2"),
    ],
)
def test_check_json_contains_keys(tmp_path, keys):
    from clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fmap import (
        check_json_contains_keys,
    )

    fmap = tmp_path / "fmap"
    fmap.mkdir()
    json_file_success = fmap / "success.json"
    with open(json_file_success, "w") as f:
        json.dump({key: [] for key in keys}, f)

    json_file_error = fmap / "error.json"
    with open(json_file_error, "w") as f:
        json.dump({"Error": []}, f)

    assert check_json_contains_keys(json_file_success, keys)
    assert not check_json_contains_keys(json_file_error, keys)


@pytest.mark.parametrize(
    "filenames, expected",
    [
        ((".toto.txt",), BIDSFMAPCase.EMPTY_FOLDER),
        (
            (
                "sub01_sesM0_fmap_e2_ph.json",
                "sub01_sesM0_fmap_e2_ph.nii.gz",
                "sub01_sesM0_fmap_e1.json",
                "sub01_sesM0_fmap_e1.nii.gz",
            ),
            BIDSFMAPCase.ONE_PHASE_TWO_MAGNITUDES,
        ),
        (
            (
                "sub01_sesM0_fmap.json",
                "sub01_sesM0_fmap.nii.gz",
                "sub01_sesM0_magnitude.json",
                "sub01_sesM0_magnitude.nii.gz",
            ),
            BIDSFMAPCase.DIRECT_FIELDMAPS,
        ),
        (
            (
                "sub01_sesM0_fmap_e2_ph.json",
                "sub01_sesM0_fmap_e2_ph.nii.gz",
                "sub01_sesM0_fmap_e1.json",
                "sub01_sesM0_fmap_e1.nii.gz",
                "sub01_sesM0_fmap_e2.json",
                "sub01_sesM0_fmap_e2.nii.gz",
            ),
            BIDSFMAPCase.ONE_PHASE_TWO_MAGNITUDES,
        ),
        (
            (
                "sub01_sesM0_fmap_e2_ph.json",
                "sub01_sesM0_fmap_e2_ph.nii.gz",
                "sub01_sesM0_fmap_e1.json",
                "sub01_sesM0_fmap_e1.nii.gz",
                "sub01_sesM0_fmap_e2.json",
                "sub01_sesM0_fmap_e2.nii.gz",
                "sub01_sesM0_fmap_e1_ph.json",
                "sub01_sesM0_fmap_e1_ph.nii.gz",
            ),
            BIDSFMAPCase.TWO_PHASES_TWO_MAGNITUDES,
        ),
        (
            (
                "sub01_sesM0_fmap_e2_ph.json",
                "sub01_sesM0_fmap_e2_ph.nii.gz",
            ),
            BIDSFMAPCase.NOT_SUPPORTED,
        ),
        (
            (
                "sub01_sesM0_fmap_e2_ph.json",
                "sub01_sesM0_fmap_e2_ph.nii.gz",
                "sub01_sesM0_fmap_e1.json",
            ),
            BIDSFMAPCase.NOT_SUPPORTED,
        ),
    ],
)
def test_infer_case_fmap(tmp_path, filenames, expected):
    from clinica.iotools.converters.adni_to_bids.adni_modalities.adni_fmap import (
        infer_case_fmap,
    )

    fmap = tmp_path / "fmap"
    fmap.mkdir()
    for f in filenames:
        (fmap / f).touch()

    assert infer_case_fmap(fmap) == expected
