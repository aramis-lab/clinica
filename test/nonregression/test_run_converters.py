"""
This file contains a set of functional tests designed to check the
correct execution of the dataset converters available in Clinica.
"""

from pathlib import Path
from test.nonregression.testing_tools import (
    compare_bids_tsv,
    compare_folders,
    configure_paths,
)

import pytest

from clinica.converters.study_models import StudyName


@pytest.mark.parametrize("study", StudyName)
def test_converters(cmdopt, tmp_path, study: StudyName):
    from clinica.converters.factory import convert, get_converter_name

    base_dir = Path(cmdopt["input"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, get_converter_name(study)
    )
    output_dir = tmp_path / "bids"

    convert(
        study,
        path_to_dataset=input_dir / "unorganized",
        bids_dir=output_dir,
        path_to_clinical=input_dir / "clinical_data",
        subjects=input_dir / "subjects.txt",
        xml_path=input_dir / "xml_metadata",
    )
    compare_folders(output_dir, ref_dir / "bids", output_dir)
    if study in (
        StudyName.AIBL,
        StudyName.IXI,
        StudyName.OASIS,
        StudyName.HABS,
        StudyName.NIFD,
        StudyName.OASIS3,
        StudyName.GENFI,
        StudyName.UKB,
    ):
        compare_bids_tsv(output_dir, ref_dir / "bids")


@pytest.mark.parametrize(
    "gif,full",
    [
        (False, False),
        (True, False),
        (True, True),
    ],
)
def test_genfi_converter_using_options(cmdopt, tmp_path, gif, full):
    from clinica.converters.factory import convert, get_converter_name

    study = StudyName.GENFI
    base_dir = Path(cmdopt["input"])
    input_dir, tmp_dir, ref_dir = configure_paths(
        base_dir, tmp_path, get_converter_name(study)
    )
    output_dir = tmp_path / "bids"
    convert(
        study,
        path_to_dataset=input_dir / "unorganized",
        bids_dir=output_dir,
        path_to_clinical=input_dir / "clinical_data",
        subjects=input_dir / "subjects.txt",
        path_to_clinical_txt=input_dir / "additional_clinical_data.txt",
        gif=gif,
        full=full,
    )

    options = "full" if full else "gif" if gif else "cdt"
    compare_folders(
        output_dir, ref_dir / "bids_with_options" / f"bids_{options}", output_dir
    )
    compare_bids_tsv(output_dir, ref_dir / "bids_with_options" / f"bids_{options}")
