"""
This file contains a set of functional tests designed to check the
correct execution of the dataset converters available in Clinica.
"""

from pathlib import Path
from test.nonregression.testing_tools import compare_folders, configure_paths

import pytest

from clinica.iotools.bids_utils import StudyName


@pytest.mark.parametrize("study", StudyName)
def test_converters(cmdopt, tmp_path, study: StudyName):
    from clinica.iotools.converters.factory import convert, get_converter_name

    # To be removed once we have built and deployed our testing dataset
    if study == StudyName.IXI:
        pytest.skip("Non regression tests are not yet implemented for ixi-to-bids.")

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
