# coding: utf8

"""
This file contains a set of functional tests designed to check the correct execution of the pipeline and the
different functions available in Clinica
"""


import warnings
from os import PathLike
from pathlib import Path
from test.nonregression.testing_tools import compare_folders

import pytest

# Determine location for working_directory
warnings.filterwarnings("ignore")


@pytest.fixture(
    params=[
        # TODO: Update NIFD reference dataset.
        "Nifd2Bids",
        "Oasis2Bids",
        "Oasis3ToBids",
        "Adni2Bids",
        "Aibl2Bids",
        "HabsToBids",
        "UkbToBids",
        "GenfiToBids",
    ]
)
def test_name(request):
    return request.param


def run_nifd2bids(input_dir: PathLike, output_dir: PathLike, ref_dir: PathLike) -> None:
    from pathlib import PurePath
    from tempfile import TemporaryDirectory

    from clinica.iotools.converters.nifd_to_bids.nifd_to_bids import convert_images

    # Convert
    input_dir = PurePath(input_dir)
    output_dir = PurePath(output_dir)
    ref_dir = PurePath(ref_dir)

    # Act
    _ = convert_images(
        path_to_clinical=input_dir / "clinical_data",
        path_to_dataset=input_dir / "unorganized",
        bids_dir=output_dir,
    )

    # Assert
    with TemporaryDirectory() as td:
        compare_folders(output_dir, ref_dir, td)


def run_oasis2bids(
    input_dir: PathLike, output_dir: PathLike, ref_dir: PathLike
) -> None:
    from pathlib import PurePath

    from clinica.iotools.converters.oasis_to_bids.oasis_to_bids import OasisToBids

    # Convert
    input_dir = PurePath(input_dir)
    output_dir = PurePath(output_dir)
    ref_dir = PurePath(ref_dir)

    # Arrange
    clinical_data_directory = input_dir / "clinical_data"
    # Act
    oasis_to_bids = OasisToBids()
    oasis_to_bids.convert_images(input_dir / "unorganized", output_dir / "bids")
    oasis_to_bids.convert_clinical_data(clinical_data_directory, output_dir / "bids")
    # Assert
    compare_folders(output_dir / "bids", ref_dir / "bids", output_dir)


def run_oasis3tobids(
    input_dir: PathLike, output_dir: PathLike, ref_dir: PathLike
) -> None:
    from pathlib import PurePath

    from clinica.iotools.converters.oasis3_to_bids.oasis3_to_bids import convert_images

    # Convert
    input_dir = PurePath(input_dir)
    output_dir = PurePath(output_dir)
    ref_dir = PurePath(ref_dir)

    # Arrange
    clinical_data_directory = input_dir / "clinical_data"

    # Act
    convert_images(
        input_dir / "unorganized", output_dir / "bids", clinical_data_directory
    )

    # Assert
    compare_folders(output_dir / "bids", ref_dir / "bids", output_dir)


def run_adni2bids(input_dir: PathLike, output_dir: PathLike, ref_dir: PathLike) -> None:
    from pathlib import PurePath

    from clinica.iotools.converters.adni_to_bids.adni_to_bids import AdniToBids

    # Convert
    input_dir = PurePath(input_dir)
    output_dir = PurePath(output_dir)
    ref_dir = PurePath(ref_dir)

    # Arrange
    clinical_data_directory = input_dir / "clinical_data"
    xml_directory = input_dir / "xml_metadata"
    dataset_directory = input_dir / "unorganized_data"
    subjects_list = input_dir / "subjects.txt"
    modalities = ["T1", "PET_FDG", "PET_AMYLOID", "PET_TAU", "DWI", "FLAIR", "fMRI"]

    # Act
    adni_to_bids = AdniToBids()
    adni_to_bids.check_adni_dependencies()
    adni_to_bids.convert_images(
        dataset_directory,
        clinical_data_directory,
        output_dir / "bids",
        subjects_list,
        modalities,
    )
    adni_to_bids.convert_clinical_data(
        clinical_data_directory,
        output_dir / "bids",
        xml_path=xml_directory,
    )

    # Assert
    compare_folders(output_dir / "bids", ref_dir / "bids", output_dir)


def run_aibl2bids(input_dir: PathLike, output_dir: PathLike, ref_dir: PathLike) -> None:
    from pathlib import PurePath

    from clinica.iotools.converters.aibl_to_bids.aibl_to_bids import (
        convert_clinical_data,
        convert_images,
    )

    # Convert
    input_dir = PurePath(input_dir)
    output_dir = PurePath(output_dir)
    ref_dir = PurePath(ref_dir)

    # Arrange
    clinical_data_directory = input_dir / "Data_extract_3.2.5"
    dataset_directory = input_dir / "unorganized_data"

    # Act
    convert_images(
        dataset_directory,
        clinical_data_directory,
        output_dir / "bids",
    )
    convert_clinical_data(output_dir / "bids", clinical_data_directory)

    # Assert
    compare_folders(output_dir / "bids", ref_dir / "bids", output_dir)


def run_habs_to_bids(
    input_dir: PathLike, output_dir: PathLike, ref_dir: PathLike
) -> None:
    from click.testing import CliRunner

    from clinica.iotools.converters.habs_to_bids.habs_to_bids_cli import cli

    runner = CliRunner()
    result = runner.invoke(cli, [str(input_dir), str(output_dir)])

    assert result.exit_code == 0
    compare_folders(output_dir, ref_dir, output_dir)


def run_ukbtobids(input_dir: PathLike, output_dir: PathLike, ref_dir: PathLike) -> None:
    from pathlib import PurePath

    from clinica.iotools.converters.ukb_to_bids.ukb_to_bids import convert_images
    from clinica.utils.check_dependency import check_dcm2niix

    # Convert
    input_dir = PurePath(input_dir)
    output_dir = PurePath(output_dir)
    ref_dir = PurePath(ref_dir)

    # Arrange
    clinical_data_directory = input_dir / "clinical_data"

    # Act
    check_dcm2niix()
    convert_images(
        input_dir / "unorganized", output_dir / "bids", clinical_data_directory
    )

    # Assert
    compare_folders(output_dir / "bids", ref_dir / "bids", output_dir / "bids")


def run_genfitobids(
    input_dir: PathLike, output_dir: PathLike, ref_dir: PathLike
) -> None:
    from pathlib import PurePath

    from clinica.iotools.converters.genfi_to_bids.genfi_to_bids import convert_images
    from clinica.utils.check_dependency import check_dcm2niix

    # Convert
    input_dir = PurePath(input_dir)
    output_dir = PurePath(output_dir)
    ref_dir = PurePath(ref_dir)

    # Arrange
    clinical_data_directory = input_dir / "clinical_data"

    # Act
    check_dcm2niix()
    convert_images(
        input_dir / "unorganized", output_dir / "bids", clinical_data_directory
    )

    # Assert
    compare_folders(output_dir / "bids", ref_dir / "bids", output_dir / "bids")


def test_run_convertors(cmdopt, tmp_path, test_name):
    base_dir = Path(cmdopt["input"])
    input_dir = base_dir / test_name / "in"
    ref_dir = base_dir / test_name / "ref"
    tmp_out_dir = tmp_path / test_name / "out"
    tmp_out_dir.mkdir(parents=True)

    if test_name == "Nifd2Bids":
        run_nifd2bids(input_dir, tmp_out_dir, ref_dir)

    elif test_name == "Oasis2Bids":
        run_oasis2bids(input_dir, tmp_out_dir, ref_dir)

    elif test_name == "Oasis3ToBids":
        run_oasis3tobids(input_dir, tmp_out_dir, ref_dir)

    elif test_name == "Adni2Bids":
        run_adni2bids(input_dir, tmp_out_dir, ref_dir)

    elif test_name == "Aibl2Bids":
        run_aibl2bids(input_dir, tmp_out_dir, ref_dir)

    elif test_name == "HabsToBids":
        run_habs_to_bids(input_dir, tmp_out_dir, ref_dir)
    elif test_name == "UkbToBids":
        run_ukbtobids(input_dir, tmp_out_dir, ref_dir)
    elif test_name == "GenfiToBids":
        run_genfitobids(input_dir, tmp_out_dir, ref_dir)

    else:
        print(f"Test {test_name} not available.")
        assert 0
