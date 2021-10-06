# coding: utf8

"""
This file contains a set of functional tests designed to check the correct execution of the pipeline and the
different functions available in Clinica
"""


import warnings
from os import PathLike, fspath
from pathlib import Path
from test.nonregression.testing_tools import (
    clean_folder,
    compare_folders,
    compare_folders_structures,
    compare_folders_with_hashes,
    create_list_hashes,
    identical_subject_list,
    same_missing_modality_tsv,
)

import pytest

# Determine location for working_directory
warnings.filterwarnings("ignore")


@pytest.fixture(
    params=[
        "Nifd2Bids",
        "Oasis2Bids",
        "Oasis3ToBids",
        "Adni2Bids",
        "Aibl2Bids",
    ]
)
def test_name(request):
    return request.param


def run_nifd2bids(input_dir: PathLike, output_dir: PathLike, ref_dir: PathLike) -> None:
    import shutil

    from clinica.iotools.converters.nifd_to_bids.nifd_to_bids import convert_images

    # Arrange
    shutil.copytree(
        input_dir / "clinical_data",
        output_dir / "clinical_data",
        copy_function=shutil.copy,
    )
    # Arrange - Data location
    clinical_data_directory = output_dir / "clinical_data"
    # Acte - Conversion
    to_convert = convert_images(
        input_dir / "unorganized", output_dir / "bids", clinical_data_directory
    )
    # Assert
    compare_folders_structures(
        fspath(output_dir / "bids"), fspath(ref_dir / "hashes_nifd.p")
    )


def run_oasis2bids(
    input_dir: PathLike, output_dir: PathLike, ref_dir: PathLike
) -> None:
    from clinica.iotools.converters.oasis_to_bids.oasis_to_bids import OasisToBids

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
    from clinica.iotools.converters.oasis3_to_bids.oasis3_to_bids import convert_images

    # Arrange
    clinical_data_directory = input_dir / "clinical_data"
    # Act
    convert_images(
        input_dir / "unorganized", output_dir / "bids", clinical_data_directory
    )
    # Assert
    compare_folders(output_dir / "bids", ref_dir / "bids", output_dir)


def run_adni2bids(input_dir: PathLike, output_dir: PathLike, ref_dir: PathLike) -> None:
    from clinica.iotools.converters.adni_to_bids.adni_to_bids import AdniToBids

    # Arrange
    clinical_data_directory = input_dir / "clinical_data"
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
    adni_to_bids.convert_clinical_data(clinical_data_directory, output_dir / "bids")
    # Assert
    compare_folders(output_dir / "bids", ref_dir / "bids", output_dir)


def run_aibl2bids(input_dir: PathLike, output_dir: PathLike, ref_dir: PathLike) -> None:
    from clinica.iotools.converters.aibl_to_bids.aibl_to_bids import (
        convert_clinical_data,
        convert_images,
    )

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


def test_run_convertors(cmdopt, tmp_path, test_name):
    import shutil

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

    else:
        print(f"Test {test_name} not available.")
        assert 0
