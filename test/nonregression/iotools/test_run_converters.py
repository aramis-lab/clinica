# coding: utf8

"""
This file contains a set of functional tests designed to check the correct execution of the pipeline and the
different functions available in Clinica
"""


import warnings
from os import fspath
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
        # "Oasis3ToBids",
        "Adni2Bids",
        "Aibl2Bids",
    ]
)
def name(request):
    return request.param


def test_run_convertors(cmdopt, tmp_path, name):
    import shutil

    base_dir = Path(cmdopt["input"])
    input_dir = base_dir / name / "in"
    ref_dir = base_dir / name / "ref"
    tmp_out_dir = tmp_path / name / "out"
    tmp_out_dir.mkdir(parents=True)
    dataset_directory = input_dir / "unorganized"
    bids_directory = tmp_out_dir / "bids"

    if name == "Nifd2Bids":
        from clinica.iotools.converters.nifd_to_bids.nifd_to_bids import convert_images

        # Arrange
        shutil.copytree(
            input_dir / "clinical_data",
            tmp_out_dir / "clinical_data",
            copy_function=shutil.copy,
        )
        # Arrange - Data location
        clinical_data_directory = tmp_out_dir / "clinical_data"
        # Acte - Conversion
        to_convert = convert_images(
            dataset_directory, bids_directory, clinical_data_directory
        )
        # Assert
        compare_folders_structures(
            fspath(bids_directory), fspath(ref_dir / "hashes_nifd.p")
        )
    elif name == "Oasis2Bids":
        from clinica.iotools.converters.oasis_to_bids.oasis_to_bids import OasisToBids

        # Arrange
        clinical_data_directory = input_dir / "clinical_data"
        # Act
        oasis_to_bids = OasisToBids()
        oasis_to_bids.convert_images(dataset_directory, bids_directory)
        oasis_to_bids.convert_clinical_data(clinical_data_directory, bids_directory)
        # Assert
        compare_folders(tmp_out_dir / "bids", ref_dir / "bids", tmp_path)
    elif name == "Oasis3ToBids":
        from clinica.iotools.converters.oasis3_to_bids.oasis3_to_bids import (
            Oasis3ToBids,
        )

        # Arrange
        clinical_data_directory = input_dir / "clinical_data"
        # Act
        oasis_to_bids = Oasis3ToBids()
        oasis_to_bids.convert_images(dataset_directory, bids_directory)
        oasis_to_bids.convert_clinical_data(clinical_data_directory, bids_directory)
        # Assert
        compare_folders(tmp_out_dir / "bids", ref_dir / "bids", tmp_path)
    elif name == "Adni2Bids":
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
            bids_directory,
            subjects_list,
            modalities,
        )
        adni_to_bids.convert_clinical_data(clinical_data_directory, bids_directory)
        # Assert
        compare_folders(tmp_out_dir / "bids", ref_dir / "bids", tmp_path)
    elif name == "Aibl2Bids":
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
            bids_directory,
        )
        convert_clinical_data(bids_directory, clinical_data_directory)
        # Assert
        compare_folders(fspath(tmp_out_dir), fspath(ref_dir), shared_folder_name="bids")
    else:
        print(f"Test {name} not available.")
        assert 0
