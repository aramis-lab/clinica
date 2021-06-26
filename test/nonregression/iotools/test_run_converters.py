# coding: utf8

"""
This file contains a set of functional tests designed to check the correct execution of the pipeline and the
different functions available in Clinica
"""


import warnings
from os import pardir
from test.nonregression.testing_tools import (
    clean_folder,
    compare_folders,
    compare_folders_structures,
    compare_folders_with_hashes,
    create_list_hashes,
    identical_subject_list,
    same_missing_modality_tsv,
)

# Determine location for working_directory
warnings.filterwarnings("ignore")


def test_run_Oasis2Bids(cmdopt):
    from os.path import abspath, dirname, join

    from clinica.iotools.converters.oasis_to_bids.oasis_to_bids import OasisToBids

    root = join(dirname(abspath(__file__)), pardir, pardir, "data", "Oasis2Bids")

    clean_folder(join(root, "out", "bids"), recreate=True)

    # Data location
    dataset_directory = join(root, "in", "unorganized")
    bids_directory = join(root, "out", "bids")
    clinical_data_directory = join(root, "in", "clinical_data")

    oasis_to_bids = OasisToBids()
    oasis_to_bids.convert_images(dataset_directory, bids_directory)
    oasis_to_bids.convert_clinical_data(clinical_data_directory, bids_directory)

    compare_folders(join(root, "out"), join(root, "ref"), shared_folder_name="bids")
    clean_folder(join(root, "out", "bids"), recreate=True)


def test_run_Oasis3ToBids(cmdopt):
    from os.path import abspath, dirname, join

    from clinica.iotools.converters.oasis3_to_bids.oasis3_to_bids import Oasis3ToBids

    root = join(dirname(abspath(__file__)), pardir, pardir, "data", "Oasis3ToBids")

    clean_folder(join(root, "out", "bids"), recreate=True)

    # Data location
    dataset_directory = join(root, "in", "unorganized")
    bids_directory = join(root, "out", "bids")
    clinical_data_directory = join(root, "in", "clinical_data")

    oasis_to_bids = Oasis3ToBids()
    oasis_to_bids.convert_images(dataset_directory, bids_directory)
    oasis_to_bids.convert_clinical_data(clinical_data_directory, bids_directory)

    compare_folders(join(root, "out"), join(root, "ref"), shared_folder_name="bids")
    clean_folder(join(root, "out", "bids"), recreate=True)


def test_run_Adni2Bids(cmdopt):
    from os.path import abspath, dirname, join

    from clinica.iotools.converters.adni_to_bids.adni_to_bids import AdniToBids

    root = join(dirname(abspath(__file__)), pardir, pardir, "data", "Adni2Bids")

    clean_folder(join(root, "out", "bids"), recreate=True)

    adni_to_bids = AdniToBids()
    adni_to_bids.check_adni_dependencies()

    dataset_directory = join(root, "in", "unorganized_data")
    clinical_data_directory = join(root, "in", "clinical_data")
    bids_directory = join(root, "out", "bids")
    subjects_list = join(root, "in", "subjects.txt")
    modalities = ["T1", "PET_FDG", "PET_AMYLOID", "PET_TAU", "DWI", "FLAIR", "fMRI"]
    adni_to_bids.convert_images(
        dataset_directory,
        clinical_data_directory,
        bids_directory,
        subjects_list,
        modalities,
    )
    adni_to_bids.convert_clinical_data(clinical_data_directory, bids_directory)
    # Generate tree of output files
    compare_folders(join(root, "out"), join(root, "ref"), shared_folder_name="bids")
    clean_folder(join(root, "out", "bids"), recreate=True)


def test_run_Aibl2Bids(cmdopt):
    from os.path import abspath, dirname, join

    from clinica.iotools.converters.aibl_to_bids.aibl_to_bids import (
        convert_clinical_data,
        convert_images,
    )

    root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
    root = join(root, "data", "Aibl2Bids")

    dataset_directory = join(root, "in", "unorganized_data")
    clinical_data_directory = join(root, "in", "Data_extract_3.2.5")
    bids_directory = join(root, "out", "bids")

    clean_folder(join(root, "out", "bids"), recreate=True)

    # Perform conversion of dataset
    convert_images(dataset_directory, clinical_data_directory, bids_directory)
    convert_clinical_data(bids_directory, clinical_data_directory)

    # Evaluate difference between ref and out
    compare_folders(join(root, "out"), join(root, "ref"), shared_folder_name="bids")
    clean_folder(join(root, "out", "bids"), recreate=True)
