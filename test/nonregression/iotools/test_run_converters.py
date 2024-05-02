"""
This file contains a set of functional tests designed to check the correct execution of the pipeline and the
different functions available in Clinica
"""
from pathlib import Path
from test.nonregression.testing_tools import compare_folders, configure_paths


def test_run_nifd_to_bids(cmdopt, tmp_path):
    from clinica.iotools.converters.nifd_to_bids.nifd_to_bids import convert

    base_dir = Path(cmdopt["input"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "Nifd2Bids")
    output_dir = tmp_path / "bids"

    convert(
        path_to_dataset=input_dir / "unorganized",
        bids_dir=output_dir,
        path_to_clinical=input_dir / "clinical_data",
    )
    compare_folders(output_dir, ref_dir, output_dir)


def test_run_oasis_to_bids(cmdopt, tmp_path):
    from clinica.iotools.converters.oasis_to_bids.oasis_to_bids import OasisToBids

    base_dir = Path(cmdopt["input"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "Oasis2Bids")
    output_dir = tmp_path / "bids"

    OasisToBids().convert(
        input_dir / "unorganized", output_dir, input_dir / "clinical_data"
    )

    compare_folders(output_dir, ref_dir / "bids", output_dir)


def test_run_oasis3_to_bids(cmdopt, tmp_path):
    from clinica.iotools.converters.oasis3_to_bids.oasis3_to_bids import convert

    base_dir = Path(cmdopt["input"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "Oasis3ToBids")
    output_dir = tmp_path / "bids"

    convert(
        input_dir / "unorganized",
        output_dir / "bids",
        input_dir / "clinical_data",
    )
    compare_folders(output_dir / "bids", ref_dir / "bids", output_dir)


def test_run_adni_to_bids(cmdopt, tmp_path):
    from clinica.iotools.converters.adni_to_bids.adni_to_bids import AdniToBids

    base_dir = Path(cmdopt["input"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "Adni2Bids")
    output_dir = tmp_path / "bids"
    clinical_data_directory = input_dir / "clinical_data"
    xml_directory = input_dir / "xml_metadata"
    dataset_directory = input_dir / "unorganized"
    subjects_list = input_dir / "subjects.txt"
    modalities = [
        "T1",
        "PET_FDG",
        "PET_AMYLOID",
        "PET_TAU",
        "DWI",
        "FLAIR",
        "fMRI",
        "FMAP",
    ]

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

    compare_folders(output_dir / "bids", ref_dir / "bids", output_dir)


def test_run_aibl_to_bids(cmdopt, tmp_path):
    from clinica.iotools.converters.aibl_to_bids.aibl_to_bids import convert

    base_dir = Path(cmdopt["input"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "Aibl2Bids")
    output_dir = tmp_path / "bids"

    convert(
        input_dir / "unorganized_data", input_dir / "Data_extract_3.2.5", output_dir
    )

    compare_folders(output_dir, ref_dir / "bids", tmp_path)


def test_run_habs_to_bids(cmdopt, tmp_path):
    from clinica.iotools.converters.habs_to_bids.habs_to_bids import convert

    base_dir = Path(cmdopt["input"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "HabsToBids")
    output_dir = tmp_path / "bids"

    convert(input_dir / "unorganized", output_dir)
    compare_folders(output_dir, ref_dir / "bids", output_dir)


def test_run_ukb_to_bids(cmdopt, tmp_path):
    from clinica.iotools.converters.ukb_to_bids.ukb_to_bids import convert

    base_dir = Path(cmdopt["input"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "UkbToBids")
    output_dir = tmp_path / "bids"

    convert(
        input_dir / "unorganized",
        output_dir / "bids",
        input_dir / "clinical_data",
    )
    compare_folders(output_dir / "bids", ref_dir / "bids", output_dir)


def test_run_genfi_to_bids(cmdopt, tmp_path):
    from clinica.iotools.converters.genfi_to_bids.genfi_to_bids import convert

    base_dir = Path(cmdopt["input"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "GenfiToBids")
    output_dir = tmp_path / "bids"

    convert(
        input_dir / "unorganized",
        output_dir,
        path_to_clinical=None,
        gif=False,
        path_to_clinical_tsv=None,
    )
    compare_folders(output_dir, ref_dir / "bids", output_dir)
