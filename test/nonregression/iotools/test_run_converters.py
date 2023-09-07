"""
This file contains a set of functional tests designed to check the correct execution of the pipeline and the
different functions available in Clinica
"""
from pathlib import Path
from test.nonregression.testing_tools import compare_folders, configure_paths


def test_run_nifd_to_bids(cmdopt, tmp_path):
    from clinica.iotools.converters.nifd_to_bids.nifd_to_bids import convert_images

    base_dir = Path(cmdopt["input"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "Nifd2Bids")
    output_dir = tmp_path / "bids"

    convert_images(
        path_to_clinical=input_dir / "clinical_data",
        path_to_dataset=input_dir / "unorganized",
        bids_dir=output_dir,
    )

    compare_folders(output_dir, ref_dir, output_dir)


def test_run_oasis_to_bids(cmdopt, tmp_path):
    from clinica.iotools.converters.oasis_to_bids.oasis_to_bids import (
        OasisToBidsConverter,
    )

    base_dir = Path(cmdopt["input"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "Oasis2Bids")

    OasisToBidsConverter(
        input_dir / "unorganized",
        tmp_path / "bids",
        input_dir / "clinical_data",
    ).convert()

    compare_folders(tmp_path / "bids", ref_dir / "bids", tmp_path / "bids")


def test_run_oasis3_to_bids(cmdopt, tmp_path):
    from clinica.iotools.converters.oasis3_to_bids.oasis3_to_bids import (
        Oasis3ToBidsConverter,
    )

    base_dir = Path(cmdopt["input"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "Oasis3ToBids")
    output_dir = tmp_path / "bids"
    clinical_data_directory = input_dir / "clinical_data"

    Oasis3ToBidsConverter(
        input_dir / "unorganized",
        output_dir / "bids",
        clinical_data_directory,
    ).convert()

    compare_folders(output_dir / "bids", ref_dir / "bids", output_dir)


def test_run_adni_to_bids(cmdopt, tmp_path):
    from clinica.iotools.converters.adni_to_bids.adni_to_bids import AdniToBidsConverter

    base_dir = Path(cmdopt["input"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "Adni2Bids")
    output_dir = tmp_path / "bids"

    AdniToBidsConverter(
        input_dir / "unorganized_data",
        output_dir / "bids",
        input_dir / "clinical_data",
        input_dir / "xml_metadata",
    ).convert(input_dir / "subjects.txt")

    compare_folders(output_dir / "bids", ref_dir / "bids", output_dir)


def test_run_aibl_to_bids(cmdopt, tmp_path):
    from pathlib import Path

    from clinica.iotools.converters.aibl_to_bids.aibl_to_bids import AiblToBidsConverter

    base_dir = Path(cmdopt["input"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "Aibl2Bids")

    AiblToBidsConverter(
        input_dir / "unorganized_data",
        tmp_path / "bids",
        input_dir / "Data_extract_3.2.5",
    ).convert()

    compare_folders(tmp_path / "bids", ref_dir / "bids", tmp_path)


def test_run_habs_to_bids(cmdopt, tmp_path):
    from click.testing import CliRunner

    from clinica.iotools.converters.habs_to_bids.habs_to_bids_cli import cli

    base_dir = Path(cmdopt["input"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "HabsToBids")
    output_dir = tmp_path / "bids"

    runner = CliRunner()
    result = runner.invoke(cli, [str(input_dir), str(output_dir)])

    assert result.exit_code == 0
    compare_folders(output_dir, ref_dir, output_dir)


def test_run_ukb_to_bids(cmdopt, tmp_path):
    from clinica.iotools.converters.ukb_to_bids.ukb_to_bids import UkbToBidsConverter

    base_dir = Path(cmdopt["input"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "UkbToBids")
    output_dir = tmp_path / "bids"

    UkbToBidsConverter(
        input_dir / "unorganized",
        output_dir / "bids",
        input_dir / "clinical_data",
    ).convert()

    compare_folders(output_dir / "bids", ref_dir / "bids", output_dir)


def test_run_genfi_to_bids(cmdopt, tmp_path):
    from clinica.iotools.converters.genfi_to_bids.genfi_to_bids import (
        GenfiToBidsConverter,
    )

    base_dir = Path(cmdopt["input"])
    input_dir, tmp_dir, ref_dir = configure_paths(base_dir, tmp_path, "GenfiToBids")

    GenfiToBidsConverter(
        input_dir / "unorganized",
        tmp_path / "bids",
        tmp_path,
        gif=False,
        clinical_data=False,
    ).convert()

    compare_folders(tmp_path / "bids", ref_dir / "bids", tmp_path / "bids")
