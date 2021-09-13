# coding: utf8

"""
This file contains a set of functional tests designed to check the correct execution of the pipeline and the
different functions available in Clinica
"""


import warnings
from os import PathLike, fspath
from pathlib import Path
from test.nonregression.testing_tools import (
    create_list_hashes,
    identical_subject_list,
    same_missing_modality_tsv,
)

import pytest

# Determine location for working_directory
warnings.filterwarnings("ignore")


@pytest.fixture(
    params=[
        "CreateSubjectSessionList",
        # "CreateMergeFile",
        "ComputeMissingModalities",
        "CenterNifti",
    ]
)
def test_name(request):
    return request.param


def run_createsubjectsession(
    input_dir: PathLike, output_dir: PathLike, ref_dir: PathLike
) -> None:
    from pathlib import Path

    from clinica.iotools.utils import data_handling as dt

    # Arrange
    tsv_name = "subject_session_list.tsv"
    # Act - Create subject_session file
    dt.create_subs_sess_list(str(Path(input_dir) / "bids"), str(output_dir), tsv_name)
    # Assert
    out_tsv = fspath(Path(output_dir) / tsv_name)
    ref_tsv = fspath(Path(ref_dir) / tsv_name)
    assert identical_subject_list(out_tsv, ref_tsv)


def run_createmergefile(
    input_dir: PathLike, output_dir: PathLike, ref_dir: PathLike
) -> None:
    import shutil
    from filecmp import cmp
    from pathlib import Path

    import pandas as pd

    from clinica.iotools.utils import data_handling as dt

    # Arrange
    out_tsv = Path(output_dir) / "output_file.tsv"
    subject_session_tsv = Path(input_dir) / "subjects_sessions.tsv"
    caps_directory = Path(output_dir) / "caps"
    shutil.copytree(Path(input_dir) / "caps", caps_directory, copy_function=shutil.copy)
    # Act
    dt.create_merge_file(
        Path(input_dir) / "bids",
        out_tsv,
        caps_dir=caps_directory,
        tsv_file=subject_session_tsv,
        pipelines=None,
        atlas_selection=None,
        pvc_restriction=None,
        group_selection=None,
    )
    # Assert
    ref_tsv = fspath(Path(ref_dir) / "output_file.tsv")
    out_df = pd.read_csv(out_tsv, sep="/")
    ref_df = pd.read_csv(ref_tsv, sep="/")
    assert out_df.equals(ref_df)
    assert cmp(out_tsv, ref_tsv)


def run_computemissingmodalities(
    input_dir: PathLike, output_dir: PathLike, ref_dir: PathLike
) -> None:
    from pathlib import Path

    from clinica.iotools.utils import data_handling as dt

    output_name = "missing_modalities"
    bids_dir = Path(input_dir) / "bids"

    dt.compute_missing_mods(bids_dir, output_dir, output_name)

    filenames = [
        "missing_modalities_ses-M00.tsv",
        "missing_modalities_ses-M03.tsv",
        "missing_modalities_ses-M06.tsv",
        "missing_modalities_ses-M12.tsv",
        "missing_modalities_ses-M24.tsv",
        "missing_modalities_ses-M48.tsv",
    ]
    for i in range(len(filenames)):
        outname = Path(output_dir) / filenames[i]
        refname = Path(ref_dir) / filenames[i]
        if not outname.exists():
            raise FileNotFoundError(
                "A file called "
                + str(outname)
                + " should have been generated, but it does not exists"
            )
        assert same_missing_modality_tsv(outname, refname)


def run_centernifti(
    input_dir: PathLike, output_dir: PathLike, ref_dir: PathLike
) -> None:
    from pathlib import Path

    from clinica.iotools.utils.data_handling import center_all_nifti

    output_dir = Path(output_dir) / "bids_centered"
    bids_dir = Path(input_dir) / "bids"

    all_modalities = [
        "t1w",
        "pet",
        "dwi",
        "magnitude",
        "bold",
        "flair",
        "t2",
        "phasediff",
    ]

    center_all_nifti(
        fspath(bids_dir),
        fspath(output_dir),
        all_modalities,
        center_all_files=True,
    )
    hashes_out = create_list_hashes(
        fspath(output_dir), extensions_to_keep=(".nii.gz", ".nii")
    )
    hashes_ref = create_list_hashes(
        fspath(Path(ref_dir) / "bids_centered"), extensions_to_keep=(".nii.gz", ".nii")
    )
    assert hashes_out == hashes_ref

    if hashes_out != hashes_ref:
        raise RuntimeError("Hashes of nii* files are different between out and ref")


def test_run_iotools(cmdopt, tmp_path, test_name):

    base_dir = Path(cmdopt["input"])
    input_dir = base_dir / test_name / "in"
    ref_dir = base_dir / test_name / "ref"
    tmp_out_dir = tmp_path / test_name / "out"
    tmp_out_dir.mkdir(parents=True)

    if test_name == "CreateSubjectSessionList":
        run_createsubjectsession(input_dir, tmp_out_dir, ref_dir)

    elif test_name == "CreateMergeFile":
        run_createmergefile(input_dir, tmp_out_dir, ref_dir)

    elif test_name == "ComputeMissingModalities":
        run_computemissingmodalities(input_dir, tmp_out_dir, ref_dir)

    elif test_name == "CenterNifti":
        run_centernifti(input_dir, tmp_out_dir, ref_dir)

    else:
        print(f"Test {test_name} not available.")
        assert 0
