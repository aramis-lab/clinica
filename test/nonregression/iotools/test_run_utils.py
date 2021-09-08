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
        "CreateMergeFile",
        "ComputeMissingModalities",
        "CenterNifti",
    ]
)
def name(request):
    return request.param


def test_run_iotools(cmdopt, tmp_path, name):
    import shutil

    base_dir = Path(cmdopt["input"])
    input_dir = base_dir / name / "in"
    ref_dir = base_dir / name / "ref"
    tmp_out_dir = tmp_path / name / "out"
    tmp_out_dir.mkdir(parents=True)
    bids_directory = input_dir / "bids"

    if name == "CreateSubjectSessionList":
        from clinica.iotools.utils import data_handling as dt

        # Arrage
        tsv_name = "subject_session_list.tsv"
        # Act - Create subject_session file
        dt.create_subs_sess_list(fspath(bids_directory), fspath(tmp_out_dir), tsv_name)
        # Assert
        out_tsv = fspath(tmp_out_dir / tsv_name)
        ref_tsv = fspath(ref_dir / tsv_name)
        assert identical_subject_list(out_tsv, ref_tsv)
    elif name == "CreateMergeFile":
        from filecmp import cmp
        from os import remove
        from os.path import abspath, dirname, join

        from clinica.iotools.utils import data_handling as dt

        # Arrange
        out_tsv = tmp_out_dir / "output_file.tsv"
        subject_session_tsv = input_dir / "subjects_sessions.tsv"
        caps_directory = tmp_out_dir / "caps"
        shutil.copytree(input_dir / "caps", caps_directory, copy_function=shutil.copy)
        # Act
        dt.create_merge_file(
            bids_directory,
            out_tsv,
            caps_dir=caps_directory,
            tsv_file=subject_session_tsv,
            pipelines=None,
            atlas_selection=None,
            pvc_restriction=None,
            group_selection=None,
        )
        # Assert
        ref_tsv = fspath(ref_dir / "output_file.tsv")
        assert cmp(out_tsv, ref_tsv)
    elif name == "ComputeMissingModalities":
        from clinica.iotools.utils import data_handling as dt

        output_name = "missing_modalities"

        dt.compute_missing_mods(
            fspath(bids_directory), fspath(tmp_out_dir), fspath(output_name)
        )

        filenames = [
            "missing_modalities_ses-M00.tsv",
            "missing_modalities_ses-M03.tsv",
            "missing_modalities_ses-M06.tsv",
            "missing_modalities_ses-M12.tsv",
            "missing_modalities_ses-M24.tsv",
            "missing_modalities_ses-M48.tsv",
        ]
        for i in range(len(filenames)):
            outname = tmp_out_dir / filenames[i]
            refname = ref_dir / filenames[i]
            if not outname.exists():
                raise FileNotFoundError(
                    "A file called "
                    + outname
                    + " should have been generated, but it does not exists"
                )
            assert same_missing_modality_tsv(fspath(outname), fspath(refname))

    elif name == "CenterNifti":
        from clinica.iotools.utils.data_handling import center_all_nifti

        output_dir = tmp_out_dir / "bids_centered"

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
            fspath(bids_directory),
            fspath(output_dir),
            all_modalities,
            center_all_files=True,
        )
        hashes_out = create_list_hashes(
            fspath(output_dir), extensions_to_keep=(".nii.gz", ".nii")
        )
        hashes_ref = create_list_hashes(
            fspath(ref_dir / "bids_centered"), extensions_to_keep=(".nii.gz", ".nii")
        )
        assert hashes_out == hashes_ref

        if hashes_out != hashes_ref:
            raise RuntimeError("Hashes of nii* files are different between out and ref")
    else:
        print(f"Test {name} not available.")
        pass
        # assert 0


def test_run_CreateMergeFile(cmdopt):
    import shutil
    from filecmp import cmp
    from os import remove
    from os.path import abspath, dirname, join

    from clinica.iotools.utils import data_handling as dt

    root = join(dirname(abspath(__file__)), pardir, pardir, "data", "CreateMergeFile")

    bids_directory = join(root, "in", "bids")
    out_tsv = join(root, "out", "output_file.tsv")
    subject_session_tsv = join(root, "in", "subjects_sessions.tsv")

    clean_folder(join(root, "out", "caps"), recreate=False)
    shutil.copytree(join(root, "in", "caps"), join(root, "out", "caps"))
    caps_directory = join(root, "out", "caps")

    dt.create_merge_file(
        bids_directory,
        out_tsv,
        caps_dir=caps_directory,
        tsv_file=subject_session_tsv,
        pipelines=None,
        atlas_selection=None,
        pvc_restriction=None,
        group_selection=None,
    )
    # Comparison step
    ref_tsv = join(root, "ref", "output_file.tsv")
    assert cmp(out_tsv, ref_tsv)
    remove(out_tsv)
    clean_folder(join(root, "out", "caps"), recreate=False)


def test_run_ComputeMissingModalities(cmdopt):
    from os import remove
    from os.path import abspath, dirname, exists, join

    from clinica.iotools.utils import data_handling as dt

    root = join(dirname(abspath(__file__)), pardir, pardir, "data", "ComputeMissingMod")

    bids_directory = join(root, "in", "bids")
    output_directory = join(root, "out")
    output_name = "missing_modalities"

    dt.compute_missing_mods(bids_directory, output_directory, output_name)

    filenames = [
        "missing_modalities_ses-M00.tsv",
        "missing_modalities_ses-M03.tsv",
        "missing_modalities_ses-M06.tsv",
        "missing_modalities_ses-M12.tsv",
        "missing_modalities_ses-M24.tsv",
        "missing_modalities_ses-M48.tsv",
    ]
    for i in range(len(filenames)):
        outname = join(root, "out", filenames[i])
        refname = join(root, "ref", filenames[i])
        if not exists(outname):
            raise FileNotFoundError(
                "A file called "
                + outname
                + " should have been generated, but it does not exists"
            )
        assert same_missing_modality_tsv(outname, refname)
        remove(join(root, "out", filenames[i]))


def test_run_CenterNifti(cmdopt):
    from os.path import abspath, dirname, join

    from clinica.iotools.utils.data_handling import center_all_nifti

    root = dirname(abspath(join(abspath(__file__), pardir, pardir)))
    root = join(root, "data", "CenterNifti")

    clean_folder(join(root, "out", "bids_centered"), recreate=True)

    bids_dir = join(root, "in", "bids")
    output_dir = join(root, "out", "bids_centered")

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

    center_all_nifti(bids_dir, output_dir, all_modalities, center_all_files=True)
    hashes_out = create_list_hashes(output_dir, extensions_to_keep=(".nii.gz", ".nii"))
    hashes_ref = create_list_hashes(
        join(root, "ref", "bids_centered"), extensions_to_keep=(".nii.gz", ".nii")
    )

    if hashes_out != hashes_ref:
        raise RuntimeError("Hashes of nii* files are different between out and ref")
    clean_folder(join(root, "out", "bids_centered"), recreate=False)
