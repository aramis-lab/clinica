"""
This file contains a set of functional tests designed to check the correct
execution of the pipeline and the different functions available in Clinica.
"""


from pathlib import Path

import pytest


@pytest.mark.fast
def test_create_subject_session(cmdopt, tmp_path):
    from test.nonregression.testing_tools import compare_subject_session_tsv

    from clinica.iotools.utils.data_handling import create_subs_sess_list

    base_dir = Path(cmdopt["input"])
    ref_dir = base_dir / "CreateSubjectSessionList" / "ref"
    input_dir = base_dir / "CreateSubjectSessionList" / "in"
    tsv_name = "subject_session_list.tsv"

    create_subs_sess_list(str(input_dir / "bids"), str(tmp_path), tsv_name)

    assert compare_subject_session_tsv(tmp_path / tsv_name, ref_dir / tsv_name)


@pytest.mark.fast
def test_create_merge_file(cmdopt, tmp_path):
    import shutil
    from filecmp import cmp

    import pandas as pd

    from clinica.iotools.utils.data_handling import create_merge_file

    base_dir = Path(cmdopt["input"])
    ref_dir = base_dir / "CreateMergeFile" / "ref"
    input_dir = base_dir / "CreateMergeFile" / "in"
    out_tsv = tmp_path / "output_file.tsv"

    shutil.copytree(input_dir / "caps", tmp_path / "caps", copy_function=shutil.copy)

    create_merge_file(
        input_dir / "bids",
        out_tsv,
        caps_dir=tmp_path / "caps",
        tsv_file=input_dir / "subjects_sessions.tsv",
        pipelines=None,
        atlas_selection=None,
        pvc_restriction=None,
        group_selection=None,
    )

    ref_tsv = ref_dir / "output_file.tsv"
    out_df = pd.read_csv(out_tsv, sep="\t")
    ref_df = pd.read_csv(ref_tsv, sep="\t")

    assert out_df.equals(ref_df)
    assert cmp(out_tsv, ref_tsv)


@pytest.mark.fast
def test_compute_missing_modalities(cmdopt, tmp_path):
    from test.nonregression.testing_tools import compare_missing_modality_tsv

    from clinica.iotools.utils.data_handling import compute_missing_mods

    base_dir = Path(cmdopt["input"])
    ref_dir = base_dir / "ComputeMissingModalities" / "ref"
    bids_dir = base_dir / "ComputeMissingModalities" / "in" / "bids"

    compute_missing_mods(bids_dir, tmp_path, "missing_modalities")

    for filename in (
        f"missing_modalities_ses-M0{i}.tsv"
        for i in ("00", "03", "06", "12", "24", "48")
    ):
        out_name = tmp_path / filename
        ref_name = ref_dir / filename
        if not out_name.exists():
            raise FileNotFoundError(
                f"A file called {out_name} should have been generated, but it does not exists"
            )
        assert compare_missing_modality_tsv(out_name, ref_name)


@pytest.mark.fast
def test_center_nifti(cmdopt, tmp_path):
    from test.nonregression.testing_tools import (
        compare_niftis,
        compare_txt_files,
        create_list_hashes,
    )

    from clinica.iotools.utils import center_nifti

    base_dir = Path(cmdopt["input"])
    output_dir = tmp_path / "bids_centered"
    ref_dir = base_dir / "CenterNifti" / "ref"

    center_nifti(
        str(base_dir / "CenterNifti" / "in" / "bids"),
        output_dir,
        centering_threshold=0,
    )
    hashes_out = create_list_hashes(output_dir, extensions_to_keep=(".nii.gz", ".nii"))
    hashes_ref = create_list_hashes(
        ref_dir / "bids_centered", extensions_to_keep=(".nii.gz", ".nii")
    )

    assert hashes_out == hashes_ref

    if hashes_out != hashes_ref:
        raise RuntimeError("Hashes of nii* files are different between out and ref")

    compare_niftis(output_dir, ref_dir / "bids_centered")
    compare_txt_files(output_dir, ref_dir / "bids_centered")
