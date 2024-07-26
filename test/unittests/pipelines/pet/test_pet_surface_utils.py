import os
import re
from pathlib import Path
from unittest.mock import patch

import nibabel as nib
import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_array_almost_equal, assert_array_equal

from clinica.utils.image import HemiSphere
from clinica.utils.pet import SUVRReferenceRegion, Tracer


def test_normalize_suvr_error(tmp_path):
    from clinica.pipelines.pet.surface.utils import normalize_suvr
    from clinica.utils.exceptions import ClinicaImageError

    nib.save(
        nib.Nifti1Image(np.random.random((10, 10, 10)), np.eye(4)),
        tmp_path / "pet.nii.gz",
    )
    nib.save(
        nib.Nifti1Image(np.zeros((10, 10, 10)), np.eye(4)), tmp_path / "mask.nii.gz"
    )

    with pytest.raises(
        ClinicaImageError,
        match=(
            f"The eroded mask located at {tmp_path / 'mask.nii.gz'} contains only zero values. "
            "A problem likely occurred when moving the eroded mask from MNI to gtmsegspace."
        ),
    ):
        normalize_suvr(tmp_path / "pet.nii.gz", tmp_path / "mask.nii.gz")


def test_normalize_suvr(tmp_path):
    from clinica.pipelines.pet.surface.utils import normalize_suvr

    nib.save(
        nib.Nifti1Image(0.5 * np.ones((20, 20, 20)), np.eye(4)), tmp_path / "pet.nii.gz"
    )

    mask = np.zeros((20, 20, 20))
    mask[5:15, 5:15, 5:15] = 1
    nib.save(nib.Nifti1Image(mask, np.eye(4)), tmp_path / "mask.nii.gz")

    normalize_suvr(tmp_path / "pet.nii.gz", tmp_path / "mask.nii.gz", tmp_path)

    assert (tmp_path / "suvr_pet.nii.gz").exists()
    suvr = nib.load(tmp_path / "suvr_pet.nii.gz")
    assert_array_equal(suvr.affine, np.eye(4))
    assert_array_equal(suvr.get_fdata(), np.ones((20, 20, 20)))


@pytest.mark.parametrize(
    "platform,prefix",
    [
        ("linux", ""),
        ("darwin", "export DYLD_LIBRARY_PATH=$FREESURFER_HOME/lib/gcc/lib && "),
    ],
)
def test_build_mri_expand_command(tmp_path, mocker, platform, prefix):
    from clinica.pipelines.pet.surface.utils import _build_mri_expand_command

    mocker.patch("platform.system", return_value=platform)
    cmd = _build_mri_expand_command(tmp_path / "lh.white", tmp_path / "lh.white_exp-")

    assert (
        cmd
        == f"{prefix}mris_expand -thickness -N 13 {tmp_path / 'lh.white'} 0.65 {tmp_path / 'lh.white_exp-'}"
    )


def run_mri_expand_as_subprocess_patch(surface: Path, out_file: Path):
    """Function used to mock clinica.pipelines.pet.surface.utils._run_mri_expand_as_subprocess.

    It simply writes the files that are expected to be written by the function without relying on mri_expand.
    """
    for file in (
        out_file.parent / f"{out_file.name}0{n}"
        for n in (
            "00",
            "01",
            "02",
            "03",
            "04",
            "05",
            "06",
            "07",
            "08",
            "09",
            "10",
            "11",
            "12",
            "13",
        )
    ):
        file.touch()


def test_run_mris_expand(tmp_path):
    from unittest.mock import patch

    from clinica.pipelines.pet.surface.utils import run_mris_expand

    with patch(
        "clinica.pipelines.pet.surface.utils._run_mri_expand_as_subprocess",
        wraps=run_mri_expand_as_subprocess_patch,
    ) as mock:
        run_mris_expand(tmp_path / "lh.white", output_dir=tmp_path / "output")

        mock.assert_called_once_with(tmp_path / "lh.white", tmp_path / "lh.white_exp-")
        assert (tmp_path / "output").exists()
        assert set([f.name for f in (tmp_path / "output").iterdir()]) == set(
            [f"lh.white_exp-0{x}" for x in ("07", "08", "09", "10", "11", "12", "13")]
        )


def run_mri_expand_as_subprocess_patch_broken(surface: Path, out_file: Path):
    """Function used to mock clinica.pipelines.pet.surface.utils._run_mri_expand_as_subprocess in order to make it fail.

    It does not generate one of the expected files.
    """
    for file in (
        out_file.parent / f"{out_file.name}0{n}"
        for n in (
            "00",
            "01",
            "02",
            "03",
            "04",
            "05",
            "06",
            "08",
            "09",
            "10",
            "11",
            "12",
            "13",
        )
    ):
        file.touch()


def test_run_mris_expand_error(tmp_path):
    from clinica.pipelines.pet.surface.utils import run_mris_expand

    with patch(
        "clinica.pipelines.pet.surface.utils._run_mri_expand_as_subprocess",
        wraps=run_mri_expand_as_subprocess_patch_broken,
    ):
        with pytest.raises(
            FileNotFoundError,
            match=(
                f"File {tmp_path / 'lh.white_exp-007'} is missing. "
                "Something wrong might have occurred prior to this step."
            ),
        ):
            run_mris_expand(tmp_path / "lh.white", output_dir=tmp_path / "output")


def test_get_new_subjects_directory(tmp_path):
    from clinica.pipelines.pet.surface.utils import _get_new_subjects_directory

    subject_dir, freesurfer_id = _get_new_subjects_directory(
        tmp_path / "caps", "sub-01", "ses-M006"
    )

    assert (
        subject_dir
        == tmp_path
        / "caps"
        / "subjects"
        / "sub-01"
        / "ses-M006"
        / "t1"
        / "freesurfer_cross_sectional"
    )
    assert freesurfer_id == "sub-01_ses-M006"


def test_get_new_subjects_directory_longitudinal_no_long_folder_error(tmp_path):
    from clinica.pipelines.pet.surface.utils import (
        _get_new_subjects_directory_longitudinal,
    )
    from clinica.utils.exceptions import ClinicaCAPSError

    folder = tmp_path / "caps" / "subjects" / "sub-01" / "ses-M006" / "t1"
    folder.mkdir(parents=True)

    with pytest.raises(
        ClinicaCAPSError,
        match=re.escape(
            f"Folder {folder} does not contains a folder labeled long-*. "
            "Have you run t1-freesurfer-longitudinal?",
        ),
    ):
        _get_new_subjects_directory_longitudinal(
            tmp_path / "caps", "sub-01", "ses-M006"
        )


def test_get_new_subjects_directory_longitudinal_multiple_long_folder_error(tmp_path):
    from clinica.pipelines.pet.surface.utils import (
        _get_new_subjects_directory_longitudinal,
    )
    from clinica.utils.exceptions import ClinicaCAPSError

    folder = tmp_path / "caps" / "subjects" / "sub-01" / "ses-M006" / "t1"
    folder.mkdir(parents=True)
    (folder / "long-foo").mkdir()
    (folder / "long-bar").mkdir()

    with pytest.raises(
        ClinicaCAPSError,
        match=re.escape(
            f"Folder {folder} contains 2 folders labeled long-*. Only 1 can exist."
        ),
    ):
        _get_new_subjects_directory_longitudinal(
            tmp_path / "caps", "sub-01", "ses-M006"
        )


def test_get_new_subjects_directory_longitudinal(tmp_path):
    from clinica.pipelines.pet.surface.utils import (
        _get_new_subjects_directory_longitudinal,
    )

    folder = tmp_path / "caps" / "subjects" / "sub-01" / "ses-M006" / "t1"
    folder.mkdir(parents=True)
    (folder / "long-foo").mkdir()

    subject_directory, freesurfer_id = _get_new_subjects_directory_longitudinal(
        tmp_path / "caps", "sub-01", "ses-M006"
    )

    assert (
        subject_directory
        == tmp_path
        / "caps"
        / "subjects"
        / "sub-01"
        / "ses-M006"
        / "t1"
        / "long-foo"
        / "freesurfer_longitudinal"
    )
    assert freesurfer_id == "sub-01_ses-M006.long.sub-01_long-foo"


@pytest.mark.parametrize("hemisphere", HemiSphere)
@pytest.mark.parametrize(
    "platform,prefix",
    [
        ("linux", ""),
        ("darwin", "export DYLD_LIBRARY_PATH=$FREESURFER_HOME/lib/gcc/lib && "),
    ],
)
def test_build_mri_surf2surf_command(tmp_path, mocker, hemisphere, platform, prefix):
    from clinica.pipelines.pet.surface.utils import _build_mri_surf2surf_command

    mocker.patch("platform.system", return_value=platform)

    command = _build_mri_surf2surf_command(
        tmp_path / f"{hemisphere.value}.surface",
        tmp_path / "registration",
        tmp_path / "gtmseg",
        "sub-01_ses-M006.long.sub-01_long-foo",
        tmp_path / "output",
    )

    assert command == (
        f"{prefix}mri_surf2surf --reg {tmp_path / 'registration'} {tmp_path / 'gtmseg'} "
        f"--sval-xyz surface --hemi {hemisphere.value} --tval-xyz {tmp_path / 'gtmseg'} "
        f"--tval {tmp_path / 'output'} --s sub-01_ses-M006.long.sub-01_long-foo"
    )


def run_mri_surf2surf_as_subprocess_patch(
    surface: Path,
    registration: Path,
    gtmsegfile: Path,
    freesurfer_id: str,
    output_file: Path,
):
    output_file.touch()


def setup_folders(tmp_path: Path, hemisphere: HemiSphere, is_longitudinal: bool):
    os.environ["SUBJECTS_DIR"] = str(tmp_path / "old_subject_directory")
    folder = tmp_path / "caps" / "subjects" / "sub-01" / "ses-M006" / "t1"
    folder.mkdir(parents=True)
    if is_longitudinal:
        surf_folder = (
            folder
            / "long-foo"
            / "freesurfer_longitudinal"
            / "sub-01_ses-M006.long.sub-01_long-foo"
            / "surf"
        )
    else:
        surf_folder = folder / "freesurfer_cross_sectional" / "sub-01_ses-M006" / "surf"
    surf_folder.mkdir(parents=True)
    surface = tmp_path / f"{hemisphere.value}.white"
    surface.touch()
    output_dir = tmp_path / "output"
    output_dir.mkdir()

    return surface, surf_folder, output_dir


@pytest.mark.parametrize("is_longitudinal", (True, False))
@pytest.mark.parametrize("hemisphere", HemiSphere)
def test_run_mri_surf2surf(tmp_path, is_longitudinal, hemisphere):
    from clinica.pipelines.pet.surface.utils import run_mri_surf2surf

    surface, surf_folder, output_dir = setup_folders(
        tmp_path, hemisphere, is_longitudinal
    )

    with patch(
        "clinica.pipelines.pet.surface.utils._run_mri_surf2surf_as_subprocess",
        wraps=run_mri_surf2surf_as_subprocess_patch,
    ) as mock:
        run_mri_surf2surf(
            surface,
            tmp_path / "registration",
            tmp_path / "gtmseg",
            "sub-01",
            "ses-M006",
            tmp_path / "caps",
            is_longitudinal,
            output_dir,
        )
        mock.assert_called_once_with(
            surface,
            tmp_path / "registration",
            tmp_path / "gtmseg",
            "sub-01_ses-M006.long.sub-01_long-foo"
            if is_longitudinal
            else "sub-01_ses-M006",
            output_dir / f"{hemisphere.value}.white_gtmsegspace",
        )
        assert (output_dir / f"{hemisphere.value}.white_gtmsegspace").exists()
        assert not (surf_folder / f"{hemisphere.value}.white").exists()
    assert os.environ["SUBJECTS_DIR"] == str(tmp_path / "old_subject_directory")


@pytest.mark.parametrize("hemisphere", HemiSphere)
@pytest.mark.parametrize(
    "platform,prefix",
    [
        ("linux", ""),
        ("darwin", "export DYLD_LIBRARY_PATH=$FREESURFER_HOME/lib/gcc/lib && "),
    ],
)
def test_build_mri_vol2surf_command(tmp_path, mocker, hemisphere, platform, prefix):
    from clinica.pipelines.pet.surface.utils import _build_mri_vol2surf_command

    mocker.patch("platform.system", return_value=platform)

    command = _build_mri_vol2surf_command(
        tmp_path / "sub-01_ses-M006_pet.nii.gz",
        tmp_path / f"{hemisphere.value}.white",
        "sub-01_ses-M006.long.sub-01_long-foo",
        tmp_path / "output",
    )

    assert command == (
        f"{prefix}mri_vol2surf --mov {tmp_path / 'sub-01_ses-M006_pet.nii.gz'} "
        f"--o {tmp_path / 'output'} --surf white --hemi {hemisphere.value} "
        f"--regheader sub-01_ses-M006.long.sub-01_long-foo --ref gtmseg.mgz --interp nearest"
    )


def _run_mri_vol2surf_as_subprocess_patch(
    pet_volume: Path,
    surface: Path,
    freesurfer_id: str,
    output_file: Path,
):
    output_file.touch()


@pytest.mark.parametrize("is_longitudinal", (True, False))
@pytest.mark.parametrize("hemisphere", HemiSphere)
def test_run_mri_vol2surf(tmp_path, is_longitudinal, hemisphere):
    from clinica.pipelines.pet.surface.utils import run_mri_vol2surf

    surface, surf_folder, output_dir = setup_folders(
        tmp_path, hemisphere, is_longitudinal
    )
    (tmp_path / "gtmseg").touch()
    (surf_folder.parent / "mri").mkdir()

    with patch(
        "clinica.pipelines.pet.surface.utils._run_mri_vol2surf_as_subprocess",
        wraps=_run_mri_vol2surf_as_subprocess_patch,
    ) as mock:
        run_mri_vol2surf(
            tmp_path / "bids" / "sub-01_ses-M006_pet.nii.gz",
            surface,
            "sub-01",
            "ses-M006",
            tmp_path / "caps",
            tmp_path / "gtmseg",
            is_longitudinal,
            output_dir,
        )
        mock.assert_called_once_with(
            tmp_path / "bids" / "sub-01_ses-M006_pet.nii.gz",
            surface,
            "sub-01_ses-M006.long.sub-01_long-foo"
            if is_longitudinal
            else "sub-01_ses-M006",
            output_dir / f"{hemisphere.value}.projection_{surface.name}.mgh",
        )
        assert (
            output_dir / f"{hemisphere.value}.projection_{surface.name}.mgh"
        ).exists()
        assert not (surf_folder / f"{hemisphere.value}.white").exists()
        assert not (surf_folder.parent / "mri" / "gtmseg.mgz").exists()
    assert os.environ["SUBJECTS_DIR"] == str(tmp_path / "old_subject_directory")


@pytest.mark.parametrize("surfaces", [(), ("foo",), ("foo", "bar"), range(6)])
def test_compute_weighted_mean_surface_error(tmp_path, surfaces):
    from clinica.pipelines.pet.surface.utils import compute_weighted_mean_surface

    with pytest.raises(
        ValueError,
        match=(
            "There should be 7 surfaces at this point of the pipeline. "
            f"However 'compute_weighted_mean_surface' received {len(surfaces)} surfaces. "
            "Something probably went wrong in prior steps of the pipeline."
        ),
    ):
        compute_weighted_mean_surface(surfaces, tmp_path / "output")


def test_get_coefficient_for_normal_repartition():
    from clinica.pipelines.pet.surface.utils import (
        _get_coefficient_for_normal_repartition,
    )

    assert _get_coefficient_for_normal_repartition() == (
        0.1034,
        0.1399,
        0.1677,
        0.1782,
        0.1677,
        0.1399,
        0.1034,
    )


@pytest.mark.parametrize("hemisphere", HemiSphere)
def test_compute_weighted_mean_surface(tmp_path, hemisphere):
    from clinica.pipelines.pet.surface.utils import (
        _get_coefficient_for_normal_repartition,
        compute_weighted_mean_surface,
    )

    (tmp_path / "output").mkdir()
    test_image = nib.MGHImage(np.ones((10, 10, 10), dtype="float32"), np.eye(4))
    for i in range(7):
        nib.save(test_image, tmp_path / f"{hemisphere.value}.surf_{i}.mgh")

    average_image = compute_weighted_mean_surface(
        [tmp_path / f"{hemisphere.value}.surf_{i}.mgh" for i in range(7)],
        tmp_path / "output",
    )

    assert (
        average_image
        == tmp_path
        / "output"
        / f"{hemisphere.value}.averaged_projection_on_cortical_surface.mgh"
    )
    average_image = nib.load(average_image)
    assert_array_equal(average_image.affine, np.eye(4))
    assert_array_almost_equal(
        average_image.get_fdata(),
        sum(_get_coefficient_for_normal_repartition()) * np.ones((10, 10, 10)),
    )


def _run_mris_preproc_as_standalone_nipype_node_patch(
    projection: Path,
    freesurfer_id: str,
    fwhm: float,
    output_file: Path,
):
    output_file.touch()


@pytest.mark.parametrize("is_longitudinal", (True, False))
@pytest.mark.parametrize("hemisphere", HemiSphere)
def test_project_onto_fsaverage(tmp_path, is_longitudinal, hemisphere):
    from clinica.pipelines.pet.surface.utils import project_onto_fsaverage

    projection, surf_folder, output_dir = setup_folders(
        tmp_path, hemisphere, is_longitudinal
    )
    (tmp_path / "old_subject_directory" / "fsaverage").mkdir(parents=True)

    with patch(
        "clinica.pipelines.pet.surface.utils._run_mris_preproc_as_standalone_nipype_node",
        wraps=_run_mris_preproc_as_standalone_nipype_node_patch,
    ) as mock:
        project_onto_fsaverage(
            projection,
            "sub-01",
            "ses-M006",
            tmp_path / "caps",
            6,
            is_longitudinal,
            output_dir,
        )
        mock.assert_called_once_with(
            projection,
            "sub-01_ses-M006.long.sub-01_long-foo"
            if is_longitudinal
            else "sub-01_ses-M006",
            6,
            output_dir / f"fsaverage_fwhm-6_{projection.name}",
        )
        assert (output_dir / f"fsaverage_fwhm-6_{projection.name}").exists()
        assert not (surf_folder / projection.name).exists()
    assert os.environ["SUBJECTS_DIR"] == str(tmp_path / "old_subject_directory")


def test_compute_average_pet_signal_based_on_annotations_error(tmp_path):
    from clinica.pipelines.pet.surface.utils import (
        compute_average_pet_signal_based_on_annotations,
    )

    with pytest.raises(
        ValueError,
        match=(
            "The compute_average_pet_signal_based_on_annotations function requires two files "
            "for the argument 'pet_projections', one for the left hemisphere, one for the right. "
            "The following 3 were received:"
        ),
    ):
        compute_average_pet_signal_based_on_annotations(
            (
                tmp_path / "left_projection.mgh",
                tmp_path / "right_projection.mgh",
                tmp_path / "fooo.mgh",
            ),  # noqa
            {},
            tmp_path / "output",
        )


def test_compute_average_pet_signal_based_on_annotations(tmp_path):
    from clinica.pipelines.pet.surface.utils import (
        compute_average_pet_signal_based_on_annotations,
    )
    from clinica.pipelines.utils import FreeSurferAnnotationImage

    data = np.zeros((20, 20, 20), dtype="float32")
    data[5:15, 5:15, 5:15] = 1.0
    # Left hemisphere has 2 values inside the "brain"
    nib.save(
        nib.MGHImage(2 * data.flatten(), np.eye(4)),
        tmp_path / "left_projection.mgh",
    )
    # Right hemisphere has 1 values inside the "brain"
    nib.save(
        nib.MGHImage(data.flatten(), np.eye(4)),
        tmp_path / "right_projection.mgh",
    )
    (tmp_path / "output").mkdir()

    for parcellation in ("destrieux", "desikan"):
        for hemi in HemiSphere:
            nib.freesurfer.io.write_annot(
                tmp_path / f"{parcellation}_{hemi.value}.annot",
                labels=data.flatten().astype("int"),
                ctab=np.array([[25, 25, 25, 0], [255, 255, 255, 255]]),
                names=["Background", "Brain"],
            )

    filenames = compute_average_pet_signal_based_on_annotations(
        (tmp_path / "left_projection.mgh", tmp_path / "right_projection.mgh"),
        {
            "destrieux": FreeSurferAnnotationImage.from_raw(
                tmp_path / "destrieux_lh.annot",
                tmp_path / "destrieux_rh.annot",
            ),
            "desikan": FreeSurferAnnotationImage.from_raw(
                tmp_path / "desikan_lh.annot",
                tmp_path / "desikan_rh.annot",
            ),
        },
        tmp_path / "output",
    )

    assert set(filenames) == {
        tmp_path / "output" / "destrieux.tsv",
        tmp_path / "output" / "desikan.tsv",
    }
    for filename in filenames:
        assert pd.read_csv(filename, sep="\t").to_dict() == {
            "index": {0: 0, 1: 1, 2: 2, 3: 3},
            "label_name": {
                0: "Background_lh",
                1: "Background_rh",
                2: "Brain_lh",
                3: "Brain_rh",
            },
            "mean_scalar": {0: 0.0, 1: 0.0, 2: 2.0, 3: 1.0},
        }


@pytest.mark.parametrize(
    "is_longitudinal,expected",
    [
        (
            True,
            (
                "(.*(sub-.*)\\/(ses-.*)\\/pet\\/(long-.*)\\/surface_longitudinal)\\/midsurface\\/.*_hemi_([a-z]+)(.*)$",
                "\\1/\\2_\\3_\\4_hemi-\\5_midcorticalsurface",
            ),
        ),
        (
            False,
            (
                "(.*(sub-.*)\\/(ses-.*)\\/pet\\/surface)\\/midsurface\\/.*_hemi_([a-z]+)(.*)$",
                "\\1/\\2_\\3_hemi-\\4_midcorticalsurface",
            ),
        ),
    ],
)
def test_get_mid_surface_substitutions(is_longitudinal: bool, expected: tuple):
    from clinica.pipelines.pet.surface.utils import _get_mid_surface_substitutions  # noqa

    assert _get_mid_surface_substitutions(is_longitudinal) == expected


@pytest.fixture
def expected_regexp_substitution_for_projection_in_native_space(
    tracer: Tracer,
    region: SUVRReferenceRegion,
    is_longitudinal: bool,
):
    if is_longitudinal:
        return (
            "(.*(sub-.*)\\/(ses-.*)\\/pet\\/(long-.*)\\/surface_longitudinal)\\/projection_native\\/.*_hemi_([a-z]+).*",
            f"\\1/\\2_\\3_\\4_trc-{tracer.value}_pet_space-native_suvr-{region.value}_pvc-iy_hemi-\\5_projection.mgh",
        )
    return (
        "(.*(sub-.*)\\/(ses-.*)\\/pet\\/surface)\\/projection_native\\/.*_hemi_([a-z]+).*",
        f"\\1/\\2_\\3_trc-{tracer.value}_pet_space-native_suvr-{region.value}_pvc-iy_hemi-\\4_projection.mgh",
    )


@pytest.mark.parametrize("tracer", Tracer)
@pytest.mark.parametrize("region", SUVRReferenceRegion)
@pytest.mark.parametrize("is_longitudinal", (True, False))
def test_get_projection_in_native_space_substitutions(
    tracer,
    region,
    is_longitudinal,
    expected_regexp_substitution_for_projection_in_native_space,
):
    from clinica.pipelines.pet.surface.utils import (
        _get_projection_in_native_space_substitutions,
    )

    assert (
        _get_projection_in_native_space_substitutions(tracer, region, is_longitudinal)
        == expected_regexp_substitution_for_projection_in_native_space
    )


@pytest.fixture
def expected_regexp_substitution_for_projection_in_fsaverage(
    tracer: Tracer,
    region: SUVRReferenceRegion,
    is_longitudinal: bool,
):
    if is_longitudinal:
        return (
            "(.*(sub-.*)\\/(ses-.*)\\/pet\\/(long-.*)\\/surface_longitudinal)\\/projection_fsaverage\\/.*_hemi_([a-z]+).*_fwhm_([0-9]+).*",
            f"\\1/\\2_\\3_\\4_trc-{tracer.value}_pet_space-fsaverage_suvr-{region.value}_pvc-iy_hemi-\\5_fwhm-\\6_projection.mgh",
        )
    return (
        "(.*(sub-.*)\\/(ses-.*)\\/pet\\/surface)\\/projection_fsaverage\\/.*_hemi_([a-z]+).*_fwhm_([0-9]+).*",
        f"\\1/\\2_\\3_trc-{tracer.value}_pet_space-fsaverage_suvr-{region.value}_pvc-iy_hemi-\\4_fwhm-\\5_projection.mgh",
    )


@pytest.mark.parametrize("tracer", Tracer)
@pytest.mark.parametrize("region", SUVRReferenceRegion)
@pytest.mark.parametrize("is_longitudinal", (True, False))
def test_get_projection_in_fsaverage_substitution(
    tracer,
    region,
    is_longitudinal,
    expected_regexp_substitution_for_projection_in_fsaverage,
):
    from clinica.pipelines.pet.surface.utils import (
        _get_projection_in_fsaverage_substitution,
    )

    assert (
        _get_projection_in_fsaverage_substitution(tracer, region, is_longitudinal)
        == expected_regexp_substitution_for_projection_in_fsaverage
    )


@pytest.fixture
def expected_regexp_substitution_for_tsv_file_for_atlas(
    tracer: Tracer,
    region: SUVRReferenceRegion,
    atlas: str,
    is_longitudinal: bool,
):
    if is_longitudinal:
        return (
            f"(.*(sub-.*)\\/(ses-.*)\\/pet\\/(long-.*)\\/surface_longitudinal)\\/{atlas}_tsv\\/{atlas}.tsv",
            f"\\1/atlas_statistics/\\2_\\3_\\4_trc-{tracer.value}_pet_space-{atlas}_pvc-iy_suvr-{region.value}_statistics.tsv",
        )
    return (
        f"(.*(sub-.*)\\/(ses-.*)\\/pet\\/surface)\\/{atlas}_tsv\\/{atlas}.tsv",
        f"\\1/atlas_statistics/\\2_\\3_trc-{tracer.value}_pet_space-{atlas}_pvc-iy_suvr-{region.value}_statistics.tsv",
    )


@pytest.mark.parametrize("tracer", Tracer)
@pytest.mark.parametrize("region", SUVRReferenceRegion)
@pytest.mark.parametrize("atlas", ("destrieux", "desikan"))
@pytest.mark.parametrize("is_longitudinal", (True, False))
def test_get_tsv_file_for_atlas(
    tracer,
    region,
    atlas,
    is_longitudinal,
    expected_regexp_substitution_for_tsv_file_for_atlas,
):
    from clinica.pipelines.pet.surface.utils import _get_tsv_file_for_atlas  # noqa

    assert (
        _get_tsv_file_for_atlas(tracer, region, atlas, is_longitudinal)
        == expected_regexp_substitution_for_tsv_file_for_atlas
    )
