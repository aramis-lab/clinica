import pytest


@pytest.mark.parametrize(
    "filename",
    [
        "foo.txt",
        "dwi.nii.gz",
        "sub-01_ses-M000_dwi.nii.gz",
        "sub-01_ses-M000_preproc.nii.gz",
        "sub-01_ses-M000_space-T1w_preproc.nii.gz",
        "sub-01_ses-M000_space-b0_preproc.nii.gz",
    ],
)
def test_get_caps_filenames_error(tmp_path, filename):
    from clinica.pipelines.dwi_connectome.dwi_connectome_utils import get_caps_filenames

    with pytest.raises(ValueError, match="is not in a CAPS compliant format."):
        get_caps_filenames(str(tmp_path / filename))


def test_get_caps_filenames(tmp_path):
    from clinica.pipelines.dwi_connectome.dwi_connectome_utils import get_caps_filenames

    dwi_caps = tmp_path / "dwi" / "preprocessing"
    dwi_caps.mkdir(parents=True)

    assert get_caps_filenames(
        str(dwi_caps / "sub-01_ses-M000_space-b0_desc-preproc_dwi.nii.gz")
    ) == (
        "sub-01_ses-M000_space-b0_desc-preproc_model-CSD_responseFunction.txt",
        "sub-01_ses-M000_space-b0_desc-preproc_model-CSD_diffmodel.nii.gz",
        "sub-01_ses-M000_space-b0_desc-preproc_model-CSD_tractography.tck",
        [
            "sub-01_ses-M000_space-b0_desc-preproc_atlas-desikan_parcellation.nii.gz",
            "sub-01_ses-M000_space-b0_desc-preproc_atlas-destrieux_parcellation.nii.gz",
        ],
        [
            "sub-01_ses-M000_model-CSD_atlas-desikan_connectivity.tsv",
            "sub-01_ses-M000_model-CSD_atlas-destrieux_connectivity.tsv",
        ],
    )
