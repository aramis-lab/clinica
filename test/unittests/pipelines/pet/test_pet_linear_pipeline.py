from packaging.version import Version

from clinica.utils.bids import Visit
from clinica.utils.testing_utils import build_bids_directory, build_caps_directory


def test_pet_linear_info_loading(tmp_path):
    from clinica.pipelines.pet.linear.pipeline import PETLinear

    bids = build_bids_directory(tmp_path / "bids", {"sub-01": ["ses-M00"]})
    pipeline = PETLinear(bids_directory=str(bids))

    assert pipeline.info == {
        "id": "aramislab/pet_linear",
        "author": "Ravi Hassanaly",
        "version": "0.1.0",
        "space_caps": "50M",
        "space_wd": "160M",
        "dependencies": [{"type": "software", "name": "ants", "version": ">=2.2.0"}],
    }


def test_pet_linear_dependencies(tmp_path, mocker):
    from packaging.specifiers import SpecifierSet

    from clinica.pipelines.pet.linear.pipeline import PETLinear
    from clinica.utils.check_dependency import SoftwareDependency, ThirdPartySoftware

    mocker.patch(
        "clinica.utils.check_dependency._get_ants_version",
        return_value=Version("2.3.1"),
    )
    bids = build_bids_directory(tmp_path / "bids", {"sub-01": ["ses-M00"]})
    pipeline = PETLinear(bids_directory=str(bids))

    assert pipeline.dependencies == [
        SoftwareDependency(
            ThirdPartySoftware.ANTS, SpecifierSet(">=2.2.0"), Version("2.3.1")
        )
    ]


def test_pet_linear_get_processed_visits_empty(tmp_path, mocker):
    from clinica.pipelines.pet.linear.pipeline import PETLinear
    from clinica.utils.pet import SUVRReferenceRegion, Tracer

    mocker.patch(
        "clinica.utils.check_dependency._get_ants_version",
        return_value=Version("2.3.1"),
    )
    bids = build_bids_directory(
        tmp_path / "bids", {"sub-01": ["ses-M000", "ses-M006"], "sub-02": ["ses-M000"]}
    )
    caps = build_caps_directory(tmp_path / "caps", {})

    pipeline = PETLinear(
        bids_directory=str(bids),
        caps_directory=str(caps),
        parameters={
            "acq_label": Tracer.FDG,
            "suvr_reference_region": SUVRReferenceRegion.PONS,
            "uncropped_image": False,
        },
    )
    assert pipeline.get_processed_visits() == []


def test_pet_linear_get_processed_visits(tmp_path, mocker):
    """Test the get_processed_visits for PETLinear.

    We build the following CAPS dataset:

    caps
    ├── dataset_description.json
    └── subjects
        ├── sub-01
        │         ├── ses-M000
        │         │         └── pet_linear
        │         │             ├── sub-01_ses-M000_trc-18FAV45_pet_space-MNI152NLin2009cSym_desc-Crop_res-1x1x1_suvr-pons2_pet.nii.gz
        │         │             └── sub-01_ses-M000_trc-18FAV45_pet_space-T1w_rigid.mat
        │         └── ses-M006
        │             └── pet_linear
        │                 ├── sub-01_ses-M006_trc-18FAV45_pet_space-MNI152NLin2009cSym_desc-Crop_res-1x1x1_suvr-pons2_pet.nii.gz
        │                 ├── sub-01_ses-M006_trc-18FAV45_pet_space-T1w_rigid.mat
        │                 ├── sub-01_ses-M006_trc-18FFDG_pet_space-MNI152NLin2009cSym_desc-Crop_res-1x1x1_suvr-pons_pet.nii.gz
        │                 └── sub-01_ses-M006_trc-18FFDG_pet_space-T1w_rigid.mat
        └── sub-02
            ├── ses-M000
            │         └── pet_linear
            │             ├── sub-02_ses-M000_trc-18FFDG_pet_space-MNI152NLin2009cSym_desc-Crop_res-1x1x1_suvr-pons_pet.nii.gz
            │             └── sub-02_ses-M000_trc-18FFDG_pet_space-T1w_rigid.mat
            └── ses-M006
                └── pet_linear
                    ├── sub-02_ses-M006_trc-18FAV45_pet_space-MNI152NLin2009cSym_desc-Crop_res-1x1x1_suvr-pons2_pet.nii.gz
                    └── sub-02_ses-M006_trc-18FAV45_pet_space-T1w_rigid.mat

    And make sure that a PETLinear pipeline with tracer 18FFDG et pons SUVR will only consider
    (sub-01, ses-M006) and (sub-02, ses-M000) as already processed. The other folders contain
    PET images for the pet-linear pipeline but they were obtained with different tracers and SUVR.
    """
    from clinica.pipelines.pet.linear.pipeline import PETLinear
    from clinica.utils.pet import SUVRReferenceRegion, Tracer

    mocker.patch(
        "clinica.utils.check_dependency._get_ants_version",
        return_value=Version("2.3.1"),
    )
    bids = build_bids_directory(
        tmp_path / "bids", {"sub-01": ["ses-M000", "ses-M006"], "sub-02": ["ses-M000"]}
    )
    caps = build_caps_directory(
        tmp_path / "caps",
        {
            "pipelines": {
                "pet_linear": {
                    "acq_label": Tracer.FDG,
                    "suvr_reference_region": SUVRReferenceRegion.PONS,
                    "uncropped_image": False,
                }
            },
            "subjects": {
                "sub-01": ["ses-M006"],
                "sub-02": ["ses-M000"],
            },
        },
    )
    caps = build_caps_directory(
        tmp_path / "caps",
        {
            "pipelines": {
                "pet_linear": {
                    "acq_label": Tracer.AV45,
                    "suvr_reference_region": SUVRReferenceRegion.PONS2,
                    "uncropped_image": False,
                }
            },
            "subjects": {
                "sub-01": ["ses-M000", "ses-M006"],
                "sub-02": ["ses-M006"],
            },
        },
    )

    pipeline = PETLinear(
        bids_directory=str(bids),
        caps_directory=str(caps),
        parameters={
            "acq_label": Tracer.FDG,
            "suvr_reference_region": SUVRReferenceRegion.PONS,
            "uncropped_image": False,
        },
    )

    assert pipeline.get_processed_visits() == [
        Visit("sub-01", "ses-M006"),
        Visit("sub-02", "ses-M000"),
    ]
