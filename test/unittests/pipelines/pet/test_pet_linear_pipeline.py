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
    from packaging.version import Version

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
