from clinica.utils.testing_utils import build_bids_directory


def test_pet_volume_info_loading(tmp_path):
    from clinica.pipelines.pet.volume.pipeline import PETVolume

    bids = build_bids_directory(tmp_path / "bids", {"sub-01": ["ses-M00"]})
    pipeline = PETVolume(bids_directory=str(bids))

    assert pipeline.info == {
        "id": "aramislab/pet-volume",
        "author": "Jorge Samper-Gonzalez",
        "version": "0.1.0",
        "space_caps": "55M",
        "space_wd": "245M",
        "dependencies": [{"type": "software", "name": "spm", "version": ">=12"}],
    }


def test_pet_volume_dependencies(tmp_path, mocker):
    from packaging.specifiers import SpecifierSet
    from packaging.version import Version

    from clinica.pipelines.pet.volume.pipeline import PETVolume
    from clinica.utils.check_dependency import SoftwareDependency, ThirdPartySoftware

    mocker.patch(
        "clinica.utils.check_dependency._get_spm_version",
        return_value=Version("12.7219"),
    )
    bids = build_bids_directory(tmp_path / "bids", {"sub-01": ["ses-M00"]})
    pipeline = PETVolume(bids_directory=str(bids))

    assert pipeline.dependencies == [
        SoftwareDependency(
            ThirdPartySoftware.SPM, SpecifierSet(">=12"), Version("12.7219")
        )
    ]
