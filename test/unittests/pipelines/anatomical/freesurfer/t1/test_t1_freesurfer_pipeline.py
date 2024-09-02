from clinica.utils.testing_utils import build_bids_directory


def test_t1_freesurfer_info_loading(tmp_path):
    from clinica.pipelines.anatomical.freesurfer.t1.pipeline import T1FreeSurfer

    bids = build_bids_directory(tmp_path / "bids", {"sub-01": ["ses-M00"]})
    pipeline = T1FreeSurfer(bids_directory=str(bids))

    assert pipeline.info == {
        "id": "aramislab/t1-freesurfer",
        "author": ["Junhao Wen", "Alexandre Routier"],
        "version": "0.1",
        "space_caps": "1G",
        "space_wd": "1G",
        "dependencies": [
            {"type": "software", "name": "freesurfer", "version": ">=6.0.0"}
        ],
    }


def test_t1_freesurfer_dependencies(tmp_path, mocker):
    from packaging.specifiers import SpecifierSet
    from packaging.version import Version

    from clinica.pipelines.anatomical.freesurfer.t1.pipeline import T1FreeSurfer
    from clinica.utils.check_dependency import SoftwareDependency, ThirdPartySoftware

    mocker.patch(
        "clinica.utils.check_dependency._get_freesurfer_version",
        return_value=Version("7.2.1"),
    )
    bids = build_bids_directory(tmp_path / "bids", {"sub-01": ["ses-M00"]})
    pipeline = T1FreeSurfer(bids_directory=str(bids))

    assert pipeline.dependencies == [
        SoftwareDependency(
            ThirdPartySoftware.FREESURFER, SpecifierSet(">=6.0.0"), Version("7.2.1")
        )
    ]
