import pytest

from clinica.converters.study_models import StudyName


@pytest.mark.parametrize(
    "study,study_id,expected",
    [
        (StudyName.ADNI, "001_S_0001", "sub-ADNI001S0001"),
        (StudyName.NIFD, "1_S_0001", "sub-NIFD1S0001"),
        (StudyName.AIBL, "10", "sub-AIBL10"),
        (StudyName.UKB, "0101001", "sub-UKB0101001"),
        (StudyName.GENFI, "MAPT009", "sub-MAPT009"),
        (StudyName.OASIS3, "OAS30001", "sub-OAS30001"),
        (StudyName.HABS, "P_INIBUB", "sub-HABSINIBUB"),
        (StudyName.OASIS, "OAS1_0001_MR1", "sub-OASIS10001"),
        (StudyName.IXI, "IXI001", "sub-IXI001"),
    ],
)
def test_study_to_bids_id_passing(study: StudyName, study_id: str, expected: str):
    from clinica.converters.study_models import bids_id_factory

    assert bids_id_factory(study).from_original_study_id(study_id) == expected


@pytest.mark.parametrize(
    "study,study_id",
    [
        (StudyName.ADNI, "001S0001"),
        (StudyName.ADNI, "001_X_0001"),
        (StudyName.ADNI, "foo_S_0001"),
        (StudyName.NIFD, "1S0001"),
        (StudyName.NIFD, "1_X_0001"),
        (StudyName.NIFD, "foo_S_0001"),
        (StudyName.AIBL, "10A"),
        (StudyName.UKB, "0101001A"),
        (StudyName.GENFI, "MAPT009?"),
        (StudyName.OASIS3, "OAS3_0001"),
        (StudyName.OASIS3, "OAS3001"),
        (StudyName.HABS, "PINIBUB"),
        (StudyName.HABS, "X_INIBUB"),
        (StudyName.HABS, "P_INIBUB?"),
        (StudyName.OASIS, "OAS10001MR1"),
        (StudyName.OASIS, "OAS1_0001_MRI1"),
        (StudyName.OASIS, "OAS1_001_MR1"),
        (StudyName.IXI, "IXI_001"),
        (StudyName.IXI, "IXI0001"),
    ],
)
def test_study_to_bids_id_value_error(study: StudyName, study_id: str):
    from clinica.converters.study_models import bids_id_factory

    with pytest.raises(ValueError):
        bids_id_factory(study).from_original_study_id(study_id)


@pytest.mark.parametrize(
    "study,source_id,bids_id",
    [
        (StudyName.ADNI, "001_S_0001", "sub-ADNI001S0001"),
        (StudyName.NIFD, "1_S_0001", "sub-NIFD1S0001"),
        (StudyName.AIBL, "10", "sub-AIBL10"),
        (StudyName.UKB, "0101001", "sub-UKB0101001"),
        (StudyName.GENFI, "MAPT009", "sub-MAPT009"),
        (StudyName.OASIS3, "OAS30001", "sub-OAS30001"),
        (StudyName.HABS, "P_INIBUB", "sub-HABSINIBUB"),
        (StudyName.OASIS, "OAS1_0001_MR1", "sub-OASIS10001"),
        (StudyName.IXI, "IXI001", "sub-IXI001"),
    ],
)
def test_bids_to_study(study, bids_id, source_id):
    from clinica.converters.study_models import bids_id_factory

    assert bids_id_factory(study)(bids_id).to_original_study_id() == source_id
