import pytest


@pytest.mark.parametrize(
    "dictionnary, expected",
    [
        (
            {
                "TotalReadoutTime": 1,
                "EstimatedTotalReadoutTime": "",
                "PhaseEncodingSteps": "",
                "PixelBandwidth": "",
                "PhaseEncodingDirection": "j+",
                "PhaseEncodingAxis": "",
            },
            [1, "j+"],
        ),
        (
            {
                "TotalReadoutTime": "",
                "EstimatedTotalReadoutTime": 1,
                "PhaseEncodingSteps": "",
                "PixelBandwidth": "",
                "PhaseEncodingDirection": "j",
                "PhaseEncodingAxis": "",
            },
            [1, "j"],
        ),
        (
            {
                "TotalReadoutTime": "",
                "EstimatedTotalReadoutTime": "",
                "PhaseEncodingSteps": 1,
                "PixelBandwidth": 0.5,
                "PhaseEncodingDirection": "j-",
                "PhaseEncodingAxis": "",
            },
            [2, "j-"],
        ),
    ],
)
def test_extract_metadata_from_dwi_json(dictionnary, expected):
    """This function tests something."""
    import json

    from clinica.utils.filemanip import extract_metadata_from_dwi_json

    # create a json with parametrize ?
    with open("metadata.json", "w") as outfile:
        json.dump(dictionnary, outfile)
    assert (
        extract_metadata_from_dwi_json(
            "metadata.json",
            [
                "TotalReadoutTime",
                "EstimatedTotalReadoutTime",
                "PhaseEncodingSteps",
                "PixelBandwidth",
                "PhaseEncodingDirection",
                "PhaseEncodingAxis",
            ],
        )
        == expected
    )
