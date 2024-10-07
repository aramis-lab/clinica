from pathlib import Path

import pandas as pd
import pytest

from clinica.iotools.converters.oasis_to_bids.oasis_to_bids_utils import (
    create_sessions_dict,
    write_sessions_tsv,
)


def build_clinical_data(tmp_path: Path) -> None:
    df = pd.DataFrame(
        {
            "ID": ["OAS1_0001_MR1", "OAS1_0002_MR1"],
            "M/F": ["F", "M"],
            "Hand": ["R", "L"],
            "Age": [74, 67],
            "Educ": [2, 2],
            "SES": [3, 3],
            "MMSE": [29, 29],
            "CDR": [0, 0],
            "eTIV": [1344, 1344],
            "nWBV": [0.704, 0.645],
            "ASF": [1.306, 1.100],
            "Delay": [float("nan"), float("nan")],
        }
    )
    df.to_csv(tmp_path / "oasis_cross-sectional.csv", index=False)

    # todo : do excel for future
    bids_ids = ["sub-OASIS10001", "sub-OASIS10002"]


def test_create_sessions_dict():
    # todo
    clinical_specifications_folder = Path("clinica/iotools/converters/specifications")
    pass


def test_write_sessions_tsv():
    # todo
    pass
