from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from clinica.pipelines.anatomical.freesurfer.utils import (
    ColumnType,  # noqa
    ImageID,  # noqa
    InfoType,  # noqa
)
from clinica.utils.image import HemiSphere


@pytest.mark.parametrize(
    "sub_id,expected",
    [
        ("sub-CLNC01_ses-M00", "sub-CLNC01_ses-M00"),
        ("sub-CLNC01_long-M00M18", "sub-CLNC01_long-M00M18"),
        (
            "sub-CLNC01_ses-M00.long.sub-CLNC01_long-M00M18",
            "sub-CLNC01_ses-M00_long-M00M18",
        ),
    ],
)
def test_get_prefix(sub_id, expected):
    from clinica.pipelines.anatomical.freesurfer.utils import _get_prefix  # noqa

    assert _get_prefix(sub_id) == expected


@pytest.mark.parametrize("column_type", ColumnType)
def test_read_stats_file_file_errors(tmp_path, column_type):
    from clinica.pipelines.anatomical.freesurfer.utils import _read_stats_file  # noqa

    with pytest.raises(
        FileNotFoundError,
        match=f"Stats file {tmp_path / 'foo.stats'} could not be found.",
    ):
        _read_stats_file(tmp_path / "foo.stats", column_type)  # noqa


@pytest.fixture
def fake_stats_dataframe(columns: list[str], n_rows: int) -> pd.DataFrame:
    """Generate a fake dummy stats Pandas dataframe for testing purposes."""
    return pd.DataFrame({k: ["foo"] * n_rows for k in columns})


@pytest.mark.parametrize(
    "columns,column_type,n_rows",
    [
        (
            [
                "StructName",
                "NumVert",
                "SurfArea",
                "GrayVol",
                "ThickAvg",
                "ThickStd",
                "MeanCurv",
                "GausCurv",
                "FoldInd",
                "CurvInd",
            ],
            ColumnType.PARCELLATION,
            6,
        ),
        (
            [
                "Index",
                "SegId",
                "NVoxels",
                "Volume_mm3",
                "StructName",
                "normMean",
                "normStdDev",
                "normMin",
                "normMax",
                "normRange",
            ],
            ColumnType.SEGMENTATION,
            6,
        ),
    ],
)
def test_read_stats_file(tmp_path, columns, column_type, n_rows, fake_stats_dataframe):
    from clinica.pipelines.anatomical.freesurfer.utils import _read_stats_file  # noqa

    fake_stats_dataframe.to_csv(
        tmp_path / "stats_file.csv", sep=" ", index=False, header=False
    )
    df = _read_stats_file(tmp_path / "stats_file.csv", column_type)
    assert len(df) == n_rows
    assert set(df.columns) == set(columns)


@pytest.mark.parametrize(
    "atlas,expected",
    [
        ("desikan", "aparc.stats"),
        ("destrieux", "aparc.a2009s.stats"),
        ("ba", "BA_exvivo.stats"),
        ("foo", "foo.stats"),
    ],
)
def test_get_stats_filename_for_atlas(tmp_path, atlas, expected):
    from clinica.pipelines.anatomical.freesurfer.utils import (
        _get_stats_filename_for_atlas,  # noqa
    )

    df = _get_stats_filename_for_atlas(tmp_path, atlas)
    for hemi in HemiSphere:
        assert df[hemi.value] == tmp_path / f"{hemi.value}.{expected}"


def test_generate_tsv_for_parcellation_errors(tmp_path):
    from clinica.pipelines.anatomical.freesurfer.utils import (
        _generate_tsv_for_parcellation,  # noqa
    )

    not_supported_atlas = "foo"
    with pytest.raises(
        FileNotFoundError,
        match=f"Stats file {tmp_path / 'lh.foo.stats'} could not be found.",
    ):
        _generate_tsv_for_parcellation(
            tmp_path,
            tmp_path,
            "sub-CLNC01_ses-M00",
            [not_supported_atlas],
        )


def generate_fake_secondary_stats_line(label: str, value: float, unit: str) -> str:
    """Generate a fake line of a freesurfer stats file."""
    unit_to_keyword = {"mm^3": "volume", "mm^2": "area", "mm": "thickness"}

    return f"# Measure {label}, foo, {unit_to_keyword[unit]}, {value}, {unit}"


def generate_fake_primary_stats(
    columns: list, regions: list, value: float, id_region: int, segmentation: bool
) -> str:
    """Generate fake primary stats lines."""
    content = f"# NTableCols {len(columns)}\n# "
    if segmentation:
        content = f"# NRows {len(regions)}\n# NTableCols {len(columns)}\n# ColHeaders "
    content += " ".join(columns)
    for i, region in enumerate(regions):
        row = [str(value) for _ in columns]
        row[id_region] = region
        row = " ".join(row)
        if segmentation:
            content += f"\n{i + 1} {row}"
        else:
            content += f"\n{row}"
    return content


COLUMNS_SEGMENTATION = [
    "Index",
    "SegId",
    "NVoxels",
    "Volume_mm3",
    "StructName",
    "normMean",
    "normStdDev",
    "normMin",
    "normMax",
    "normRange",
]
COLUMNS_PARCELLATION = [
    "StructName",
    "NumVert",
    "SurfArea",
    "GrayVol",
    "ThickAvg",
    "ThickStd",
    "MeanCurv",
    "GausCurv",
    "FoldInd",
    "CurvInd",
]
TEST_SEGMENTATION_STRUCT_NAMES = [
    "Left-Lateral-Ventricle",
    "Left-Inf-Lat-Vent",
    "Left-Cerebellum-White-Matter",
    "BrainSeg",
    "BrainSegNotVent",
    "BrainSegNotVentSurf",
]
TEST_PARCELLATION_STRUCT_NAMES = [
    "bankssts",
    "caudalanteriorcingulate",
    "caudalmiddlefrontal",
    "cuneus",
    "entorhinal",
    "fusiform",
    "inferiorparietal",
    "inferiortemporal",
    "isthmuscingulate",
    "lateraloccipital",
]


def generate_fake_segmentation_stats_file_content(dummy_number: float = 42.0) -> str:
    """Generate a portion of the content of a Freesurfer segmentation stats file.

    The portion contains the so-called secondary stats line (lines starting with
    '# Measure') as well as the table of primary stats. All other lines are ignored
    by the tested functions and not added here. All numbers will be set to the
    provided 'dummy_number'. This function should only be used for testing purposes.
    """
    secondary_stats_lines = "\n".join(
        [
            generate_fake_secondary_stats_line(label, dummy_number, unit)
            for label, unit in zip(
                ["region1", "region2", "region3", "region4"],
                ["mm^3", "mm^3", "mm^3", "mm^2"],
            )
        ]
    )
    primary_stats = generate_fake_primary_stats(
        COLUMNS_SEGMENTATION,
        TEST_SEGMENTATION_STRUCT_NAMES,
        dummy_number,
        id_region=4,
        segmentation=True,
    )
    return "\n".join([secondary_stats_lines, primary_stats])


def generate_fake_parcellation_stats_file_content(dummy_number: float = 42.0) -> str:
    """Generate a portion of the content of a Freesurfer parcellation stats file.

    The portion contains the so-called secondary stats line (lines starting with
    '# Measure') as well as the table of primary stats. All other lines are ignored
    by the tested functions and not added here. All numbers will be set to the
    provided 'dummy_number'. This function should only be used for testing purposes.
    """
    secondary_stats_lines = "\n".join(
        [
            generate_fake_secondary_stats_line(label, dummy_number, unit)
            for label, unit in zip(
                [
                    "volumetric_region1",
                    "volumetric_region2",
                    "area_region1",
                    "thickness_region1",
                    "volumetric_region3",
                    "area_region2",
                ],
                ["mm^3", "mm^3", "mm^2", "mm", "mm^3", "mm^2"],
            )
        ]
    )
    primary_stats = generate_fake_primary_stats(
        COLUMNS_PARCELLATION,
        TEST_PARCELLATION_STRUCT_NAMES,
        dummy_number,
        id_region=0,
        segmentation=False,
    )
    return "\n".join([secondary_stats_lines, primary_stats])


@pytest.mark.parametrize(
    "prefix",
    ["sub-CLNC01_ses-M00", "sub-CLNC01_long-M00M18", "sub-CLNC01_ses-M00_long-M00M18"],
)
def test_generate_tsv_for_segmentation(tmp_path, prefix):
    from clinica.pipelines.anatomical.freesurfer.utils import (
        _generate_tsv_for_segmentation,  # noqa
    )

    stats_folder = tmp_path / "stats"
    stats_folder.mkdir()
    for filename in ("aseg", "wmparc"):
        with open(stats_folder / f"{filename}.stats", "w") as fp:
            fp.write(generate_fake_segmentation_stats_file_content(dummy_number=42.0))
    _generate_tsv_for_segmentation(stats_folder, tmp_path, prefix)
    for filename in ("segmentationVolumes.tsv", "parcellation-wm_volume.tsv"):
        tsv_file = tmp_path / f"{prefix}_{filename}"
        assert tsv_file.exists()
        df = pd.read_csv(tsv_file, sep="\t")
        assert np.all(df["label_value"].values == 42.0)
        assert set(df["label_name"].values) == set(
            [f"region{i}" for i in range(1, 4)] + TEST_SEGMENTATION_STRUCT_NAMES
        )


@pytest.mark.parametrize(
    "prefix",
    ["sub-CLNC01_ses-M00", "sub-CLNC01_long-M00M18", "sub-CLNC01_ses-M00_long-M00M18"],
)
def test_generate_tsv_for_parcellation(tmp_path, prefix):
    from clinica.pipelines.anatomical.freesurfer.utils import (
        _generate_tsv_for_parcellation,  # noqa
    )

    atlases = ["desikan", "destrieux", "ba"]
    stats_folder = tmp_path / "stats"
    stats_folder.mkdir()
    for filename in ("aparc", "aparc.a2009s", "BA_exvivo"):
        for hemi in HemiSphere:
            with open(stats_folder / f"{hemi.value}.{filename}.stats", "w") as fp:
                fp.write(
                    generate_fake_parcellation_stats_file_content(dummy_number=42.0)
                )
    _generate_tsv_for_parcellation(stats_folder, tmp_path, prefix, atlases)
    for atlas in atlases:
        for info_type in (InfoType.VOLUME, InfoType.AREA, InfoType.THICKNESS):
            tsv_file = tmp_path / f"{prefix}_parcellation-{atlas}_{info_type.value}.tsv"
            assert tsv_file.exists()
            df = pd.read_csv(tsv_file, sep="\t")
            assert np.all(df["label_value"].values == 42.0)
            if info_type == InfoType.VOLUME:
                expected_secondary_labels = [
                    f"volumetric_region{i}" for i in range(1, 4)
                ]
            elif info_type == InfoType.AREA:
                expected_secondary_labels = [f"area_region{i}" for i in range(1, 3)]
            elif info_type == InfoType.THICKNESS:
                expected_secondary_labels = [
                    f"thickness_region{i}" for i in range(1, 2)
                ]
            assert set(df["label_name"].values) == set(
                expected_secondary_labels + TEST_PARCELLATION_STRUCT_NAMES
            )


def test_generate_regional_measures_file_not_found_error(tmp_path):
    from clinica.pipelines.anatomical.freesurfer.utils import generate_regional_measures

    with pytest.raises(
        FileNotFoundError,
        match="Image sub-CLNC01 | ses-M00 does not contain FreeSurfer segmentation",
    ):
        generate_regional_measures(tmp_path, "sub-CLNC01_ses-M00", ["destrieux"])


def test_generate_regional_measures_type_error(tmp_path):
    from clinica.pipelines.anatomical.freesurfer.utils import generate_regional_measures

    with pytest.raises(
        ValueError,
        match="atlases should be a list of strings. <class 'str'> was provided instead.",
    ):
        generate_regional_measures(tmp_path, "sub-CLNC01_ses-M00", "destrieux")  # noqa


@pytest.mark.parametrize("sub_id", ["", "fooo"])
def test_extract_image_id_from_freesurfer_id_error(sub_id):
    from clinica.pipelines.anatomical.freesurfer.utils import (
        extract_image_id_from_freesurfer_id,
    )

    with pytest.raises(
        ValueError,
        match=f"The provided Freesurfer ID {sub_id} could not be parsed.",
    ):
        extract_image_id_from_freesurfer_id(sub_id)


@pytest.mark.parametrize(
    "sub_id,expected",
    [
        ("sub-CLNC01_ses-M00", ImageID("sub-CLNC01", "ses-M00", "")),
        ("sub-CLNC01_long-M00M18", ImageID("sub-CLNC01", "", "long-M00M18")),
        (
            "sub-CLNC01_ses-M00.long.sub-CLNC01_long-M00M18",
            ImageID("sub-CLNC01", "ses-M00", "long-M00M18"),
        ),
    ],
)
def test_extract_image_id_from_longitudinal_segmentation(sub_id, expected):
    from clinica.pipelines.anatomical.freesurfer.utils import (
        extract_image_id_from_freesurfer_id,
    )

    assert extract_image_id_from_freesurfer_id(sub_id) == expected


@pytest.mark.parametrize(
    "line,info,expected",
    [
        ("foo", InfoType.VOLUME, None),
        ("", InfoType.VOLUME, None),
        ("# Measure, foo", InfoType.THICKNESS, None),
        (
            "# Measure BrainSeg, BrainSegVol, Brain Segmentation Volume, 42.00, mm^3",
            InfoType.VOLUME,
            ("BrainSeg", "42.00"),
        ),
        (
            "# Measure BrainSeg, BrainSegVol, Brain Segmentation Volume, 42.00, mm^3",
            InfoType.AREA,
            None,
        ),
        (
            "# Measure Cortex, MeanThickness, Mean Thickness, 2.6, mm",
            InfoType.THICKNESS,
            ("Cortex", "2.6"),
        ),
        (
            "# Measure Cortex, MeanThickness, Mean Thickness, 2.6, mm",
            InfoType.VOLUME,
            None,
        ),
        (
            "# Measure Cortex, WhiteSurfArea, White Surface Total Area, 82.1, mm^2",
            InfoType.AREA,
            ("Cortex", "82.1"),
        ),
        (
            "# Measure Cortex, WhiteSurfArea, White Surface Total Area, 82.1, mm^2",
            InfoType.THICKNESS,
            None,
        ),
    ],
)
def test_extract_region_and_stat_value(line, info, expected):
    from clinica.pipelines.anatomical.freesurfer.utils import (
        _extract_region_and_stat_value,  # noqa
    )

    assert _extract_region_and_stat_value(line, info) == expected


@pytest.mark.parametrize(
    "line,info,expected",
    [
        ("", InfoType.AREA, False),
        ("foo", InfoType.VOLUME, False),
        ("#Measure", InfoType.VOLUME, False),
        ("# Measure", InfoType.THICKNESS, False),
        ("# Measure foo mm^3", InfoType.VOLUME, True),
        ("# Measure foo mm^2", InfoType.AREA, True),
        ("# Measure foo mm", InfoType.THICKNESS, True),
        ("# Measure foo, bar baz, mm^3", InfoType.VOLUME, True),
        ("#Measure foo mm^3", InfoType.VOLUME, False),
    ],
)
def test_filter_stats_line(line, info, expected):
    from clinica.pipelines.anatomical.freesurfer.utils import _filter_stats_line  # noqa

    assert _filter_stats_line(line, info) == expected


@pytest.fixture
def stats_file(tmp_path) -> Path:
    stats_file_content = (
        "# Measure BrainSeg, BrainSegVol, Brain Segmentation Volume, 42.00, mm^3\n"
        "# foo\n#bar\n"
        "# Measure Cortex, MeanThickness, Mean Thickness, 2.6, mm\n"
        "# Measure BrainSegNotVentSurf, BrainSegVolNotVentSurf, Brain Segmentation "
        "Volume Without Ventricles from Surf, 111.726, mm^3\n"
        "# Measure Cortex, WhiteSurfArea, White Surface Total Area, 82.1, mm^2\n"
        "# foo\n#bar\n"
        "foo bar baz"
    )
    with open(tmp_path / "stats_file.stats", "w") as fp:
        fp.write(stats_file_content)

    return tmp_path / "stats_file.stats"


def test_get_secondary_stats_mean_curvature(tmp_path, stats_file):
    from clinica.pipelines.anatomical.freesurfer.utils import _get_secondary_stats  # noqa

    assert _get_secondary_stats(stats_file, InfoType.MEANCURV) is None


def test_get_secondary_stats_volume(tmp_path, stats_file):
    from clinica.pipelines.anatomical.freesurfer.utils import _get_secondary_stats  # noqa

    df = _get_secondary_stats(stats_file, InfoType.VOLUME)

    assert list(df["label_name"]) == ["BrainSeg", "BrainSegNotVentSurf"]
    assert list(df["label_value"]) == ["42.00", "111.726"]


def test_get_secondary_stats_area(tmp_path, stats_file):
    from clinica.pipelines.anatomical.freesurfer.utils import _get_secondary_stats  # noqa

    df = _get_secondary_stats(stats_file, InfoType.AREA)

    assert list(df["label_name"]) == ["Cortex"]
    assert list(df["label_value"]) == ["82.1"]


def test_get_secondary_stats_thickness(tmp_path, stats_file):
    from clinica.pipelines.anatomical.freesurfer.utils import _get_secondary_stats  # noqa

    df = _get_secondary_stats(stats_file, InfoType.THICKNESS)

    assert list(df["label_name"]) == ["Cortex"]
    assert list(df["label_value"]) == ["2.6"]
