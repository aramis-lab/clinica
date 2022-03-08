import os
import pytest
import pandas as pd
from pathlib import Path
import xml.etree.ElementTree as ET


def _get_xml_templates():
    suffix = "_template.xml"
    current_dir = os.path.dirname(os.path.realpath(__file__))
    return [
        _.rstrip(suffix)
        for _ in os.listdir(Path(current_dir) / "data")
        if _.endswith(suffix)
    ]


@pytest.fixture
def basic_xml_tree():
    """Basic XML tree used for testing."""
    xml_tree = ET.Element("root")
    xml_tree.text = "100"
    xml_level_1 = [ET.SubElement(xml_tree, "leaf1"), ET.SubElement(xml_tree, "leaf2")]
    return xml_tree


def test_check_xml_tag():
    """Test function `_check_xml_tag`."""
    from clinica.iotools.converters.adni_to_bids.adni_json import _check_xml_tag
    _check_xml_tag("foo", "foo")
    with pytest.raises(ValueError, match="Bad tag"):
        _check_xml_tag("foo", "bar")


def test_check_xml_nb_children(basic_xml_tree):
    """Test function `_check_xml_nb_children`."""
    from clinica.iotools.converters.adni_to_bids.adni_json import _check_xml_nb_children
    _check_xml_nb_children(basic_xml_tree, 2)
    _check_xml_nb_children(basic_xml_tree, [1, 2, 6])
    with pytest.raises(ValueError,
                       match="Bad number of children for <root>: got 2 != 3"):
        _check_xml_nb_children(basic_xml_tree, 3)
    with pytest.raises(ValueError,
                       match="Bad number of children for <root>: got 2, not in"):
        _check_xml_nb_children(basic_xml_tree, [1, 3])


def test_check_xml(basic_xml_tree):
    """Test function `_check_xml`."""
    from clinica.iotools.converters.adni_to_bids.adni_json import _check_xml
    assert _check_xml(basic_xml_tree, "root", 2) == basic_xml_tree
    assert _check_xml(basic_xml_tree, "root", [1, 2]) == basic_xml_tree
    assert _check_xml(basic_xml_tree, "root") == basic_xml_tree
    with pytest.raises(ValueError, match="Bad tag"):
        _check_xml(basic_xml_tree, "foo")


def test_get_text(basic_xml_tree):
    """Test function `_get_text`."""
    from clinica.iotools.converters.adni_to_bids.adni_json import _get_text
    assert _get_text(basic_xml_tree) == "100"
    assert _get_text(basic_xml_tree, cast=int) == 100


def test_check_xml_and_get_text(basic_xml_tree):
    """Test function `_check_xml_and_get_text`."""
    from clinica.iotools.converters.adni_to_bids.adni_json import _check_xml_and_get_text
    xml_leaf = ET.Element("leaf")
    xml_leaf.text = "12"
    assert _check_xml_and_get_text(xml_leaf,  "leaf") == "12"
    assert _check_xml_and_get_text(xml_leaf,  "leaf", cast=int) == 12
    with pytest.raises(ValueError, match="Bad number of children for <root>"):
        _check_xml_and_get_text(basic_xml_tree,  "root") 


def test_read_xml_files(tmp_path):
    """Test function `_read_xml_files`."""
    from clinica.iotools.converters.adni_to_bids.adni_json import _read_xml_files
    with pytest.raises(IndexError, match="No ADNI xml files"):
        _read_xml_files()
    xml_path = tmp_path / "xml_files"
    xml_path.mkdir()
    os.chdir(xml_path)
    clinica_path = xml_path / "Clinica_processed_metadata"
    clinica_path.mkdir()
    dummy_file = clinica_path / "ADNI_1234.xml"
    with dummy_file.open("w") as fp:
        fp.write("foo")
    assert os.listdir() == ['Clinica_processed_metadata']
    assert _read_xml_files() == ['Clinica_processed_metadata/ADNI_1234.xml']
    subjects = ["01", "02", "06", "12"]
    for subj in subjects:
        with open(xml_path / f"ADNI_{subj}.xml", "w") as fp:
            fp.write(f"foo {subj}")
    assert _read_xml_files(subjects, xml_path) == [str(xml_path / f"ADNI_{subj}.xml") for subj in subjects]
    os.chdir(os.path.dirname(__file__))


def _load_xml_from_template(
        template_id,
        project="ADNI",
        modality="MRI",
        acq_time=pd.Timestamp(2017, 1, 1, 12),
):
    from string import Template
    other_substitutes = {
            "study_id": 100,
            "series_id": 200,
            "image_id": 300,
            "proc_id": 3615,
    }
    temp = Path(f"./data/{template_id}_template.xml").read_text()
    temp = Template(temp.replace("\n", ""))
    return temp.safe_substitute(
            project=project, modality=modality,
            acq_time=acq_time, **other_substitutes
    )


def _write_xml_example(base_path, template_id="ADNI_123_S_4567", **kwargs):
    suffix = kwargs.pop("suffix", None)
    xml_content = _load_xml_from_template(template_id, **kwargs)
    filename = f"{template_id}.xml"
    if suffix is not None:
        filename = f"{template_id}_{suffix}.xml"
    xml_file_path = base_path / filename
    with open(xml_file_path, "w") as fp:
        fp.write(xml_content)
    return xml_file_path


@pytest.mark.parametrize("template_id", _get_xml_templates())
def test_get_root_from_xml_path(tmp_path, template_id):
    """Test function `_get_root_from_xml_path`."""
    from clinica.iotools.converters.adni_to_bids.adni_json import _get_root_from_xml_path
    xml_file = _write_xml_example(tmp_path, template_id=template_id)
    assert isinstance(_get_root_from_xml_path(xml_file), ET.Element)


@pytest.fixture
def expected_image_metadata(template_id):
    expected = {}
    if template_id == "ADNI_234_S_5678":
        expected = {
                'image_proc_pipe': 'Grinder Pipeline',
                'image_proc_id': 3615,
                'image_proc_desc': 'MT1; GradWarp; N3m'
        }
    elif template_id == "ADNI_345_S_6789":
        expected = {
                'image_proc_pipe': 'UCSD ADNI Pipeline',
                'image_proc_id': 3615,
                'image_proc_desc': 'MPR; GradWarp; N3; Scaled'
        }
    return expected


@pytest.mark.parametrize("template_id", _get_xml_templates())
def test_parsing(tmp_path, template_id, expected_image_metadata):
    """Test function `_get_root_from_xml_path`."""
    from clinica.iotools.converters.adni_to_bids.adni_json import (
            _get_root_from_xml_path, _parse_project, _parse_subject,
            _parse_study, _parse_series, _parse_images,
    )
    expected_subject_id = template_id[5:]
    xml_files = {
            "correct": _write_xml_example(
                tmp_path, template_id=template_id, suffix="correct"),
            "bad_project": _write_xml_example(
                tmp_path, project="foo",
                template_id=template_id, suffix="bad_project"),
            "bad_study": _write_xml_example(
                tmp_path, modality="bar",
                template_id=template_id, suffix="bad_study"),
    }
    roots = {k: _get_root_from_xml_path(v) for k, v in xml_files.items()}

    # _parse_project
    bad = roots.pop("bad_project")
    with pytest.raises(ValueError, match="Not ADNI cohort"):
        _parse_project(bad)
    with pytest.raises(ValueError, match="XML root should have only one child."):
        _parse_project(ET.Element("root"))
    projects = {k: _parse_project(root) for k, root in roots.items()}
    assert all([isinstance(project, ET.Element) for project in projects.values()])

    # _parse_subject
    subjects = {k: _parse_subject(project) for k, project in projects.items()}
    assert all([_[0] == expected_subject_id for _ in subjects.values()])
    assert all([isinstance(_[1], ET.Element) for _ in subjects.values()])
    bad_project = ET.Element("project")
    bad_subject = ET.SubElement(bad_project, "subject")
    [ET.SubElement(bad_subject, f"{i}") for i in range(6)]
    with pytest.raises(ValueError,
                       match="Bad number of children for <subject>: got 6"):
        _parse_subject(bad_project)

    # _parse_study
    bad = subjects.pop("bad_study")
    with pytest.raises(ValueError, match="Unexpected modality bar"):
        _parse_study(bad[1])
    studies = {k: _parse_study(subject[1]) for k, subject in subjects.items()}
    assert all([isinstance(study, ET.Element) for study in studies.values()])

    # _parse_series
    series = {k: _parse_series(study) for k, study in studies.items()}
    assert all([isinstance(serie, ET.Element) for serie in series.values()])
    for length_series, length_study in zip([3, 4], [6, 7]):
        bad_study = ET.Element("study")
        bad_series = [ET.SubElement(bad_study, "series") for i in range(length_study)]
        [ET.SubElement(bad_series[5], f"{i}") for i in range(length_series)]
        with pytest.raises(ValueError,
                           match=f"Bad number of children for <series>: got {length_series}"):
            _parse_series(bad_study)

    # _parse_images
    images = {k: _parse_images(study) for k, study in studies.items()}
    assert all([isinstance(_[0], ET.Element) for _ in images.values()])
    if template_id == "ADNI_345_S_6789":
        assert all([_[1] == 2 for _ in images.values()])
    else:
        assert all([_[1] is None for _ in images.values()])
    assert all([_[2] == expected_image_metadata for _ in images.values()])


@pytest.fixture
def expected_mprage(template_id):
    expected = {
            "ADNI_123_S_4567": "Accelerated Sagittal MPRAGE",
            "ADNI_234_S_5678": "MPRAGE GRAPPA2",
            "ADNI_345_S_6789": "MP-RAGE",
    }
    return expected[template_id]


@pytest.mark.parametrize("template_id", _get_xml_templates())
def test_parse_xml_file(template_id, tmp_path, expected_mprage):
    """Test function `_parse_xml_file`."""
    from clinica.iotools.converters.adni_to_bids.adni_json import _parse_xml_file
    xml_file = _write_xml_example(tmp_path, template_id=template_id)
    scan_metadata = _parse_xml_file(xml_file)
    expected_subject_id = template_id[5:]
    assert scan_metadata["id"] == expected_subject_id
    assert scan_metadata["acq_time"] == pd.Timestamp(2017, 1, 1, 12)
    assert scan_metadata["image_orig_id"] == 300
    assert scan_metadata["image_orig_seq"] == expected_mprage
    assert scan_metadata["MRAcquisitionType"] == "3D"
    expected_pulse = "GR/IR"
    expected_manufacturer = "SIEMENS"
    expected_strength = "3.0"
    if template_id == "ADNI_345_S_6789":
        expected_pulse = "RM"
        expected_manufacturer = "GE MEDICAL SYSTEMS"
        expected_strength = "1.5"
    assert scan_metadata["PulseSequenceType"] == expected_pulse
    assert scan_metadata["Manufacturer"] == expected_manufacturer
    assert scan_metadata["MagneticFieldStrength"] == expected_strength

