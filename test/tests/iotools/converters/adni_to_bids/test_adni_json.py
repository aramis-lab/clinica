import os
import pytest
import pandas as pd
from pathlib import Path
import xml.etree.ElementTree as ET


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
        template_id="ADNI_123_S_4567",
        project="ADNI",
        modality="MRI",
        acq_time=pd.Timestamp(2017, 1, 1, 12),
):
    """TODO."""
    from string import Template
    other_substitutes = {
            "study_id": 100,
            "series_id": 200,
            "image_id": 300,
    }
    temp = Path(f"./data/{template_id}_template.xml").read_text()
    temp = Template(temp.replace("\n", ""))
    return temp.safe_substitute(
            project=project, modality=modality,
            acq_time=acq_time, **other_substitutes
    )


def _write_xml_example(base_path, template_id="ADNI_123_S_4567", **kwargs):
    """TODO."""
    xml_content = _load_xml_from_template(**kwargs)
    xml_file_path = base_path / f"{template_id}.xml"
    with open(xml_file_path, "w") as fp:
        fp.write(xml_content)
    return xml_file_path


def test_get_root_from_xml_path(tmp_path):
    """Test function `_get_root_from_xml_path`."""
    from clinica.iotools.converters.adni_to_bids.adni_json import _get_root_from_xml_path
    xml_file = _write_xml_example(tmp_path)
    assert isinstance(_get_root_from_xml_path(xml_file), ET.Element)


def test_parsing(tmp_path):
    """Test function `_get_root_from_xml_path`."""
    from clinica.iotools.converters.adni_to_bids.adni_json import (
            _get_root_from_xml_path, _parse_project, _parse_subject,
            _parse_study, _parse_series,
    )
    template_id = "123_S_4567"
    xml_files = {
            "correct": _write_xml_example(
                tmp_path, template_id=f"ADNI_{template_id}_correct"),
            "bad_project": _write_xml_example(
                tmp_path, project="foo",
                template_id=f"ADNI_{template_id}_bad_project"),
            "bad_study": _write_xml_example(
                tmp_path, modality="bar",
                template_id=f"ADNI_{template_id}_bad_study"),
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
    assert all([_[0] == template_id for _ in subjects.values()])
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
