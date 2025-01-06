import xml.etree.ElementTree
from pathlib import Path
from typing import Iterable, Optional, Tuple, Union

import pandas as pd

from clinica.utils.exceptions import ClinicaXMLParserError
from clinica.utils.stream import cprint, log_and_raise

__all__ = ["create_json_metadata"]


LOGGING_HEADER = "[ADNI JSON]"


def _log_and_raise_xml_parser_error(message: str):
    log_and_raise(f"{LOGGING_HEADER} {message}", ClinicaXMLParserError)


# Map from names of extracted metadata to their proper BIDS names.
METADATA_NAME_MAPPING = {
    "acquisition_type": "MRAcquisitionType",
    "pulse_sequence": "PulseSequenceType",
    "manufacturer": "Manufacturer",
    "field_strength": "MagneticFieldStrength",
}


def _bids_id_to_loni(bids_id: str) -> Optional[str]:
    """Convert a subject id of the form sub-ADNI000S0000
    back to original format 000_S_0000
    """
    from clinica.iotools.bids_utils import StudyName, bids_id_factory

    try:
        return bids_id_factory(StudyName.ADNI)(bids_id).to_original_study_id()
    except ValueError:
        return None


def _read_xml_files(
    subject_ids: Optional[Iterable[str]] = None,
    xml_path: Optional[Path] = None,
) -> Iterable[Path]:
    """Return the XML files in the folder `xml_path` for the provided `subj_ids`.
    This function assumes that file are named "ADNI_{sub_ids}.xml".
    If no files were found, an `IndexError` is raised.
    """
    xml_path = xml_path or Path.cwd()
    if subject_ids:
        xml_files = []
        xml_regex = [f"ADNI_{e}*.xml" for e in subject_ids]
        for subj_files in xml_regex:
            xml_files.extend([f for f in xml_path.glob(subj_files)])
    else:
        xml_files = [f for f in xml_path.glob("Clinica_processed_metadata/ADNI_*.xml")]
    if len(xml_files) == 0:
        msg = f"{LOGGING_HEADER} No ADNI xml files were found for reading the metadata."
        cprint(msg, lvl="error")
        raise FileNotFoundError(msg)
    return xml_files


def _check_xml_tag(element_tag: str, expected_tag: str):
    """Check that the XML element tagged matches the expected tag."""
    if element_tag != expected_tag:
        _log_and_raise_xml_parser_error(
            f"Bad tag: expected {expected_tag}, got {element_tag}"
        )


def _check_xml_nb_children(
    xml_el: xml.etree.ElementTree.Element,
    expected_nb_children: Union[int, Iterable[int]],
):
    """Check that the given XML element has the expected number of children.
    Raise a `ValueError` otherwise.
    """
    nb_children = len(xml_el)
    if isinstance(expected_nb_children, int):
        if nb_children != expected_nb_children:
            _log_and_raise_xml_parser_error(
                f"Bad number of children for <{xml_el.tag}>: "
                f"got {nb_children} != {expected_nb_children}"
            )
    else:
        if nb_children not in expected_nb_children:
            _log_and_raise_xml_parser_error(
                f"Bad number of children for <{xml_el.tag}>: "
                f"got {nb_children}, not in {expected_nb_children}"
            )


def _check_xml(
    xml_el: xml.etree.ElementTree.Element,
    exp_tag: str,
    exp_nb_children: Optional[Union[int, Iterable[int]]] = None,
) -> xml.etree.ElementTree.Element:
    """Check that the given xml element is as expected (tag name + nb of children).
    Raise a ValueError if not and return the element if it is valid.
    """
    _check_xml_tag(xml_el.tag, exp_tag)
    if exp_nb_children is not None:
        _check_xml_nb_children(xml_el, exp_nb_children)
    return xml_el


def _parse_project(
    root: xml.etree.ElementTree.Element,
) -> xml.etree.ElementTree.Element:
    """Check that root has 1 child: <project> having 4 children,
    extract project from them.
    Check that the project identifier contains ADNI.
    """
    if len(root) != 1:
        _log_and_raise_xml_parser_error(
            "XML root should have only one child. " f"{len(root)} children found.",
        )
    project = _check_xml(root[0], "project", 4)
    _check_xml_project_identifier(project, "ADNI")
    return project


def _check_xml_and_get_text(
    xml_element: xml.etree.ElementTree.Element, expected_tag: str, *, cast=None
) -> str:
    """Check xml element and return its text (leaf)"""
    _check_xml(xml_element, expected_tag, 0)  # leaf
    return _get_text(xml_element, cast=cast)


def _get_text(xml_element: xml.etree.ElementTree.Element, cast=None) -> str:
    """Get text of xml element (leaf)."""
    if cast is not None:
        return cast(xml_element.text)
    return xml_element.text


def _check_derived_image_xml(derived: xml.etree.ElementTree.Element) -> None:
    """Perform sanity checks on derived images."""
    if len(derived) < 9:
        _log_and_raise_xml_parser_error("Derived image does not have enough field.")
    for idx, field, expected in zip(
        [2, 3, 4, 5, 6],
        ["imageType", "tissue", "hemisphere", "anatomicStructure", "registration"],
        ["image volume", "All", "Both", "Brain", "native"],
    ):
        _check_xml_field(derived[idx], field, expected)
    idx = -1 if derived[-1].tag == "creationDate" else -2
    if (
        image_create_date := _check_xml_and_get_text(derived[idx], "creationDate")
    ) != "0000-00-00":
        try:
            _validate_date_iso_format(image_create_date)
        except ValueError as e:
            _log_and_raise_xml_parser_error(
                f"The creationDate for the derived image is not valid: {e}"
            )


def _check_xml_field(
    xml_element: xml.etree.ElementTree.Element, field: str, expected_value: str
) -> None:
    """Check that the given field of the given XML element has the expected value."""
    if (found_value := _check_xml_and_get_text(xml_element, field)) != expected_value:
        _log_and_raise_xml_parser_error(
            f"The {field} for the derived image should be '{expected_value}'. "
            f"{found_value} was found instead."
        )


def _validate_date_iso_format(date: str) -> None:
    import datetime

    try:
        datetime.date.fromisoformat(date)
    except ValueError:
        raise ValueError(f"Incorrect data format {date}, should be YYYY-MM-DD")


def _get_derived_image_metadata(derived: xml.etree.ElementTree.Element) -> dict:
    """Return metadata for derived images after
    performing some sanity checks.
    """
    _check_derived_image_xml(derived)
    return {
        "image_proc_id": _check_xml_and_get_text(derived[0], "imageUID", cast=int),
        "image_proc_desc": _check_xml_and_get_text(
            derived[1], "processedDataLabel"
        ).strip(),
    }


def _check_xml_project_identifier(
    project: xml.etree.ElementTree.Element, expected: str
):
    """Check the project identifier."""
    if _check_xml_and_get_text(project[0], "projectIdentifier") != expected:
        _log_and_raise_xml_parser_error(f"Not {expected} cohort")


def _get_original_image_metadata(original_image: xml.etree.ElementTree.Element) -> dict:
    """Get original image metadata."""
    if not isinstance(original_image, list):
        return _get_image_metadata(original_image)
    # only keep common metadata (remaining will be None) and
    original_images_metadata = [_get_image_metadata(_) for _ in original_image]
    original_images_metadata = {
        _.pop("image_orig_id"): _ for _ in original_images_metadata
    }
    original_image_metadata = {
        "image_orig_id": "|".join(map(str, original_images_metadata.keys()))
    }
    for k in next(iter(original_images_metadata.values())).keys():
        # loop on original images metadata and only add consistent metadata
        set_v = set(_[k] for _ in original_images_metadata.values())
        if len(set_v) == 1:
            original_image_metadata[k] = set_v.pop()
    return original_image_metadata


def _check_processed_image_rating(
    processed_image_rating: Union[int, None],
    original_image_metadata: dict,
    derived_image_metadata: dict,
):
    """Check, if possible, that the rating of the processed images are
    consistent with the rating of the original image.
    """
    from clinica.utils.stream import cprint

    rating = original_image_metadata.get("image_orig_rating", None)
    if processed_image_rating is not None:
        if processed_image_rating != rating:
            cprint(
                msg=(
                    f"{LOGGING_HEADER} Image rating for processed image "
                    f"{derived_image_metadata.get('image_proc_id')} not "
                    "consistent with rating of original image"
                ),
                lvl="info",
            )


def _get_root_from_xml_path(xml_path: Path) -> xml.etree.ElementTree.Element:
    """Return the root XML element from the XML file path."""
    from xml.etree import ElementTree

    try:
        tree = ElementTree.parse(xml_path)
        root = tree.getroot()
        return root
    except Exception as e:
        _log_and_raise_xml_parser_error(
            f"Error parsing XML file {xml_path.name} :\n{e}"
        )


def _parse_series(
    study: xml.etree.ElementTree.Element,
) -> xml.etree.ElementTree.Element:
    """Get series xml element from study xml element."""
    expected_length = 3 if len(study) == 7 else 4
    return _check_xml(study[5], "series", expected_length)


def _parse_derived_images(
    series: xml.etree.ElementTree.Element,
) -> Tuple[Union[xml.etree.ElementTree.Element, list], Union[int, None], dict]:
    """Parse the series xml element when derived images are present.
    Return the original image, the rating of the processed images,
    and the derived image metadata as a dict.
    """
    # we can have 1 or 2 related imgs...
    series_meta = _check_xml(series[-1], "seriesLevelMeta", {3, 4})
    related_image = _check_xml(series_meta[2], "relatedImageDetail", 1)
    original_image = _check_xml(related_image[0], "originalRelatedImage", {3, 4})
    if len(series_meta) == 4:
        # 2 original images
        related_image_bis = _check_xml(series_meta[-1], "relatedImageDetail", 1)
        original_image_bis = _check_xml(
            related_image_bis[0], "originalRelatedImage", len(original_image)
        )
        original_image = [original_image, original_image_bis]
    # specific
    pipeline_name = _check_xml_and_get_text(
        _check_xml(series_meta[0], "annotation", 1)[0], "text"
    )
    deriv = _check_xml(series_meta[1], "derivedProduct")
    derived_metadata = {"image_proc_pipe": pipeline_name.strip()}
    derived_metadata.update(_get_derived_image_metadata(deriv))
    processed_image_rating = _get_image_rating(deriv[-1])
    return original_image, processed_image_rating, derived_metadata


def _check_image(img: xml.etree.ElementTree.Element):
    """Check that the image XML element is valid."""
    if img.tag not in {
        "imagingProtocol",  # for original image metadata
        "originalRelatedImage",  # for processed image
    }:
        _log_and_raise_xml_parser_error(
            f"Bad image tag <{img.tag}>. "
            "Should be either 'imagingProtocol' or 'originalRelatedImage'."
        )
    if len(img) not in {3, 4}:
        _log_and_raise_xml_parser_error(
            f"Image XML element has {len(img)} children. " "Expected either 3 or 4.",
        )


def _clean_protocol_metadata(protocol_metadata: dict) -> dict:
    """Replace confusing '|' (reserved for separation btw elements) with space.
    Only 1 manufacturer with this apparently but...
    """
    return {
        k: v.replace("|", " ")
        if isinstance(v, str)
        else v  # there are some None (not str)
        for k, v in protocol_metadata.items()
    }


def _filter_metadata(metadata: dict) -> dict:
    """Filter and clean a given metadata dictionary according to the
    METADATA_NAME_MAPPING dictionary."""
    filtered = dict()
    for k, v in metadata.items():
        if k in METADATA_NAME_MAPPING:
            filtered[METADATA_NAME_MAPPING[k]] = v
    return filtered


def _get_image_metadata(img: xml.etree.ElementTree.Element) -> dict:
    """Return information on original image as dict of metadata."""
    _check_image(img)
    protocol = _check_xml(img[2], "protocolTerm")  # do not check children
    assert all(p.tag == "protocol" and "term" in p.attrib.keys() for p in protocol)
    protocol_metadata = {
        "_".join(p.attrib["term"].split()).lower(): p.text for p in protocol
    }
    protocol_metadata = _clean_protocol_metadata(protocol_metadata)
    filtered_protocol_metadata = _filter_metadata(protocol_metadata)
    return {
        "image_orig_id": _check_xml_and_get_text(img[0], "imageUID", cast=int),
        "image_orig_seq": _check_xml_and_get_text(
            img[1], "description"
        ),  # .replace('|', ' '),
        **filtered_protocol_metadata,
    }


def _get_image_rating(image_rating: xml.etree.ElementTree.Element) -> Optional[str]:
    """Get the image rating value as an integer from the xml element."""
    if image_rating.tag != "imageRating":
        return None
    _check_xml(image_rating, "imageRating", 2)
    image_rating_desc = _check_xml_and_get_text(image_rating[0], "ratingDescription")
    image_rating_val = _check_xml_and_get_text(image_rating[1], "value", cast=int)
    assert image_rating_desc == str(image_rating_val)
    return image_rating_val


def _check_modality(study: xml.etree.ElementTree.Element, expected_modality: str):
    """Check that the modality of the given study is the expected one."""
    series = _parse_series(study)
    modality = _check_xml_and_get_text(series[1], "modality")
    if modality != expected_modality:
        _log_and_raise_xml_parser_error(
            f"Unexpected modality {modality}, expected {expected_modality}."
        )


def _parse_subject(
    project: xml.etree.ElementTree.Element,
) -> Tuple[str, xml.etree.ElementTree.Element]:
    """From the project xml element, parse the subject and subject id."""
    subject = _check_xml(project[-1], "subject", {5, 7})  # with APOE or not
    subject_id = _check_xml_and_get_text(subject[0], "subjectIdentifier")
    return subject_id, subject


def _parse_study(
    subject: xml.etree.ElementTree.Element,
) -> xml.etree.ElementTree.Element:
    """Parse the study from the subject xml element
    and check that it is MRI.
    """
    study = _check_xml(subject[-1], "study", {6, 7})
    _check_modality(study, "MRI")
    return study


def _parse_images(
    study: xml.etree.ElementTree.Element,
) -> Tuple[Union[xml.etree.ElementTree.Element, list], Union[int, None], dict]:
    """Parse the images, original and derived if available."""
    if len(study) == 7:
        # Original only
        original_image = _check_xml(study[6], "imagingProtocol", {3, 4})
        return original_image, None, dict()
    else:
        series = _parse_series(study)
        return _parse_derived_images(series)


def _parse_xml_file(xml_path: Path) -> dict:
    """Parse the given XML file and return the desired metadata as a dict."""
    project = _parse_project(_get_root_from_xml_path(xml_path))
    subject_id, subject = _parse_subject(project)
    study = _parse_study(subject)
    series = _parse_series(study)
    original_image, processed_image_rating, derived_image_metadata = _parse_images(
        study
    )
    original_image_metadata = _get_original_image_metadata(original_image)
    _check_processed_image_rating(
        processed_image_rating, original_image_metadata, derived_image_metadata
    )
    acquisition_date = pd.Timestamp(_check_xml_and_get_text(series[2], "dateAcquired"))
    scan_metadata = {
        "id": subject_id,
        "acq_time": acquisition_date.isoformat(),
        **original_image_metadata,
        **derived_image_metadata,
    }
    if "image_proc_id" not in scan_metadata and "image_orig_id" not in scan_metadata:
        _log_and_raise_xml_parser_error(
            f"Scan metadata for subject {subject_id} has no image ID."
        )
    return scan_metadata


class FuncWithException:
    def __init__(self, f):
        self.f = f

    def __call__(self, *args, **kwargs):
        """Return a tuple: (result, exception)"""
        try:
            return self.f(*args, **kwargs), None
        except Exception as e:
            return None, e


def _run_parsers(xml_files: Iterable[Path]) -> Tuple[Iterable[dict], dict]:
    """Run the parser `_parse_xml_file` on the list of files `xml_files`.
    Returns a tuple consisting if parsed images
    metadata and captured exceptions.
    """
    import os

    parser = FuncWithException(_parse_xml_file)
    images_with_exceptions = dict(
        zip(
            xml_files,  # tqdm is buggy when chunksize > 1
            [parser(xml_file) for xml_file in xml_files],
        )
    )
    images = [
        image
        for image, exception in images_with_exceptions.values()
        if exception is None
    ]
    exceptions = {
        os.path.basename(xml_p): exception
        for xml_p, (_, exception) in images_with_exceptions.items()
        if exception is not None
    }
    return images, exceptions


def _create_mri_meta_df(images: Iterable[dict]) -> pd.DataFrame:
    """Create dataframe from images metadata parsed from the XML files."""
    mri_meta = pd.DataFrame(images).convert_dtypes()
    # Scan ID as expected by Clinica (derived if so, original otherwise)
    mri_meta["T1w_scan_id"] = mri_meta["image_proc_id"].combine_first(
        mri_meta["image_orig_id"]
    )
    # Index
    mri_meta.set_index(
        [
            "T1w_scan_id",
            "id",
            # "GENDER",
            # "VISIT_DESC",
            # "VISIT_AGE",
            # "SCAN_DATE",
            # "STUDY_ID",
            # "SERIES_ID",
        ],
        inplace=True,
    )
    mri_meta.reset_index(inplace=True)
    return mri_meta


def _get_existing_scan_dataframe(
    subject_path: Path, session: str
) -> Tuple[Optional[pd.DataFrame], Path]:
    """Retrieve existing scan dataframe at the given
    `subj_path`, and the given `session`.
    If no existing scan file is found, a warning is given.
    """
    import warnings

    from clinica.utils.stream import cprint

    scans_tsv_path = subject_path / session / f"{subject_path.name}_{session}_scans.tsv"
    if scans_tsv_path.exists():
        df_scans = pd.read_csv(scans_tsv_path, sep="\t")
        df_scans["scan_id"] = df_scans["scan_id"].astype("Int64")
        return df_scans, scans_tsv_path
    msg = f"{LOGGING_HEADER} No scan tsv file for subject {subject_path.name} and session {session}"
    cprint(msg, lvl="warning")
    warnings.warn(msg)
    return None, scans_tsv_path


def _merge_scan_and_metadata(
    df_scans: pd.DataFrame, df_meta: pd.DataFrame, strategy: dict
) -> pd.DataFrame:
    """Perform a merge between the two provided dataframe
    according to the strategy.
    """
    return pd.merge(df_scans, df_meta, **strategy)


def _get_json_filename_from_scan_filename(scan_filename: Path) -> Path:
    """Replace the given `scan_filename` extensions with '.json'."""
    extensions = "".join(scan_filename.suffixes)
    return Path(str(scan_filename).replace(extensions, ".json"))


def _add_json_scan_metadata(
    json_path: Path,
    metadata: Union[dict, str],
    indent: int = 4,
    keep_none: bool = False,
):
    """Load existing metadata at the provided `json_path`, merge with new
    provided `metadata`, and finally write to disk at `json_path`.

    .. warning::
        This overwrites the existing json file.
    """
    import json

    existing_metadata = dict()
    if json_path.exists():
        with open(json_path, "r") as fp:
            existing_metadata = json.load(fp)
    if isinstance(metadata, str):
        metadata = json.loads(metadata)
    filtered_metadata = {
        k: v for k, v in metadata.items() if k in METADATA_NAME_MAPPING.values()
    }
    updated_metadata = {**existing_metadata, **filtered_metadata}
    if not keep_none:
        updated_metadata = {k: v for k, v in updated_metadata.items() if v is not None}
    if len(updated_metadata) > 0:
        with open(json_path, "w") as fp:
            json.dump(updated_metadata, fp, indent=indent)


def _add_metadata_to_scans(df_meta: pd.DataFrame, bids_subjects_paths: Iterable[Path]):
    """Add the metadata to the appropriate tsv and json files."""
    from clinica.iotools.bids_utils import get_bids_sess_list

    merge_strategy = {"how": "left", "left_on": "scan_id", "right_on": "T1w_scan_id"}

    for subject_path in bids_subjects_paths:
        if (sess_list := get_bids_sess_list(subject_path)) is not None:
            for sess in sess_list:
                df_scans, scans_tsv_path = _get_existing_scan_dataframe(
                    subject_path, sess
                )
                if df_scans is not None:
                    columns_to_keep = list(df_scans.columns) + ["acq_time"]
                    df_merged = _merge_scan_and_metadata(
                        df_scans, df_meta, merge_strategy
                    )
                    for _, scan_row in df_merged.iterrows():
                        _add_json_scan_metadata(
                            _get_json_filename_from_scan_filename(
                                subject_path / sess / scan_row["filename"]
                            ),
                            scan_row.to_json(),
                        )
                    df_merged[columns_to_keep].to_csv(
                        scans_tsv_path, sep="\t", index=False
                    )


def create_json_metadata(
    bids_subjects_paths: Iterable[Path], bids_ids: Iterable[str], xml_path: Path
):
    """Create json metadata dictionary and add the metadata to the
    appropriate files in the BIDS hierarchy."""
    loni_ids = [_bids_id_to_loni(bids_id) for bids_id in bids_ids]
    xml_files = _read_xml_files(loni_ids, xml_path)
    images, exe = _run_parsers(xml_files)
    df_meta = _create_mri_meta_df(images)
    _add_metadata_to_scans(df_meta, bids_subjects_paths)
