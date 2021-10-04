import xml.etree.ElementTree
from typing import Tuple
import pandas as pd


def read_xml_files(subj_ids=[], xml_path="") -> list:
    from glob import glob
    from os import path

    if subj_ids:
        xml_files = []
        xml_regex = [path.join(xml_path, ("ADNI_" + e + "*.xml")) for e in subj_ids]
        for subj_files in xml_regex:
            xml_files.extend(glob(subj_files))
    else:
        xml_files = glob("Clinica_processed_metadata/ADNI_*.xml")

    if len(xml_files) > 0:
        return xml_files
    else:
        raise IndexError("No ADNI xml files were found for reading the metadata.")


def xml_check(
    xml_el: xml.etree.ElementTree.Element, exp_tag: str, exp_nb_children=None
) -> xml.etree.ElementTree.Element:
    """
    Check that xml element is as expected (tag name + nb of children)
    """

    el_tag = xml_el.tag
    assert el_tag == exp_tag, f"Bad tag: expected {exp_tag}, got {el_tag}"

    if exp_nb_children is not None:  # no check otherwise
        nb_children = len(xml_el)
        if isinstance(exp_nb_children, int):
            assert (
                nb_children == exp_nb_children
            ), f"Bad number of children for <{el_tag}>: got {nb_children} != {exp_nb_children}"
        else:
            # container
            assert (
                nb_children in exp_nb_children
            ), f"Bad number of children for <{el_tag}>: got {nb_children}, not in {exp_nb_children}"

    return xml_el


def xml_check_and_get_text(
    xml_el: xml.etree.ElementTree.Element, exp_tag: str, *, cast=None
) -> str:
    """
    Check xml element and return its text (leaf)
    """
    xml_check(xml_el, exp_tag, 0)  # leaf
    val = xml_el.text
    if cast is not None:
        val = cast(val)
    return val


def get_img_rating(img_rating: xml.etree.ElementTree.Element) -> int:

    xml_check(img_rating, "imageRating", 2)

    # value or description seem to be seemed things....
    img_rating_desc = xml_check_and_get_text(img_rating[0], "ratingDescription")
    img_rating_val = xml_check_and_get_text(img_rating[1], "value", cast=int)
    assert img_rating_desc == str(img_rating_val)

    return img_rating_val


def get_img_metadata(img):
    """
    Return info of original image as dict
    """

    assert img.tag in {
        "imagingProtocol",  # for original image metadata
        "originalRelatedImage",  # for processed image
    }
    assert len(img) in {3, 4}  # children check

    proto = xml_check(img[2], "protocolTerm")  # do not check children

    img_rating_val = get_img_rating(img[3]) if len(img) == 4 else None

    assert all(p.tag == "protocol" and "term" in p.attrib.keys() for p in proto)

    proto_d = {
        "MRI_" + "_".join(p.attrib["term"].split()).upper(): p.text for p in proto
    }

    # replace confusing '|' (reserved for separation btw elements) with space (only 1 manufacturer with this apparently but...)
    proto_d = {
        k: v.replace("|", " ")
        if isinstance(v, str)
        else v  # there are some None (not str)
        for k, v in proto_d.items()
    }

    return {
        "IMAGE_ORIG_ID": xml_check_and_get_text(img[0], "imageUID", cast=int),
        "IMAGE_ORIG_SEQ": xml_check_and_get_text(
            img[1], "description"
        ),  # .replace('|', ' '),
        # "IMAGE_ORIG_RATING": img_rating_val,  # desc is same thing
        **proto_d,
    }


def parse_xml_file(xml_path: str) -> dict:
    import os
    import xml.etree.ElementTree as ET
    from clinica.utils.stream import cprint

    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
    except Exception as e:
        raise ValueError(os.path.basename(xml_path) + ": " + str(e))

    # check that root has 1 child: <project> having 4 children
    assert len(root) == 1
    proj = xml_check(root[0], "project", 4)

    # check project identifier
    assert (
        xml_check_and_get_text(proj[0], "projectIdentifier") == "ADNI"
    ), "Not ADNI cohort"
    # site id (can be directly extracted from subj ID rather...)

    subj = xml_check(proj[-1], "subject", {5, 7})  # with APOE or not

    vis = xml_check(subj[-2], "visit", None)  # can have cog scores also here
    # not formatted as a visit code but description of visit instead...
    vis_desc = xml_check_and_get_text(vis[0], "visitIdentifier")

    study = xml_check(subj[-1], "study", {6, 7})

    derived_d = {}
    proc_img_rating = None

    if len(study) == 7:
        # original img
        series = xml_check(study[5], "series", 3)
        orig_img = xml_check(study[6], "imagingProtocol", {3, 4})

    else:
        # derived img
        series = xml_check(study[5], "series", 4)
        series_meta = xml_check(
            series[-1], "seriesLevelMeta", {3, 4}
        )  # we can have 1 or 2 related imgs...

        rel_img = xml_check(series_meta[2], "relatedImageDetail", 1)
        orig_img = xml_check(rel_img[0], "originalRelatedImage", {3, 4})

        if len(series_meta) == 4:
            # 2 original images
            rel_img_bis = xml_check(series_meta[-1], "relatedImageDetail", 1)
            orig_img_bis = xml_check(
                rel_img_bis[0], "originalRelatedImage", len(orig_img)
            )
            orig_img = [orig_img, orig_img_bis]

        # specific
        pipe_name = xml_check_and_get_text(
            xml_check(series_meta[0], "annotation", 1)[0], "text"
        )
        # assert pipe_name == 'Grinder Pipeline', f'Wrong pipeline name: {pipe_name}'

        deriv = xml_check(series_meta[1], "derivedProduct")
        assert len(deriv) >= 9  # multiple <provenanceDetail> nodes possible

        derived_d = {
            "IMAGE_PROC_ID": xml_check_and_get_text(deriv[0], "imageUID", cast=int),
            "IMAGE_PROC_PIPE": pipe_name.strip(),
            "IMAGE_PROC_DESC": xml_check_and_get_text(
                deriv[1], "processedDataLabel"
            ).strip(),
        }

        assert (
            xml_check_and_get_text(deriv[2], "imageType") == "image volume"
        ), "imageType"
        assert xml_check_and_get_text(deriv[3], "tissue") == "All", "tissue"
        assert xml_check_and_get_text(deriv[4], "hemisphere") == "Both", "hemis"
        assert (
            xml_check_and_get_text(deriv[5], "anatomicStructure") == "Brain"
        ), "anat struct"
        assert (
            xml_check_and_get_text(deriv[6], "registration") == "native"
        ), "registration not native"
        # skip check relatedImage, original ID consistence + "derived from"...
        # skip processing steps with names, versions, dates and so on

        # 2 cases
        if deriv[-1].tag == "creationDate":
            assert (
                xml_check_and_get_text(deriv[-1], "creationDate") == "0000-00-00"
            ), "creaDate"
        else:
            assert (
                xml_check_and_get_text(deriv[-2], "creationDate") == "0000-00-00"
            ), "creaDate"
            proc_img_rating = get_img_rating(deriv[-1])

    assert xml_check_and_get_text(series[1], "modality") == "MRI", "Unexpected modality"

    if isinstance(orig_img, list):
        # multiple original images

        # only keep common metadata (remaining will be None) and
        orig_imgs_d = [get_img_metadata(oi) for oi in orig_img]
        orig_imgs_d = {oid.pop("IMAGE_ORIG_ID"): oid for oid in orig_imgs_d}
        orig_img_d = {
            "IMAGE_ORIG_ID": "|".join(map(str, orig_imgs_d.keys())),
        }
        for k in next(iter(orig_imgs_d.values())).keys():
            # loop on original images metadata and only add consistent metadata
            set_v = set(oid[k] for oid in orig_imgs_d.values())
            if len(set_v) == 1:
                orig_img_d[k] = set_v.pop()
            # else:
            # orig_img_d[k] = '|'.join(map(str, set_v)) # keep all values with '|'

        # assert {**get_img_metadata(orig_img), 'IMAGE_ORIG_ID': -1} == \
        #       {**get_img_metadata(orig_img_bis), 'IMAGE_ORIG_ID': -1}, 'Not consistent multiple derived scans'

    else:
        # standard case: only one original image
        orig_img_d = get_img_metadata(orig_img)

    if proc_img_rating is not None:
        if proc_img_rating != orig_img_d.get("IMAGE_ORIG_RATING", None):
            cprint(
                msg=f"Image rating for processed image {derived_d.get('IMAGE_PROC_ID')} not consistent with rating of original image",
                lvl="info",
            )

    row_d = {
        "ID": xml_check_and_get_text(subj[0], "subjectIdentifier"),
        # subject cofactor to check
        # "GENDER": xml_check_and_get_text(subj[2], "subjectSex"),
        # include apoe as well?
        # "VISIT_DESC": vis_desc,  # to join visit ?
        # "VISIT_AGE": xml_check_and_get_text(
        #    study[1], "subjectAge", cast=float
        # ),  # to join visit ?
        # "SCAN_DATE": xml_check_and_get_text(series[2], "dateAcquired"),
        # "STUDY_ID": xml_check_and_get_text(study[0], "studyIdentifier", cast=int),
        # series
        # "SERIES_ID": xml_check_and_get_text(series[0], "seriesIdentifier", cast=int),
        # (original) imagingProtocol
        **orig_img_d,
        # derived image props
        **derived_d,
    }

    return row_d


def create_mri_meta_df(imgs: list) -> pd.DataFrame:

    mri_meta = pd.DataFrame(imgs).convert_dtypes()

    # Scan ID as expected by Clinica (derived if so, original otherwise)
    mri_meta["T1w_scan_id"] = mri_meta["IMAGE_PROC_ID"].combine_first(
        mri_meta["IMAGE_ORIG_ID"]
    )

    # Index
    mri_meta.set_index(
        [
            "T1w_scan_id",
            "ID",
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


class func_with_exception:
    def __init__(self, f):
        self.f = f

    def __call__(self, *args, **kwargs):
        """
        Return a tuple: (result, exception)
        """
        try:
            return self.f(*args, **kwargs), None
        except Exception as e:
            return None, e


def run_parsers(xml_files: list) -> Tuple[list, dict]:

    import os

    # from concurrent.futures import ProcessPoolExecutor, as_completed #ThreadPoolExecutor
    from math import ceil

    parser = func_with_exception(parse_xml_file)
    imgs_with_excep = dict(
        zip(
            xml_files,  # tqdm is buggy when chunksize > 1
            [parser(xml_file) for xml_file in xml_files],
        )
    )

    imgs = [img for img, e in imgs_with_excep.values() if e is None]
    exceps = {
        os.path.basename(xml_p): e
        for xml_p, (_, e) in imgs_with_excep.items()
        if e is not None
    }

    return imgs, exceps


def create_json_metadata(bids_subjs_paths: str, bids_ids: list, xml_path: str) -> None:
    """
    Create json metadata dictionary and add the columns to the appropriate scans.tsv files
    """

    from clinica.iotools.converters.adni_to_bids.adni_utils import bids_id_to_loni

    loni_ids = [bids_id_to_loni(bids_id) for bids_id in bids_ids]
    xml_files = read_xml_files(loni_ids, xml_path)
    imgs, _ = run_parsers(xml_files)
    df_meta = create_mri_meta_df(imgs)

    add_metadata_to_scans(df_meta, bids_subjs_paths)
    return


def add_metadata_to_scans(df_meta, bids_subjs_paths: list) -> None:
    """
    Add the metadata to the appropriate scans.tsv file
    """
    from pathlib import Path
    from clinica.iotools.bids_utils import get_bids_sess_list

    for subj_path in bids_subjs_paths:
        sess_list = get_bids_sess_list(subj_path)
        subj_id = Path(subj_path).name
        if sess_list:
            for sess in sess_list:
                scans_tsv_path = Path(subj_path) / sess / f"{subj_id}_{sess}_scans.tsv"
                if scans_tsv_path.exists():
                    df_scans = pd.read_csv(scans_tsv_path, sep="\t")
                    df_merged = pd.merge(
                        df_scans,
                        df_meta,
                        how="left",
                        left_on="scan_id",
                        right_on="T1w_scan_id",
                    )
                    df_merged.to_csv(scans_tsv_path, sep="\t")
    return
