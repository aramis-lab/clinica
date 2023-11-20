from enum import Enum, auto
from typing import List, Optional, Tuple, Union

import pandas as pd


class ADNIStudy(Enum):
    """Possible versions of ADNI studies."""

    ADNI1 = auto()
    ADNI2 = auto()
    ADNI3 = auto()
    ADNIGO = auto()

    @classmethod
    def from_string(cls, study_name: str):
        if study_name == "ADNI1":
            return cls.ADNI1
        if study_name == "ADNI2":
            return cls.ADNI2
        if study_name == "ADNI3":
            return cls.ADNI3
        if study_name == "ADNIGO":
            return cls.ADNIGO
        raise ValueError(
            f"Invalid study name for ADNI: {study_name}. "
            f"Possible values are {list(ADNIStudy)}"
        )


def visits_to_timepoints(
    subject,
    mri_list_subj,
    adnimerge_subj,
    modality,
    visit_field="VISIT",
    scandate_field="SCANDATE",
):
    """
    A correspondence is established between ADNIMERGE and MRILIST visits, since they have different notations.

    Args:
        subject: Subject identifier
        mri_list_subj: Dataframe containing list of MRI scans for the subject
        adnimerge_subj: Dataframe containing visits data for the subject
        modality: Imaging modality
        visit_field: field name corresponding to the visit
        scandate_field: field name corresponding to the scan date

    Returns:
        [Returns]
    """
    from clinica.utils.stream import cprint

    if modality == "T1":
        mri_list_subj = mri_list_subj[mri_list_subj[visit_field] != "ADNI Baseline"]

    visits = dict()
    unique_visits = list(mri_list_subj[visit_field].unique())
    pending_timepoints = []

    # We try to obtain the corresponding image Visit for a given VISCODE
    for adni_row in adnimerge_subj.iterrows():
        visit = adni_row[1]
        study = ADNIStudy.from_string(visit.ORIGPROT)
        preferred_visit_name = _get_preferred_visit_name(study, visit.VISCODE)

        if preferred_visit_name in unique_visits:
            key_preferred_visit = (visit.VISCODE, visit.COLPROT, visit.ORIGPROT)
            if key_preferred_visit not in visits.keys():
                visits[key_preferred_visit] = preferred_visit_name
            elif visits[key_preferred_visit] != preferred_visit_name:
                cprint(
                    f"[{modality}] Subject {subject} has multiple visits for one timepoint."
                )

            unique_visits.remove(preferred_visit_name)
            continue

        pending_timepoints.append(visit)

    # Then for images.Visit non matching the expected labels we find the closest date in visits list
    for visit in unique_visits:
        image = (mri_list_subj[mri_list_subj[visit_field] == visit]).iloc[0]

        closest_visit = _get_closest_visit(
            image[scandate_field],
            image[visit_field],
            pending_timepoints,
            subject,
        )
        if closest_visit is None:
            cprint(image, lvl="debug")
            continue
        key_min_visit = (
            closest_visit.VISCODE,
            closest_visit.COLPROT,
            closest_visit.ORIGPROT,
        )

        if key_min_visit not in visits.keys():
            visits[key_min_visit] = image[visit_field]
        elif visits[key_min_visit] != image[visit_field]:
            cprint(
                f"[{modality}] Subject {subject} has multiple visits for one timepoint."
            )

    return visits


def _get_preferred_visit_name(study: ADNIStudy, visit_code: str) -> str:
    """Return the expected visit name for a given visit depending
    on the study (ADNI1, ADNI2, ADNI3, ADNIGO).

    Parameters
    ----------
    study : ADNIStudy
        The study for this visit.
    visit_code : str
        The visit code to convert to a new form.

    Returns
    -------
    str :
        The expected visit name.
    """
    if study == ADNIStudy.ADNI3:
        return _get_preferred_visit_name_adni3(visit_code)
    if study == ADNIStudy.ADNI2:
        return _get_preferred_visit_name_adni2(visit_code)
    if study == ADNIStudy.ADNI1:
        return _get_preferred_visit_name_adni1(visit_code)
    if study == ADNIStudy.ADNIGO:
        return _get_preferred_visit_name_adnigo(visit_code)


def _get_preferred_visit_name_adni3(visit_code: str) -> str:
    if visit_code == "bl":
        return "ADNI Screening"
    return f"ADNI3 Year {_parse_year_from_visit_code(visit_code)} Visit"


def _parse_year_from_visit_code(visit_code: str) -> float:
    """Return the year corresponding to the visit code.
    Assumes a visit code of the form 'mXXX'.
    """
    return _parse_month_from_visit_code(visit_code) / 12


def _parse_month_from_visit_code(visit_code: str) -> int:
    """Return the month corresponding to the visit code.
    Assumes a visit code of the form 'mXXX'.
    """
    try:
        return int(visit_code[1:])
    except Exception:
        raise ValueError(
            f"Cannot extract month from visit code {visit_code}."
            "Expected a code of the form 'mXXX' where XXX can be casted to an integer."
        )


def _get_preferred_visit_name_adni2(visit_code: str) -> str:
    if visit_code == "bl":
        return "ADNI2 Screening MRI-New Pt"
    if visit_code == "m03":
        return "ADNI2 Month 3 MRI-New Pt"
    if visit_code == "m06":
        return "ADNI2 Month 6-New Pt"
    return f"ADNI2 Year {_parse_year_from_visit_code(visit_code)} Visit"


def _get_preferred_visit_name_adni1(visit_code: str) -> str:
    if visit_code == "bl":
        return "ADNI Screening"
    if visit_code == "m03":
        return "ADNIGO Month 3 MRI"  # does not make any sense...
    month = _parse_month_from_visit_code(visit_code)
    return f"ADNI{'1/' if month < 54 else ''}GO Month {month}"


def _get_preferred_visit_name_adnigo(visit_code: str) -> str:
    if visit_code == "bl":
        return "ADNIGO Screening MRI"
    return _get_preferred_visit_name_adni1(visit_code)


def _get_closest_visit(
    image_acquisition_date: str,
    image_visit: str,
    visits: List[pd.Series],
    subject: str,
) -> Optional[pd.Series]:
    """Choose the visit with the closest date to a given image acquisition date.

    Parameters
    ----------
    image_acquisition_date : str
        The date in string format when the image was acquired.

    image_visit : str
        The string identifier of the visit.

    visits : list of pd.Series
        List of visit entries among which we have to find the closest one.

    subject : str
        Subject identifier.

    Returns
    -------
    None or pd.Series
        If the list of visits provided is empty, this function returns None.
        Otherwise, it returns a pd.Series corresponding to the closest visit found.
    """
    from clinica.utils.stream import cprint

    if len(visits) == 0:
        cprint(
            "No corresponding timepoint in ADNIMERGE for "
            f"subject {subject} in visit {image_visit}"
        )
        return None

    if len(visits) == 1:
        return visits[0]

    (
        closest_visit,
        second_closest_visit,
        smallest_time_difference,
        second_smallest_time_difference,
    ) = _get_closest_and_second_closest_visits(visits, image_acquisition_date)

    return (
        _get_closest_visit_for_large_time_difference(
            subject,
            image_visit,
            image_acquisition_date,
            closest_visit,
            second_closest_visit,
            smallest_time_difference,
            second_smallest_time_difference,
        )
        if smallest_time_difference > 90
        else closest_visit
    )


def _get_closest_visit_for_large_time_difference(
    subject: str,
    image_visit: str,
    image_acquisition_date: str,
    closest_visit: pd.Series,
    second_closest_visit: pd.Series,
    smallest_time_difference: int,
    second_smallest_time_difference: int,
) -> pd.Series:
    """Special case when the closest visit found is more than 90 days
    apart from the image acquisition date.

    In this situation, if the image acquisition date is in between the
    closest visit and the second-closest visit, and the second-closest
    visit is close enough from the image acquisition date, then select
    the second-closest visit instead of the closest one.
    """
    from clinica.utils.stream import cprint

    cprint(
        _get_smallest_time_difference_too_large_message(
            subject,
            image_visit,
            image_acquisition_date,
            closest_visit,
            second_closest_visit,
            smallest_time_difference,
            second_smallest_time_difference,
        ),
        lvl="debug",
    )
    if _second_closest_visit_is_better(
        closest_visit.EXAMDATE,
        image_acquisition_date,
        second_closest_visit.EXAMDATE,
    ):
        closest_visit = second_closest_visit

    cprint(msg=f"We prefer {closest_visit.VISCODE}", lvl="debug")

    return closest_visit


def _second_closest_visit_is_better(
    closest_visit_date: str,
    image_acquisition_date: str,
    second_closest_visit_date: str,
) -> bool:
    """If image is too close to the date between two visits we prefer the earlier visit."""
    from datetime import datetime

    date_format = "%Y-%m-%d"
    if (
        datetime.strptime(closest_visit_date, date_format)
        > datetime.strptime(image_acquisition_date, date_format)
        > datetime.strptime(second_closest_visit_date, date_format)
    ):
        diff = days_between(closest_visit_date, second_closest_visit_date)
        smallest_time_difference = days_between(
            image_acquisition_date, closest_visit_date
        )
        if abs((diff / 2.0) - smallest_time_difference) < 30:
            return True
    return False


def _get_closest_and_second_closest_visits(
    visits: List[pd.Series],
    reference_date: str,
) -> Tuple[pd.Series, pd.Series, int, int]:
    """Return the closest and second-closest visit from the reference date.

    Parameters
    ----------
    visits : list of pd.Series
        The list of visits to analyze. Should have a length greater than two.

    reference_date : str
        The reference date used to compute the time differences.

    Returns
    -------
    closest_visit : pd.Series
        The visit whose EXAMDATE is the closest to the reference date.

    second_closest_visit : pd.Series
        The second closest visit to the reference date.

    min_time_difference : int
        The number of days between the reference date and the closest visit.

    second_min_time_difference : int
        The number of days between the reference date and the second closest visit.
    """
    import numpy as np

    if len(visits) < 2:
        raise ValueError(f"Expected at least two visits, received {len(visits)}.")
    time_differences_between_reference_and_visits_in_days = [
        days_between(reference_date, visit_date)
        for visit_date in [visit.EXAMDATE for visit in visits]
    ]
    idx = np.argsort(time_differences_between_reference_and_visits_in_days)
    min_time_difference = time_differences_between_reference_and_visits_in_days[idx[0]]
    second_min_time_difference = time_differences_between_reference_and_visits_in_days[
        idx[1]
    ]

    return (
        visits[idx[0]],
        visits[idx[1]],
        min_time_difference,
        second_min_time_difference,
    )


def _get_smallest_time_difference_too_large_message(
    subject: str,
    image_visit: str,
    image_acquisition_date: str,
    closest_visit: pd.Series,
    second_closest_visit: pd.Series,
    smallest_time_difference: int,
    second_smallest_time_difference: int,
) -> str:
    """Build the debug message when the closest found visit is further than 90 days in time."""
    msg = (
        f"More than 90 days for corresponding timepoint in ADNIMERGE for "
        f"subject {subject} in visit {image_visit} on {image_acquisition_date}.\n"
    )
    idx = 1
    for visit, difference in zip(
        [closest_visit, second_closest_visit],
        [smallest_time_difference, second_smallest_time_difference],
    ):
        msg += (
            f"Timepoint {idx}: {visit.VISCODE} - {visit.ORIGPROT} "
            f"on {visit.EXAMDATE} (Distance: {difference} days)\n"
        )
        idx += 1
    return msg


def days_between(d1, d2):
    """Calculate the days between two dates.

    Args:
        d1: date 1
        d2: date 2

    Returns: number of days between date 2 and date 1
    """
    from datetime import datetime

    d1 = datetime.strptime(d1, "%Y-%m-%d")
    d2 = datetime.strptime(d2, "%Y-%m-%d")
    return abs((d2 - d1).days)


def select_image_qc(id_list, mri_qc_subj):
    """Select image from several scans according to QC.

    Args:
        id_list: List of images identifiers to choose among
        mri_qc_subj: Dataframe containing list of QC of scans for the subject

    Returns: int - Chosen image identifier
    """
    import numpy as np

    if len(id_list) == 0:
        return None

    image_ids = ["I" + str(imageuid) for imageuid in id_list]
    int_ids = [int(imageuid) for imageuid in id_list]
    images_qc = mri_qc_subj[mri_qc_subj.loni_image.isin(image_ids)]

    if images_qc.empty:
        return max(int_ids)

    if np.sum(images_qc.series_selected) == 1:
        selected_image = (
            images_qc[images_qc.series_selected == 1].iloc[0].loni_image[1:]
        )
    else:
        images_not_rejected = images_qc[images_qc.series_quality < 4]

        if images_not_rejected.empty:
            # There are no images that passed the qc,
            # so we'll try to see if there are other images without qc.
            # Otherwise, return None.
            qc_ids = set([int(qc_id[1:]) for qc_id in images_qc.loni_image.unique()])
            no_qc_ids = list(set(int_ids) - qc_ids)

            if len(no_qc_ids) == 0:
                return None
            else:
                return max(no_qc_ids)

        # We select the image with the best (lower) QC.
        # If no positive QC available we choose image with -1 (not performed)
        series_quality = [
            q if q > 0 else 4 for q in list(images_not_rejected.series_quality)
        ]
        best_q = np.amin(series_quality)
        if best_q == 4:
            best_q = -1
        images_best_qc = images_not_rejected[
            images_not_rejected.series_quality == best_q
        ]
        if len(images_best_qc) == 1:
            selected_image = images_best_qc.iloc[0].loni_image[1:]
        else:
            best_ids = [int(x[1:]) for x in images_best_qc.loni_image.unique()]
            selected_image = max(best_ids)

    return int(selected_image)


def get_images_pet(
    subject,
    pet_qc_subj,
    subject_pet_meta,
    df_cols,
    modality,
    sequences_preprocessing_step,
    viscode_field="VISCODE2",
):
    """Selection of scans passing QC and at the chosen preprocessing stage is performed.

    Args:
        subject: Subject identifier
        pet_qc_subj: Dataframe containing QC for scans for the subject
        subject_pet_meta: Dataframe containing metadata for scans for the subject
        df_cols: Columns of output dataframe
        modality: Imaging modality
        sequences_preprocessing_step: List of sequence names that correspond to desired preprocessing stage
        viscode_field: Name of the field in the pet_qc_subj dataframe that provides to the visit code

    Returns: Dataframe containing images metadata
    """
    import pandas as pd

    from clinica.utils.stream import cprint

    subj_dfs_list = []

    for visit in list(pet_qc_subj[viscode_field].unique()):
        if pd.isna(visit):
            continue

        pet_qc_visit = pet_qc_subj[pet_qc_subj[viscode_field] == visit]

        if pet_qc_visit.empty:
            continue

        # If there are several scans for a timepoint we start with image acquired last (higher LONIUID)
        pet_qc_visit = pet_qc_visit.sort_values("LONIUID", ascending=False)

        original_pet_meta = pd.DataFrame(columns=subject_pet_meta.columns)
        qc_visit = pet_qc_visit.iloc[0]
        for qc_index in range(len(pet_qc_visit)):
            qc_visit = pet_qc_visit.iloc[qc_index]

            # We are looking for FDG PET metadata of Original images, that passed QC,
            # acquired at the same date as the current scan that passed QC for the current visit,
            # not containing ‘early’ in the sequence name

            original_pet_meta = subject_pet_meta[
                (subject_pet_meta["Orig/Proc"] == "Original")
                & (subject_pet_meta["Image ID"] == int(qc_visit.LONIUID[1:]))
                & (subject_pet_meta["Scan Date"] == qc_visit.EXAMDATE)
                & ~subject_pet_meta.Sequence.str.contains("early", case=False, na=False)
            ]
            # Check if we found a matching image. If yes, we stop looking for it.
            if not original_pet_meta.empty:
                break

        if original_pet_meta.empty:
            cprint(
                f"No {modality} images metadata for subject {subject} and visit {qc_visit[viscode_field]}"
            )
            continue

        original_image = original_pet_meta.iloc[0]

        # Co-registered and Averaged image with the same Series ID of the original image
        averaged_pet_meta = subject_pet_meta[
            subject_pet_meta["Sequence"].isin(sequences_preprocessing_step)
            & (subject_pet_meta["Series ID"] == original_image["Series ID"])
        ]

        # If an explicit Co-registered, Averaged image does not exist,
        # the original image is already in that preprocessing stage.

        if averaged_pet_meta.empty:
            sel_image = original_image
            original = True
        else:
            sel_image = averaged_pet_meta.iloc[0]
            original = False

        phase = "ADNI1" if modality == "PIB-PET" else qc_visit.Phase
        visit = sel_image.Visit
        sequence = replace_sequence_chars(sel_image.Sequence)
        date = sel_image["Scan Date"]
        study_id = sel_image["Study ID"]
        series_id = sel_image["Series ID"]
        image_id = sel_image["Image ID"]

        # If it is an amyloid PET we need to find which is the tracer of the scan and add it to the
        if modality == "Amyloid-PET":
            if "av45" in sel_image.Sequence.lower():
                tracer = "AV45"
            elif "fbb" in sel_image.Sequence.lower():
                tracer = "FBB"
            else:
                cprint(
                    msg=(
                        f"Unknown tracer for Amyloid PET image metadata for subject {subject} "
                        f"for visit {qc_visit[viscode_field]}"
                    ),
                    lvl="warning",
                )
                continue

            scan_data = [
                [
                    phase,
                    subject,
                    qc_visit[viscode_field],
                    str(visit),
                    sequence,
                    date,
                    str(study_id),
                    str(series_id),
                    str(image_id),
                    original,
                    tracer,
                ]
            ]
        else:
            scan_data = [
                [
                    phase,
                    subject,
                    qc_visit[viscode_field],
                    str(visit),
                    sequence,
                    date,
                    str(study_id),
                    str(series_id),
                    str(image_id),
                    original,
                ]
            ]

        row_to_append = pd.DataFrame(scan_data, columns=df_cols)
        subj_dfs_list.append(row_to_append)

    return subj_dfs_list


def correct_diagnosis_sc_adni3(clinical_data_dir, participants_df):
    """Add missing diagnosis at screening in ADNI3 due to changes in DX_bl.

    See https://groups.google.com/g/adni-data/c/whYhKafN_0Q.
    The identification of SMC, EMCI, LMCI is lost in ADNI3.

    Args:
        clinical_data_dir: path to the input folder containing clinical data.
        participants_df: DataFrame containing metadata at the participant level.

    Returns:
        Corrected participants_df.
    """
    from os import path

    import pandas as pd

    from clinica.utils.stream import cprint

    diagnosis_dict = {1: "CN", 2: "MCI", 3: "AD"}
    dxsum_df = load_clinical_csv(clinical_data_dir, "DXSUM_PDXCONV_ADNIALL").set_index(
        ["PTID", "VISCODE2"]
    )
    missing_sc = participants_df[participants_df.original_study == "ADNI3"]
    participants_df.set_index("alternative_id_1", drop=True, inplace=True)
    for alternative_id in missing_sc.alternative_id_1.values:
        try:
            diagnosis_sc = diagnosis_dict[
                dxsum_df.loc[(alternative_id, "sc"), "DIAGNOSIS"].values[0]
            ]
            participants_df.loc[alternative_id, "diagnosis_sc"] = diagnosis_sc
        except KeyError:
            cprint(
                msg=(f"Unknown screening diagnosis for subject {alternative_id}."),
                lvl="warning",
            )
            participants_df.loc[alternative_id, "diagnosis_sc"] = "n/a"

    participants_df.reset_index(inplace=True, drop=False)
    return participants_df


def replace_sequence_chars(sequence_name):
    """Replace special characters in the sequence by underscores (as done for corresponding folder names in ADNI).

    Args:
        sequence_name: sequence to process

    Returns: the new string
    """
    import re

    return re.sub("[ /;*()<>:]", "_", sequence_name)


def write_adni_sessions_tsv(df_subj_sessions, bids_subjs_paths):
    import pandas as pd

    """Write the result of method create_session_dict into several TSV files.

    Args:
        df_subj_sessions: global dataframe containing clinical sessions data for all subjects
        bids_subjs_paths: a list with the path to all bids subjects
    """
    import os
    from os import path

    def compute_amyloid_status(row: pd.DataFrame) -> str:
        if (
            pd.isnull(row["adni_av45"])
            and pd.isnull(row["adni_pib"])
            and isinstance(row["adni_abeta"], float)
        ):
            return "Au"
        elif row["adni_av45"] > 1.1 or row["adni_pib"] > 1.5 or row["adni_abeta"] < 192:
            return "A+"
        elif (
            (pd.isnull(row["adni_av45"]) or row["adni_av45"] < 1.1)
            and (pd.isnull(row["adni_pib"]) or row["adni_pib"] < 1.5)
            and (isinstance(row["adni_abeta"], float) or row["adni_abeta"] > 192)
        ):
            return "A-"
        else:
            return "n/a"

    def compute_ptau_status(row: pd.DataFrame) -> str:
        if pd.isnull(row["adni_ptau"]):
            return "Tu"
        elif row["adni_ptau"] > 23:
            return "T+"
        else:
            return "T-"

    df_subj_sessions["adas_memory"] = (
        df_subj_sessions["adas_Q1"]
        + df_subj_sessions["adas_Q4"]
        + df_subj_sessions["adas_Q7"]
        + df_subj_sessions["adas_Q8"]
        + df_subj_sessions["adas_Q9"]
    )  # / 45
    df_subj_sessions["adas_language"] = (
        df_subj_sessions["adas_Q2"]
        + df_subj_sessions["adas_Q5"]
        + df_subj_sessions["adas_Q10"]
        + df_subj_sessions["adas_Q11"]
        + df_subj_sessions["adas_Q12"]
    )  # / 25
    df_subj_sessions["adas_praxis"] = (
        df_subj_sessions["adas_Q3"] + df_subj_sessions["adas_Q6"]
    )  # / 10
    df_subj_sessions["adas_concentration"] = df_subj_sessions["adas_Q13"]  # / 5

    df_subj_sessions = df_subj_sessions.fillna("n/a")

    # compute the amyloid and ptau status
    df_subj_sessions["a_stat"] = df_subj_sessions.apply(
        lambda x: compute_amyloid_status(
            x[["adni_av45", "adni_pib", "adni_abeta"]].apply(
                pd.to_numeric, errors="coerce"
            )
        ),
        axis=1,
    )
    df_subj_sessions["tau_stat"] = df_subj_sessions.apply(
        lambda x: compute_ptau_status(
            x[["adni_ptau"]].apply(pd.to_numeric, errors="coerce")
        ),
        axis=1,
    )

    for sp in bids_subjs_paths:
        os.makedirs(sp, exist_ok=True)

        bids_id = sp.split(os.sep)[-1]

        if bids_id == "conversion_info":
            continue
        else:
            df_tmp = df_subj_sessions[df_subj_sessions["RID"] == bids_id[12:]]

            df_tmp.to_csv(
                path.join(sp, f"{bids_id}_sessions.tsv"),
                sep="\t",
                index=False,
                encoding="utf-8",
            )


def remove_fields_duplicated(bids_fields):
    """Remove duplicated fields in a list.

    Args:
        bids_fields: List of fields

    Returns: list
    """
    seen = set()
    seen_add = seen.add
    return [x for x in bids_fields if not (x in seen or seen_add(x))]


def bids_id_to_loni(bids_id: str) -> Union[str, None]:
    """Convert a subject id of the form sub-ADNI000S0000
    back to original format 000_S_0000
    """
    import re

    ids = re.findall(r"\d+", bids_id)
    if len(ids) == 2:
        return ids[0] + "_S_" + ids[1]
    return None


def filter_subj_bids(df_files, location, bids_ids):
    import clinica.iotools.bids_utils as bids

    # Depending on the file that needs to be open, identify and
    # preprocess the column that contains the subjects ids.
    bids_ids = [x[8:] for x in bids_ids if "sub-ADNI" in x]
    if location == "ADNIMERGE.csv":
        df_files["RID"] = df_files["PTID"].apply(
            lambda x: bids.remove_space_and_symbols(x)[4:]
        )
    else:
        df_files["RID"] = df_files["RID"].apply(lambda x: pad_id(x))
    return df_files.loc[df_files["RID"].isin([x[4:] for x in bids_ids])]


def update_age(row):
    """Update age with time passed since bl to current visit"""
    from datetime import datetime

    if row["session_id"] != "ses-M000":
        examdate = datetime.strptime(row["EXAMDATE"], "%Y-%m-%d")
        examdate_bl = datetime.strptime(row["EXAMDATE_bl"], "%Y-%m-%d")
        delta = examdate - examdate_bl

        updated_age = round(
            float(row["AGE"]) + (delta.days / 365.25),
            1,
        )
    else:
        updated_age = row["AGE"]
    return updated_age


def pad_id(ref_id):
    """Add leading zeros to the RID to keep à length of 4"""
    rid = str(ref_id)
    # Fill the rid with the needed number of zero
    if 4 - len(rid) > 0:
        zeros_to_add = 4 - len(rid)
        rid = "0" * zeros_to_add + rid
    return rid


def _compute_session_id(df: pd.DataFrame, csv_filename: str) -> pd.DataFrame:
    """Compute session ID from visit code.

    The CSV file name is used to determine in which column of
    the dataframe the visit code should be found.
    """
    visit_code_column = _get_visit_code_column_name(csv_filename)
    if visit_code_column not in df:
        raise ValueError(
            f"DataFrame does not contain a column named '{visit_code_column}', "
            "which is supposed to encode the visit code. The columns present in "
            f"the DataFrame are: {df.columns}."
        )
    return df.assign(
        session_id=lambda _df: _df[visit_code_column].apply(
            lambda x: _get_session_id_from_visit_code(x)
        )
    )


def _get_visit_code_column_name(csv_filename: str) -> str:
    """Return the name of the column containing the visit code."""
    if _is_a_visit_code_2_type(csv_filename):
        return "VISCODE2"
    if _is_a_time_point_type(csv_filename):
        return "Timepoint"
    return "VISCODE"


def _is_a_visit_code_2_type(csv_filename: str) -> bool:
    """If the csv file is among these files, then the visit code column is 'VISCODE2'."""
    return csv_filename in {
        "ADAS_ADNIGO2.csv",
        "DXSUM_PDXCONV_ADNIALL.csv",
        "CDR.csv",
        "NEUROBAT.csv",
        "GDSCALE.csv",
        "MODHACH.csv",
        "MOCA.csv",
        "NPIQ.csv",
        "MEDHIST.csv",
        "VITALS.csv",
        "UWNPSYCHSUM_03_07_19.csv",
        "UWNPSYCHSUM_03_26_20.csv",
        "ECOGPT.csv",
        "ECOGSP.csv",
        "FCI.csv",
        "CCI.csv",
        "NPIQ.csv",
        "NPI.csv",
    }


def _is_a_time_point_type(csv_filename: str) -> bool:
    """If the csv file is among these files, then the visit code column is 'Timepoint'."""
    return csv_filename in {
        "BHR_EVERYDAY_COGNITION.csv",
        "BHR_BASELINE_QUESTIONNAIRE.csv",
        "BHR_LONGITUDINAL_QUESTIONNAIRE.csv",
    }


def _get_session_id_from_visit_code(visit_code: str) -> Optional[str]:
    """Checks that the visit code found in the data is supported by
    the converter utility function `viscode_to_session` before using it.
    There are a few special cases that are handled here.
    """
    from clinica.iotools.converter_utils import viscode_to_session

    if _is_visit_code_not_supported(visit_code):
        return None
    if visit_code == "sc":
        return "sc"
    return viscode_to_session(visit_code)


def _is_visit_code_not_supported(visit_code: str) -> bool:
    """Return True is the visit code is not supported by the converter.

    There are a few known values like "f" or "uns1" which are present in
    ADNI data that are not supported for a mapping to a session ID.
    """
    unsupported_values = {"f", "uns1"}
    return pd.isnull(visit_code) or visit_code in unsupported_values


def create_adni_sessions_dict(
    bids_ids, clinic_specs_path, clinical_data_dir, bids_subjs_paths
):
    """Extract all the data required for the sessions files and organize them in a dictionary.

    Args:
        bids_ids: list of subject IDs to process
        clinic_specs_path: path to the specifications for converting the clinical data
        clinical_data_dir: path to the clinical data folder
        bids_subjs_paths: a list with the path to all the BIDS subjects
    """

    import pandas as pd

    from clinica.utils.stream import cprint

    # Load data
    df_sessions = pd.read_csv(clinic_specs_path + "_sessions.tsv", sep="\t")

    files = list(df_sessions["ADNI location"].unique())
    files = [x for x in files if not pd.isnull(x)]
    df_subj_session = pd.DataFrame()
    # write line to get field_bids = sessions['BIDS CLINICA'] without the null values

    # Iterate over the metadata files
    for location in files:
        location = location.split("/")[0]
        try:
            df_file = load_clinical_csv(clinical_data_dir, location.split(".")[0])
        except IOError:
            continue
        df_filtered = filter_subj_bids(df_file, location, bids_ids).copy()
        if not df_filtered.empty:
            df_filtered = _compute_session_id(df_filtered, location)

            # Filter rows with invalid session IDs.
            df_filtered.dropna(subset="session_id", inplace=True)

            if location == "ADNIMERGE.csv":
                df_filtered["AGE"] = df_filtered.apply(lambda x: update_age(x), axis=1)
            df_subj_session = update_sessions_df(
                df_subj_session, df_filtered, df_sessions, location
            )
        else:
            cprint(
                f"Clinical dataframe extracted from {location} is empty after filtering."
            )
            dict_column_correspondence = dict(
                zip(df_sessions["ADNI"], df_sessions["BIDS CLINICA"])
            )
            df_filtered.rename(columns=dict_column_correspondence, inplace=True)
            df_filtered = df_filtered.loc[
                :, (~df_filtered.columns.isin(df_subj_session.columns))
            ]
            df_subj_session = pd.concat([df_subj_session, df_filtered], axis=1)
    if df_subj_session.empty:
        raise ValueError("Empty dataset detected. Clinical data cannot be extracted.")

    # Nv/None refer to sessions whose session is undefined. "sc" is the screening session with unreliable (incomplete)
    # data.
    df_subj_session = df_subj_session[
        df_subj_session.session_id.str.contains(r"ses-M\d+")
    ]
    write_adni_sessions_tsv(df_subj_session, bids_subjs_paths)


def update_sessions_df(df_subj_session, df_filtered, df_sessions, location):
    """Update the sessions dataframe with data of current subject.

    Args:
        df_subj_session: dataframe containing aggregate sessions
        df_filtered: dataframe with current subject sessions data
        df_sessions: dataframe with the the metadata information
        location: the clinica data filename
    """

    import pandas as pd

    df_columns_to_add = df_sessions[
        df_sessions["ADNI location"].str.contains(location, na=False)
    ][["BIDS CLINICA", "ADNI"]]

    df_columns_to_add = df_columns_to_add[
        (df_columns_to_add["ADNI"].isin(df_filtered.columns))
        & (~df_columns_to_add["BIDS CLINICA"].isin(df_subj_session.columns))
    ]

    df_temp = df_filtered[
        ["RID", "session_id"] + list(df_columns_to_add["ADNI"])
    ].copy()
    df_temp.columns = ["RID", "session_id"] + list(df_columns_to_add["BIDS CLINICA"])

    # if error in adni data (duplicate session id), keep only the first row
    df_temp.drop_duplicates(subset=["RID", "session_id"], keep="first", inplace=True)

    try:
        df_temp["diagnosis"] = df_temp["diagnosis"].apply(
            lambda x: convert_diagnosis_code(x)
        )
    except KeyError:
        pass

    if df_subj_session.empty:
        df_subj_session = df_temp
    else:
        df_subj_session = pd.merge(
            df_temp, df_subj_session, on=["RID", "session_id"], how="outer"
        )
    return df_subj_session


def convert_diagnosis_code(diagnosis_code):
    """Convert the string field DX contained in ADNIMERGE.csv into a code that identify the diagnosis.

    Args:
        diagnosis_code: a string that represents the diagnosis

    Returns:
        A code that identify a diagnosis
    """
    from pandas import isna

    # Legend
    diagnosis = {"CN": "CN", "MCI": "MCI", "Dementia": "AD"}

    if isna(diagnosis_code):
        return diagnosis_code
    else:
        return diagnosis[diagnosis_code]


def create_adni_scans_files(conversion_path, bids_subjs_paths):
    """Create scans.tsv files for ADNI.

    Args:
        conversion_path: path to the conversion_info folder
        bids_subjs_paths: list of bids subject paths
    """
    import os
    from glob import glob
    from os import path
    from os.path import normpath

    import pandas as pd

    from clinica.utils.stream import cprint

    scans_fields_bids = ["filename", "scan_id", "mri_field"]

    conversion_versions = [
        folder
        for folder in os.listdir(conversion_path)
        if folder[0] == "v" and path.isdir(path.join(conversion_path, folder))
    ]
    conversion_versions = sorted(
        conversion_versions, key=lambda x: int(x.split("v")[1])
    )
    older_version = conversion_versions[-1]
    converted_dict = dict()
    for tsv_path in os.listdir(path.join(conversion_path, older_version)):
        modality = tsv_path.split("_paths")[0]
        df = pd.read_csv(path.join(conversion_path, older_version, tsv_path), sep="\t")
        df.set_index(["Subject_ID", "VISCODE"], inplace=True, drop=True)
        converted_dict[modality] = df

    for bids_subj_path in bids_subjs_paths:
        # Create the file
        bids_id = os.path.basename(normpath(bids_subj_path))
        subject_id = "_S_".join(bids_id[8::].split("S"))

        sessions_paths = glob(path.join(bids_subj_path, "ses-*"))
        for session_path in sessions_paths:
            session_name = session_path.split(os.sep)[-1]
            viscode = session_label_to_viscode(session_name[4::])
            tsv_name = f"{bids_id}_{session_name}_scans.tsv"

            # If the file already exists, remove it
            if os.path.exists(path.join(session_path, tsv_name)):
                os.remove(path.join(session_path, tsv_name))

            scans_tsv = open(path.join(session_path, tsv_name), "a")
            scans_df = pd.DataFrame(columns=scans_fields_bids)
            scans_df.to_csv(scans_tsv, sep="\t", index=False, encoding="utf-8")

            # Extract modalities available for each subject
            mod_available = glob(path.join(session_path, "*"))
            for mod in mod_available:
                mod_name = os.path.basename(mod)
                files = glob(path.join(mod, "*"))
                for file in files:
                    scans_df = pd.DataFrame(index=[0], columns=scans_fields_bids)
                    file_name = os.path.basename(file)
                    scans_df["filename"] = path.join(mod_name, file_name)
                    converted_mod = find_conversion_mod(file_name)
                    conversion_df = converted_dict[converted_mod]
                    try:
                        scan_id = conversion_df.loc[(subject_id, viscode), "Image_ID"]
                        scans_df["scan_id"] = scan_id
                        if "Field_Strength" in conversion_df.columns.values:
                            field_strength = conversion_df.loc[
                                (subject_id, viscode), "Field_Strength"
                            ]
                            scans_df["mri_field"] = field_strength
                    except KeyError:
                        cprint(
                            msg=f"No information found for file {file_name}",
                            lvl="warning",
                        )

                    scans_df = scans_df.fillna("n/a")
                    scans_df.to_csv(
                        scans_tsv, header=False, sep="\t", index=False, encoding="utf-8"
                    )


def find_conversion_mod(file_name):
    from clinica.utils.pet import Tracer

    suffix = file_name.split("_")[-1].split(".")[0]
    if suffix == "T1w":
        return "t1"
    elif suffix == "FLAIR":
        return "flair"
    elif suffix == "dwi":
        return "dwi"
    elif suffix == "bold":
        return "fmri"
    elif suffix == "pet":
        tracer = file_name.split("trc-")[1].split("_")[0]
        if tracer in (Tracer.AV45, Tracer.FBB):
            return "amyloid_pet"
        else:
            return f"{tracer}_pet"
    else:
        raise ValueError(f"Conversion modality could not be found for file {file_name}")


def find_image_path(images, source_dir, modality, prefix, id_field):
    """
    For each image, the path to an existing image file or folder is created from image metadata.

    This function adds two columns to the input dataframe: 'Is_Dicom', and 'Path'.

    Args:
        images: List of images metadata
        source_dir: path to the ADNI directory
        modality: Imaging modality
        prefix: Prefix to prepend to image identifier to create name of folder containing the image
        id_field: Name of the field in images metadata dataframe containing the image identifier to use

    Returns: Dataframe containing metadata and existing paths
    """
    from pathlib import Path

    import pandas as pd

    from clinica.utils.stream import cprint

    is_dicom = []
    image_folders = []

    for _, image in images.iterrows():
        # Base directory where to look for image files.
        seq_path = Path(source_dir) / str(image["Subject_ID"])

        # Generator finding files containing the image ID.
        find_file = filter(
            lambda p: p.is_file(),
            seq_path.rglob(f"*_I{image['Image_ID']}.*"),
        )

        try:
            # Grab the first file from the generator.
            next_file: Path = next(find_file)

            # Whether the image data are DICOM or not.
            dicom = "dcm" in next_file.suffix

            # Compute the image path (DICOM folder or NIfTI path).
            image_path = str(next_file.parent if dicom else next_file)
        except StopIteration:
            # In case no file is found.
            image_path = ""
            dicom = True

        is_dicom.append(dicom)
        image_folders.append(image_path)
        if image_path == "":
            cprint(
                msg=(
                    f"No {modality} image path found for subject {image.Subject_ID} in visit {image.VISCODE} "
                    f"with image ID {image.Image_ID}"
                ),
                lvl="info",
            )

    images.loc[:, "Is_Dicom"] = pd.Series(is_dicom, index=images.index)
    images.loc[:, "Path"] = pd.Series(image_folders, index=images.index)

    return images


def paths_to_bids(
    images, bids_dir, modality, mod_to_update=False, n_procs: Optional[int] = 1
):
    """Images in the list are converted and copied to directory in BIDS format.

    Args:
        images: List of images metadata and paths
        bids_dir: Path to the output BIDS directory
        modality: Imaging modality
        mod_to_update: If True, pre-existing images in the BIDS directory will be erased and extracted agai n.
        n_procs: int, optional. The number of processes to use for conversion. Default=1.
    """
    from functools import partial
    from multiprocessing import Pool

    if modality.lower() not in [
        "t1",
        "dwi",
        "flair",
        "fmri",
        "fdg",
        "fdg_uniform",
        "pib",
        "av45_fbb",
        "tau",
    ]:
        # This should never be reached
        raise RuntimeError(
            modality.lower() + " is not supported for conversion in paths_to_bids"
        )

    images_list = list([data for _, data in images.iterrows()])

    with Pool(processes=n_procs) as pool:
        create_file_ = partial(
            create_file,
            modality=modality,
            bids_dir=bids_dir,
            mod_to_update=mod_to_update,
        )
        output_file_treated = pool.map(create_file_, images_list)

    return output_file_treated


def create_file(image, modality, bids_dir, mod_to_update):
    """
    Image file is created at the corresponding output folder
    as result of image conversion (DICOM to NIFTI) and centering,
    as specified for each modality

    Args:
        image: Image metadata
        modality: Imaging modality
        bids_dir: Path to the output BIDS directory
        mod_to_update: If True, pre-existing images in the BIDS directory will be erased and extracted again.

    Returns:
        [Returns]
    """
    import os
    import re
    import shutil
    from glob import glob
    from os import path

    from numpy import nan

    from clinica.iotools.bids_utils import run_dcm2niix
    from clinica.iotools.converter_utils import viscode_to_session
    from clinica.iotools.utils.data_handling import center_nifti_origin
    from clinica.utils.pet import ReconstructionMethod, Tracer
    from clinica.utils.stream import cprint

    modality_specific = {
        "t1": {
            "output_path": "anat",
            "output_filename": "_T1w",
            "to_center": True,
            "json": "n",
        },
        "dwi": {
            "output_path": "dwi",
            "output_filename": "_dwi",
            "to_center": False,
            "json": "y",
        },
        "flair": {
            "output_path": "anat",
            "output_filename": "_FLAIR",
            "to_center": True,
            "json": "y",
        },
        "fmri": {
            "output_path": "func",
            "output_filename": "_task-rest_bold",
            "to_center": False,
            "json": "y",
        },
        "fdg": {
            "output_path": "pet",
            "output_filename": f"_trc-{Tracer.FDG}_rec-{ReconstructionMethod.CO_REGISTERED_AVERAGED}_pet",
            "to_center": True,
            "json": "n",
        },
        "fdg_uniform": {
            "output_path": "pet",
            "output_filename": f"_trc-{Tracer.FDG}_rec-{ReconstructionMethod.COREGISTERED_ISOTROPIC}_pet",
            "to_center": True,
            "json": "n",
        },
        "pib": {
            "output_path": "pet",
            "output_filename": f"_trc-{Tracer.PIB}_pet",
            "to_center": True,
            "json": "n",
        },
        "av45": {
            "output_path": "pet",
            "output_filename": f"_trc-{Tracer.AV45}_pet",
            "to_center": True,
            "json": "n",
        },
        "fbb": {
            "output_path": "pet",
            "output_filename": f"_trc-{Tracer.FBB}_pet",
            "to_center": True,
            "json": "n",
        },
        "tau": {
            "output_path": "pet",
            "output_filename": f"_trc-{Tracer.AV1451}_pet",
            "to_center": True,
            "json": "n",
        },
    }

    subject = image.Subject_ID
    viscode = image.VISCODE

    if modality == "av45_fbb":
        modality = image.Tracer.lower()

    if not image.Path:
        cprint(
            msg=(
                f"[{modality.upper()}] No path specified for {subject} "
                f"in session {viscode}"
            ),
            lvl="info",
        )
        return nan
    else:
        cprint(
            msg=(
                f"[{modality.upper()}] Processing subject {subject}"
                f"in session {viscode}"
            ),
            lvl="info",
        )

    session = viscode_to_session(viscode)

    image_path = image.Path
    image_id = image.Image_ID
    # If the original image is a DICOM, check if contains two DICOM
    # inside the same folder
    if image.Is_Dicom:
        image_path = check_two_dcm_folder(image_path, bids_dir, image_id)
    bids_subj = subject.replace("_", "")
    output_path = path.join(
        bids_dir,
        "sub-ADNI" + bids_subj,
        session,
        modality_specific[modality]["output_path"],
    )
    output_filename = (
        "sub-ADNI"
        + bids_subj
        + "_"
        + session
        + modality_specific[modality]["output_filename"]
    )

    # If updated mode is selected, check if an old image is existing and remove it
    existing_im = glob(path.join(output_path, output_filename + "*"))
    if mod_to_update and len(existing_im) > 0:
        print("Removing the old image...")
        for im in existing_im:
            os.remove(im)

    if mod_to_update or not len(existing_im) > 0:
        try:
            os.makedirs(output_path)
        except OSError:
            # Folder already created with previous instance
            pass

        generate_json = modality_specific[modality]["json"]
        if modality_specific[modality]["to_center"]:
            zip_image = "n"
        else:
            zip_image = "y"

        if image.Is_Dicom:
            run_dcm2niix(
                input_dir=image_path,
                output_dir=output_path,
                output_fmt=output_filename,
                compress=True if zip_image == "y" else False,
                bids_sidecar=True if generate_json == "y" else False,
            )

            # If "_t" - the trigger delay time - exists in dcm2niix output filename, we remove it
            exception_t = glob(path.join(output_path, output_filename + "_t[0-9]*"))
            for trigger_time in exception_t:
                res = re.search(r"_t\d+\.", trigger_time)
                no_trigger_time = trigger_time.replace(
                    trigger_time[res.start() : res.end()], "."
                )
                os.rename(trigger_time, no_trigger_time)

            # Removing ADC images if one is generated
            adc_image = path.join(output_path, output_filename + "_ADC.nii.gz")
            if os.path.exists(adc_image):
                os.remove(adc_image)

            nifti_file = path.join(output_path, output_filename + ".nii")
            output_image = nifti_file + ".gz"

            # Conditions to check if output NIFTI files exists,
            # and, if DWI, if .bvec and .bval files are also present
            nifti_exists = path.isfile(nifti_file) or path.isfile(output_image)
            dwi_bvec_and_bval_exist = not (modality == "dwi") or (
                path.isfile(path.join(output_path, output_filename + ".bvec"))
                and path.isfile(path.join(output_path, output_filename + ".bval"))
            )

            # Check if conversion worked (output files exist)
            if not nifti_exists or not dwi_bvec_and_bval_exist:
                cprint(
                    msg=f"Conversion with dcm2niix failed for subject {subject} and session {session}",
                    lvl="warning",
                )
                return nan

            # Case when JSON file was expected, but not generated by dcm2niix
            elif generate_json == "y" and not os.path.exists(
                path.join(output_path, output_filename + ".json")
            ):
                cprint(
                    msg=f"JSON file not generated by dcm2niix for subject {subject} and session {session}",
                    lvl="warning",
                )

            if modality_specific[modality]["to_center"]:
                center_nifti_origin(nifti_file, output_image)
                os.remove(nifti_file)

        else:
            output_image = path.join(output_path, output_filename + ".nii.gz")
            if modality_specific[modality]["to_center"]:
                center_nifti_origin(image_path, output_image)
                if not output_image:
                    cprint(
                        msg=(
                            f"For subject {subject} in session {session}, "
                            f"an error occurred whilst recentering Nifti image: {image_path}"
                        ),
                        lvl="error",
                    )
            else:
                shutil.copy(image_path, output_image)

        # Check if there is still the folder tmp_dcm_folder and remove it
        remove_tmp_dmc_folder(bids_dir, image_id)
        return output_image

    else:
        return nan


def session_label_to_viscode(session_name: str) -> str:
    """Replace the session name passed as input with the session label 'bl' or 'mXXX'.

    Parameters
    ----------
    session_name: str
        Name of the session (MXXX).

    Returns
    -------
    str:
        'bl' if is the baseline session or the original session name.
    """
    if session_name == "M000":
        return "bl"
    else:
        return f"m{(int(session_name[1:])):02d}"


def check_two_dcm_folder(dicom_path, bids_folder, image_uid):
    """[summary].

    Check if a folder contains more than one DICOM and if yes, copy the DICOM related to
    image id passed as parameter into a temporary folder called tmp_dicom_folder.

    Args:
        dicom_path (str): path to the DICOM folder
        bids_folder (str): path to the BIDS folder where the dataset will be stored
        image_uid (str): image id of the fMRI

    Returns:
        str: path to the original DICOM folder or the path to a temporary folder called
            tmp_dicom_folder where only the DICOM to convert is copied
    """
    import os
    import shutil
    from glob import glob
    from os import path
    from shutil import copy

    temp_folder_name = f"tmp_dcm_folder_{str(image_uid).strip(' ')}"
    dest_path = path.join(bids_folder, temp_folder_name)

    # Check if there are dicom files inside the folder not belonging to the desired image
    dicom_list = glob(path.join(dicom_path, "*.dcm"))
    image_list = glob(path.join(dicom_path, f"*{image_uid}.dcm"))
    if len(dicom_list) != len(image_list):
        # Remove the precedent tmp_dcm_folder if present.
        if os.path.exists(dest_path):
            shutil.rmtree(dest_path)
        os.mkdir(dest_path)
        dmc_to_conv = glob(path.join(dicom_path, f"*{str(image_uid)}.dcm*"))
        for d in dmc_to_conv:
            copy(d, dest_path)
        return dest_path
    else:
        return dicom_path


def remove_tmp_dmc_folder(bids_dir, image_id):
    """Remove the folder tmp_dmc_folder created by the method check_two_dcm_folder (if existing).

    Args:
        bids_dir: path to the BIDS directory
        image_id:
    """
    from os.path import exists, join
    from shutil import rmtree

    tmp_dcm_folder_path = join(bids_dir, f"tmp_dcm_folder_{str(image_id).strip(' ')}")
    if exists(tmp_dcm_folder_path):
        rmtree(tmp_dcm_folder_path)


def load_clinical_csv(clinical_dir: str, filename: str) -> pd.DataFrame:
    """Load the clinical csv from ADNI. This function is able to find the csv in the
    different known format available, the old format with just the name, and the new
    format with the name and the date of download.

    Parameters
    ----------
    clinical_dir: str
        Directory containing the csv.

    filename: str
        name of the file without the suffix.

    Returns
    -------
    pd.DataFrame:
        Dataframe corresponding to the filename.
    """
    import re
    from pathlib import Path

    pattern = filename + "(_\d{1,2}[A-Za-z]{3}\d{4})?.csv"
    files_matching_pattern = [
        f for f in Path(clinical_dir).rglob("*.csv") if re.search(pattern, (f.name))
    ]
    if len(files_matching_pattern) != 1:
        raise IOError(
            f"Expecting to find exactly one file in folder {clinical_dir} "
            f"matching pattern {pattern}. {len(files_matching_pattern)} "
            f"files were found instead : \n{'- '.join(str(files_matching_pattern))}"
        )
    try:
        return pd.read_csv(files_matching_pattern[0], sep=",", low_memory=False)
    except:
        raise ValueError(
            f"File {str(files_matching_pattern[0])} was found but could not "
            "be loaded as a DataFrame. Please check your data."
        )
