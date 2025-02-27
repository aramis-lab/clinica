from typing import Optional

import pandas as pd

from clinica.utils.stream import cprint

from .._utils import ADNIStudy

__all__ = ["visits_to_timepoints"]


def visits_to_timepoints(
    subject: str,
    mri_list_subj: pd.DataFrame,
    adnimerge_subj: pd.DataFrame,
    modality: str,
    visit_field: str = "VISIT",
    scandate_field: str = "SCANDATE",
) -> dict:
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
    if modality == "T1":
        mri_list_subj = mri_list_subj[mri_list_subj[visit_field] != "ADNI Baseline"]

    visits = dict()
    unique_visits = list(mri_list_subj[visit_field].unique())
    pending_timepoints = []

    # We try to obtain the corresponding image Visit for a given VISCODE
    for adni_row in adnimerge_subj.iterrows():
        visit = adni_row[1]
        study = ADNIStudy(visit.ORIGPROT)
        preferred_visit_name = _get_preferred_visit_name(study, visit.VISCODE)

        if preferred_visit_name in unique_visits:
            key_preferred_visit = (visit.VISCODE, visit.COLPROT, visit.ORIGPROT)
            if key_preferred_visit not in visits:
                visits[key_preferred_visit] = preferred_visit_name
            elif visits[key_preferred_visit] != preferred_visit_name:
                cprint(
                    f"[{modality}] Subject {subject} has multiple visits for one timepoint.",
                    lvl="info",
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
            cprint(f"No closest visit found for image {image}", lvl="debug")
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
                f"[{modality}] Subject {subject} has multiple visits for one timepoint.",
                lvl="debug",
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
    visits: list[pd.Series],
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
    if len(visits) == 0:
        cprint(
            "No corresponding timepoint in ADNIMERGE for "
            f"subject {subject} in visit {image_visit}",
            lvl="info",
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
        diff = _compute_number_of_days_between(
            closest_visit_date, second_closest_visit_date
        )
        smallest_time_difference = _compute_number_of_days_between(
            image_acquisition_date, closest_visit_date
        )
        if abs((diff / 2.0) - smallest_time_difference) < 30:
            return True
    return False


def _get_closest_and_second_closest_visits(
    visits: list[pd.Series],
    reference_date: str,
) -> tuple[pd.Series, pd.Series, int, int]:
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
        _compute_number_of_days_between(reference_date, visit_date)
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


def _compute_number_of_days_between(date_1: str, date_2: str) -> int:
    """Calculate the days between two dates."""
    from datetime import datetime

    date_1 = datetime.strptime(date_1, "%Y-%m-%d")
    date_2 = datetime.strptime(date_2, "%Y-%m-%d")
    return abs((date_2 - date_1).days)
