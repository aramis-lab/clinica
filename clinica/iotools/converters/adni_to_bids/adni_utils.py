# coding: utf-8


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

        preferred_visit_name = get_preferred_visit_name(visit)

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

        key_min_visit = get_closest_visit(
            image, pending_timepoints, subject, visit_field, scandate_field
        )

        if not key_min_visit:
            continue

        if key_min_visit not in visits.keys():
            visits[key_min_visit] = image[visit_field]
        elif visits[key_min_visit] != image[visit_field]:
            cprint(
                f"[{modality}] Subject {subject} has multiple visits for one timepoint."
            )

    return visits


def get_preferred_visit_name(visit):
    """Provide the expected visit name for a given visit.

    Args:
        visit: A visit entry from ADNIMERGE

    Returns:
        string: expected visit name
    """
    if visit.ORIGPROT == "ADNI3":
        if visit.VISCODE == "bl":
            preferred_visit_name = "ADNI Screening"
        else:
            year = str(int(visit.VISCODE[1:]) / 12)
            preferred_visit_name = f"ADNI3 Year {year} Visit"
    elif visit.ORIGPROT == "ADNI2":
        if visit.VISCODE == "bl":
            preferred_visit_name = "ADNI2 Screening MRI-New Pt"
        elif visit.VISCODE == "m03":
            preferred_visit_name = "ADNI2 Month 3 MRI-New Pt"
        elif visit.VISCODE == "m06":
            preferred_visit_name = "ADNI2 Month 6-New Pt"
        else:
            year = str(int(visit.VISCODE[1:]) / 12)
            preferred_visit_name = f"ADNI2 Year {year} Visit"
    else:
        if visit.VISCODE == "bl":
            if visit.ORIGPROT == "ADNI1":
                preferred_visit_name = "ADNI Screening"
            else:  # ADNIGO
                preferred_visit_name = "ADNIGO Screening MRI"
        elif visit.VISCODE == "m03":  # Only for ADNIGO Month 3
            preferred_visit_name = "ADNIGO Month 3 MRI"
        else:
            month = int(visit.VISCODE[1:])
            if month < 54:
                preferred_visit_name = f"ADNI1/GO Month {str(month)}"
            else:
                preferred_visit_name = f"ADNIGO Month {str(month)}"
    return preferred_visit_name


def get_closest_visit(image, pending_timepoints, subject, visit_field, scandate_field):
    """Choose the visit with the closest date to a given image acquisition date.

    Args:
        image: An image entry from MPRAGEMETA
        pending_timepoints: List of visit entries from ADNIMERGE
        subject: Subject identifier
        visit_field:
        scandate_field:

    Returns:
        [Returns]
    """
    from datetime import datetime

    from clinica.utils.stream import cprint

    min_db = 100000
    min_db2 = 0
    min_visit = None
    min_visit2 = None

    for timepoint in pending_timepoints:
        db = days_between(image[scandate_field], timepoint.EXAMDATE)
        if db < min_db:
            min_db2 = min_db
            min_visit2 = min_visit

            min_db = db
            min_visit = timepoint

    if min_visit is None:
        cprint(
            f"No corresponding timepoint in ADNIMERGE for subject {subject} in visit {image[visit_field]}"
        )
        cprint(image, lvl="debug")
        return None

    if min_visit2 is not None and min_db > 90:
        cprint(
            msg=(
                f"More than 60 days for corresponding timepoint in ADNIMERGE for "
                f"subject {subject} in visit {image[visit_field]} on {image[scandate_field]}"
            ),
            lvl="debug",
        )
        cprint(
            msg=(
                f"Timepoint 1: {min_visit.VISCODE} - {min_visit.ORIGPROT} "
                f"on {min_visit.EXAMDATE} (Distance: {min_db} days)"
            ),
            lvl="debug",
        )
        cprint(
            msg=(
                f"Timepoint 2: {min_visit2.VISCODE} - {min_visit2.ORIGPROT} "
                f"on {min_visit2.EXAMDATE} (Distance: {min_db2} days)"
            ),
            lvl="debug",
        )

        # If image is too close to the date between two visits we prefer the earlier visit
        if (
            datetime.strptime(min_visit.EXAMDATE, "%Y-%m-%d")
            > datetime.strptime(image[scandate_field], "%Y-%m-%d")
            > datetime.strptime(min_visit2.EXAMDATE, "%Y-%m-%d")
        ):
            dif = days_between(min_visit.EXAMDATE, min_visit2.EXAMDATE)
            if abs((dif / 2.0) - min_db) < 30:
                min_visit = min_visit2

        cprint(msg=f"We prefer {min_visit.VISCODE}", lvl="debug")

    key_min_visit = (min_visit.VISCODE, min_visit.COLPROT, min_visit.ORIGPROT)

    return key_min_visit


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

    selected_image = None
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

            # There are no images that passed the qc
            # so we'll try to see if there are other images without qc,
            # otherwise return None
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

    diagnosis_dict = {1: "CN", 2: "MCI", 3: "AD"}
    dxsum_df = pd.read_csv(
        path.join(clinical_data_dir, "DXSUM_PDXCONV_ADNIALL.csv")
    ).set_index(["PTID", "VISCODE2"])
    missing_sc = participants_df[participants_df.original_study == "ADNI3"]
    participants_df.set_index("alternative_id_1", drop=True, inplace=True)
    for alternative_id in missing_sc.alternative_id_1.values:
        diagnosis_sc = diagnosis_dict[
            dxsum_df.loc[(alternative_id, "sc"), "DIAGNOSIS"].values[0]
        ]
        participants_df.loc[alternative_id, "diagnosis_sc"] = diagnosis_sc

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
    """Write the result of method create_session_dict into several TSV files.

    Args:
        df_subj_sessions: global dataframe containing clinical sessions data for all subjects
        bids_subjs_paths: a list with the path to all bids subjects
    """
    import os
    from os import path

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

    for sp in bids_subjs_paths:
        if not path.exists(sp):
            os.makedirs(sp)

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


def filter_subj_bids(df_files, location, bids_ids):
    import clinica.iotools.bids_utils as bids

    # Depending of the file that needs to be open, identify and
    # do needed preprocessing on the column that contains the
    # subjects ids
    bids_ids = [x[8:] for x in bids_ids if "sub-ADNI" in x]
    if location == "ADNIMERGE.csv":
        df_files["RID"] = df_files["PTID"].apply(
            lambda x: bids.remove_space_and_symbols(x)
        )
        df_files["RID"] = df_files["RID"].apply(lambda x: x[4:])
        df_ret = df_files[df_files["RID"].isin([x[4:] for x in bids_ids])]
    else:
        df_files["RID"] = df_files["RID"].apply(lambda x: pad_id(x))
        df_ret = df_files[df_files["RID"].isin([x[4:] for x in bids_ids])]
    return df_ret


def update_age(row):
    """Update age with time passed since bl to current visit"""
    from datetime import datetime

    if row["session_id"] != "ses-M00":
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


def get_visit_id(row, location):
    """Return a common visit ID across different files"""
    import pandas as pd

    locations_visicode2 = [
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
    ]

    if location in locations_visicode2:
        if pd.isnull(row["VISCODE2"]) or row["VISCODE2"] == "f":
            return False
        if row["VISCODE2"] == "sc":
            return "sc"  # visit_id = "bl"
        else:
            visit_id = row["VISCODE2"]
    elif location in [
        "BHR_EVERYDAY_COGNITION.csv",
        "BHR_BASELINE_QUESTIONNAIRE.csv",
        "BHR_LONGITUDINAL_QUESTIONNAIRE.csv",
    ]:
        visit_id = row["Timepoint"]
    else:
        visit_id = row["VISCODE"]
    return viscode_to_session(visit_id)


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

    from os import path

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
        if path.exists(path.join(clinical_data_dir, location)):

            file_to_read_path = path.join(clinical_data_dir, location)
            cprint(f"\tReading clinical data file: {location}")

            df_file = pd.read_csv(file_to_read_path, dtype=str)
            df_filtered = filter_subj_bids(df_file, location, bids_ids).copy()

            if not df_filtered.empty:
                df_filtered["session_id"] = df_filtered.apply(
                    lambda x: get_visit_id(x, location), axis=1
                )
                if location == "ADNIMERGE.csv":
                    df_filtered["AGE"] = df_filtered.apply(
                        lambda x: update_age(x), axis=1
                    )
                df_subj_session = update_sessions_df(
                    df_subj_session, df_filtered, df_sessions, location
                )
            else:
                dict_column_correspondance = dict(
                    zip(df_sessions["ADNI"], df_sessions["BIDS CLINICA"])
                )
                df_filtered.rename(columns=dict_column_correspondance, inplace=True)
                df_filtered = df_filtered.loc[
                    :, (~df_filtered.columns.isin(df_subj_session.columns))
                ]
                df_subj_session = pd.concat([df_subj_session, df_filtered], axis=1)

    df_subj_session.drop(
        df_subj_session[df_subj_session.session_id == "sc"].index, inplace=True
    )
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
    conversion_versions.sort()
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
            viscode = session_to_viscode(session_name[4::])
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
        tracer = file_name.split("acq-")[1].split("_")[0]
        if tracer == "av45" or tracer == "fbb":
            return "amyloid_pet"
        else:
            return f"{tracer}_pet"
    else:
        raise ValueError(f"Conversion modality could not be found for file {file_name}")


def find_image_path(images, source_dir, modality, prefix, id_field):
    """
    For each image, the path to an existing image file or folder is created from image metadata.

    Args:
        images: List of images metadata
        source_dir: path to the ADNI directory
        modality: Imaging modality
        prefix: Prefix to prepend to image identifier to create name of folder containing the image
        id_field: Name of the field in images metadata dataframe containing the image identifier to use

    Returns: Dataframe containing metadata and existing paths
    """
    from os import path, walk

    import pandas as pd

    from clinica.utils.stream import cprint

    is_dicom = []
    image_folders = []

    for row in images.iterrows():
        image = row[1]
        seq_path = path.join(source_dir, str(image.Subject_ID), image.Sequence)
        image_path = ""
        for (dirpath, dirnames, filenames) in walk(seq_path):
            found = False
            for d in dirnames:
                if d == prefix + str(image[id_field]):
                    image_path = path.join(dirpath, d)
                    found = True
                    break
            if found:
                break

        dicom = True
        for (dirpath, dirnames, filenames) in walk(image_path):
            for f in filenames:
                if f.endswith(".nii"):
                    dicom = False
                    image_path = path.join(dirpath, f)
                    break

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


def paths_to_bids(images, bids_dir, modality, mod_to_update=False):
    """Images in the list are converted and copied to directory in BIDS format.

    Args:
        images: List of images metadata and paths
        bids_dir: Path to the output BIDS directory
        modality: Imaging modality
        mod_to_update: If True, pre-existing images in the BIDS directory will be erased and extracted agai n.
    """
    from functools import partial
    from multiprocessing import Pool, Value, cpu_count

    if modality.lower() not in [
        "t1",
        "dwi",
        "flair",
        "fmri",
        "fdg",
        "pib",
        "av45_fbb",
        "tau",
    ]:
        # This should never be reached
        raise RuntimeError(
            modality.lower() + " is not supported for conversion in paths_to_bids"
        )

    counter = None

    def init(args):
        """Counter is stored for later use.

        Args:
            args:
        """
        global counter
        counter = args

    counter = Value("i", 0)
    total = images.shape[0]
    # Reshape inputs to give it as a list to the workers
    images_list = []
    for i in range(total):
        images_list.append(images.iloc[i])

    # Handle multiargument for create_file functions
    partial_create_file = partial(
        create_file,
        modality=modality,
        total=total,
        bids_dir=bids_dir,
        mod_to_update=mod_to_update,
    )

    # intializer are used with the counter variable to keep track of how many
    # files have been processed
    poolrunner = Pool(max(cpu_count() - 1, 1), initializer=init, initargs=(counter,))
    output_file_treated = poolrunner.map(partial_create_file, images_list)
    del counter
    return output_file_treated


def create_file(image, modality, total, bids_dir, mod_to_update):
    """
    Image file is created at the corresponding output folder
    as result of image conversion (DICOM to NIFTI) and centering,
    as specified for each modality

    Args:
        image: Image metadata
        modality: Imaging modality
        total: Total number of images to convert
        bids_dir: Path to the output BIDS directory
        mod_to_update: If True, pre-existing images in the BIDS directory will be erased and extracted again.

    Returns:
        [Returns]
    """
    import os
    import re
    import shutil
    import subprocess
    from glob import glob
    from os import path

    from numpy import nan

    from clinica.iotools.utils.data_handling import center_nifti_origin
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
            "output_filename": "_acq-axial_dwi",
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
            "output_filename": "_task-rest_acq-fdg_pet",
            "to_center": True,
            "json": "n",
        },
        "pib": {
            "output_path": "pet",
            "output_filename": "_task-rest_acq-pib_pet",
            "to_center": True,
            "json": "n",
        },
        "av45": {
            "output_path": "pet",
            "output_filename": "_task-rest_acq-av45_pet",
            "to_center": True,
            "json": "n",
        },
        "fbb": {
            "output_path": "pet",
            "output_filename": "_task-rest_acq-fbb_pet",
            "to_center": True,
            "json": "n",
        },
        "tau": {
            "output_path": "pet",
            "output_filename": "_task-rest_acq-tau_pet",
            "to_center": True,
            "json": "n",
        },
    }

    global counter

    subject = image.Subject_ID

    if modality == "av45_fbb":
        modality = image.Tracer.lower()

    with counter.get_lock():
        counter.value += 1
        if image.Path == "":
            cprint(
                msg=(
                    f"[{modality.upper()}] No path specified for {image.Subject_ID} "
                    f"in session {image.VISCODE} {counter.value}/{total}"
                ),
                lvl="info",
            )
        else:
            cprint(
                msg=(
                    f"[{modality.upper()}] Processing subject {subject} - session {image.VISCODE},"
                    f" {counter.value}/{total}"
                ),
                lvl="info",
            )

    if image.Path == "":
        return nan

    session = viscode_to_session(image.VISCODE)
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
            command = f"dcm2niix -b {generate_json} -z {zip_image} -o {output_path} -f {output_filename} {image_path}"
            subprocess.run(
                command,
                shell=True,
                stderr=subprocess.DEVNULL,
                stdout=subprocess.DEVNULL,
            )

            # If "_t" - the trigger delay time - exists in dcm2niix output filename, we remove it
            exception_t = glob(path.join(output_path, output_filename + "_t[0-9]*"))
            for trigger_time in exception_t:
                res = re.search("_t\d+\.", trigger_time)
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
                    msg=(
                        "Conversion with dcm2niix failed, trying with dcm2nii "
                        f"for subject {subject} and session {session}"
                    ),
                    lvl="warning",
                )
                command = f"dcm2nii -a n -d n -e n -i y -g {zip_image} -p n -m n -r n -x n -o {output_path} {image_path}"
                subprocess.run(
                    command,
                    shell=True,
                    stdout=subprocess.DEVNULL,
                    stderr=subprocess.DEVNULL,
                )

                # If modality is DWI we check if .bvec and .bval files are also present
                if modality == "dwi":
                    bvec_file = path.join(
                        output_path, subject.replace("_", "") + ".bvec"
                    )
                    bval_file = path.join(
                        output_path, subject.replace("_", "") + ".bval"
                    )

                    if os.path.exists(bvec_file) and os.path.exists(bval_file):
                        os.rename(
                            bvec_file, path.join(output_path, output_filename + ".bvec")
                        )
                        os.rename(
                            bval_file, path.join(output_path, output_filename + ".bval")
                        )
                    else:
                        cprint(
                            msg=(
                                "Expected bvec and bval not generated by dcm2nii "
                                f"for subject {subject} and session {session}"
                            ),
                            lvl="warning",
                        )

                nifti_files = glob(
                    path.join(output_path, subject.replace("_", "") + ".nii*")
                )
                output_image = path.join(output_path, output_filename + ".nii.gz")

                # If no NIFTI files were converted then conversion failed
                if nifti_files:
                    nifti_file = nifti_files[0]
                    conversion_failed = False
                else:
                    conversion_failed = True

                # If image is not going to be centered, then it is already a .nii.gz file
                # and we only need to rename it to the output image filename
                if (
                    not conversion_failed
                    and not modality_specific[modality]["to_center"]
                ):
                    os.rename(nifti_file, output_image)

                # We check if neither the NIFTI file to be centered
                # nor the final compressed NIFTI file exist. In this case, conversion failed
                if conversion_failed or (
                    not path.isfile(nifti_file) and not path.isfile(output_image)
                ):
                    cprint(
                        msg=(
                            "Conversion from DICOM failed "
                            f"for subject {subject} and session {session}. Image path: {image_path}"
                        ),
                        lvl="warning",
                    )
                    # If conversion failed we remove the folder if it is empty
                    if not os.listdir(output_path):
                        os.rmdir(output_path)
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
                if output_image is None:
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


def viscode_to_session(viscode):
    """Replace the session label 'bl' with 'M00' or capitalize the session name passed as input.

    Args:
        viscode: session name

    Returns:
        M00 if is the baseline session or the original session name capitalized
    """
    if viscode == "bl":
        return "ses-M00"
    else:
        return "ses-" + viscode.capitalize()


def session_to_viscode(session_name):
    """Replace the session label 'bl' with 'M00' or capitalize the session name passed as input.

    Args:
        session_name: MXX

    Returns:
        M00 if is the baseline session or the original session name capitalized
    """
    if session_name == "M00":
        return "bl"
    else:
        return session_name.lower()


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

        # Remove the precedent tmp_dcm_folder if is existing
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
