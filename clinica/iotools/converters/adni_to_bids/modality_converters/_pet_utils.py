import pandas as pd

from clinica.utils.stream import cprint

__all__ = ["get_images_pet"]


def get_images_pet(
    subject: str,
    pet_qc_subj: pd.DataFrame,
    subject_pet_meta: pd.DataFrame,
    df_cols: list[str],
    modality: str,
    sequences_preprocessing_step: list[str],
    viscode_field: str = "VISCODE2",
) -> list[pd.DataFrame]:
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
    from clinica.iotools.converter_utils import replace_sequence_chars
    from clinica.utils.pet import Tracer

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
                f"No {modality} images metadata for subject {subject} and visit {qc_visit[viscode_field]}",
                lvl="info",
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
                tracer = Tracer.AV45.value
            elif "fbb" in sel_image.Sequence.lower():
                tracer = Tracer.FBB.value
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
