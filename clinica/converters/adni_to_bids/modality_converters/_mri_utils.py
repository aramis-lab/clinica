from pathlib import Path

import pandas as pd

__all__ = ["load_all_images_metadata"]


def load_all_images_metadata(csv_dir: Path) -> pd.DataFrame:
    """All Images and MRIQC csv files are loaded. Columns from MRIQC file are merged into All Images.

    Parameters
    ----------
    csv_dir : PathLike
        Path to the clinical data directory.

    Returns: DataFrame containing main clinical data for PET images.
    """
    from clinica.converters._utils import load_clinical_csv, load_manifest_metadata

    all_images = load_clinical_csv(csv_dir, "All_Images")
    mriqc = load_clinical_csv(csv_dir, "MRIQC")
    adnimerge = load_clinical_csv(csv_dir, "ADNIMERGE")
    manifest = load_manifest_metadata(csv_dir)

    # 1. Add the series identifier from Manifest.
    all_images = all_images.merge(
        manifest[["image_id", "series_id"]].drop_duplicates(),
        on="image_id",
        how="left",
        validate="many_to_one",
    )

    # 2. Get one unambiguous field strength per subject and study.
    mriqc = mriqc.rename(
        columns={"ParticipantID": "subject_id", "LONIStudy": "study_id"}
    )

    study_fields = (
        mriqc.groupby(["subject_id", "study_id"])["MagneticFieldStrength"]
        .agg(
            lambda values: values.dropna().iloc[0]
            if values.nunique() == 1
            else float("nan")
        )
        .reset_index()
    )

    all_images = all_images.merge(
        study_fields,
        on=["subject_id", "study_id"],
        how="left",
        validate="many_to_one",
    )

    # 3. Add ADNIMERGE's visit-level field strength.
    visit_column = "VISCODE2" if "VISCODE2" in adnimerge.columns else "VISCODE"

    adnimerge = adnimerge.rename(
        columns={"PTID": "subject_id", visit_column: "image_visit"}
    ).copy()

    # Treat these visit codes as equivalent for the fallback merge.
    baseline_codes = {"bl": "bl", "sc": "bl", "m0": "bl", "m00": "bl"}

    for df in (all_images, adnimerge):
        df["visit_key"] = (
            df["image_visit"]
            .astype("string")
            .str.strip()
            .str.lower()
            .replace(baseline_codes)
        )

    adnimerge["FLDSTRENG"] = pd.to_numeric(
        adnimerge["FLDSTRENG"]
        .astype("string")
        .str.extract(r"(\d+(?:\.\d+)?)", expand=False),
        errors="coerce",
    )

    # Collapsing visits can produce multiple values: keep only agreement.
    visit_fields = (
        adnimerge.groupby(["subject_id", "visit_key"])["FLDSTRENG"]
        .agg(
            lambda values: values.dropna().iloc[0]
            if values.nunique() == 1
            else float("nan")
        )
        .reset_index()
    )

    all_images = all_images.merge(
        visit_fields,
        on=["subject_id", "visit_key"],
        how="left",
        validate="many_to_one",
    ).drop(columns="visit_key")

    fallback = all_images["FLDSTRENG"]

    # Record provenance, then fill only missing MRIQC values.
    all_images["field_strength_source"] = pd.NA
    all_images.loc[
        all_images["MagneticFieldStrength"].notna(), "field_strength_source"
    ] = "MRIQC"

    use_fallback = all_images["MagneticFieldStrength"].isna() & fallback.notna()
    all_images.loc[use_fallback, "MagneticFieldStrength"] = fallback[use_fallback]
    all_images.loc[use_fallback, "field_strength_source"] = "ADNIMERGE"

    # Field strength applies only to MR modalities.
    is_mr = all_images["modality"].isin(["MRI", "fMRI", "DTI"])
    all_images.loc[~is_mr, "MagneticFieldStrength"] = float("nan")
    all_images.loc[~is_mr, "field_strength_source"] = pd.NA

    all_images = all_images.drop(columns="FLDSTRENG")

    all_images = all_images.rename(columns={"LONISeries": "series_id"})

    return all_images
