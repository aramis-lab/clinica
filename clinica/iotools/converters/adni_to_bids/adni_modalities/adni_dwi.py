# coding: utf-8

"""Module for converting DWI of ADNI."""


def convert_adni_dwi(
    source_dir, csv_dir, dest_dir, conversion_dir, subjs_list=None, mod_to_update=False
):
    """Convert DW images of ADNI into BIDS format.

    Args:
        source_dir: path to the ADNI directory
        csv_dir: path to the clinical data directory
        dest_dir: path to the destination BIDS directory
        conversion_dir: path to the TSV files including the paths to original images
        subjs_list: subjects list
        mod_to_update: If True, pre-existing images in the BIDS directory will be erased and extracted again.
    """
    from os import path

    import pandas as pd

    from clinica.iotools.converters.adni_to_bids.adni_utils import paths_to_bids
    from clinica.utils.stream import cprint

    if not subjs_list:
        adni_merge_path = path.join(csv_dir, "ADNIMERGE.csv")
        adni_merge = pd.read_csv(adni_merge_path, sep=",", low_memory=False)
        subjs_list = list(adni_merge.PTID.unique())

    cprint(
        f"Calculating paths of DWI images. Output will be stored in {conversion_dir}."
    )
    images = compute_dwi_paths(source_dir, csv_dir, subjs_list, conversion_dir)
    cprint("Paths of DWI images found. Exporting images into BIDS ...")
    # dwi_paths_to_bids(images, dest_dir)
    paths_to_bids(images, dest_dir, "dwi", mod_to_update=mod_to_update)
    cprint(msg="DWI conversion done.", lvl="debug")


def compute_dwi_paths(source_dir, csv_dir, subjs_list, conversion_dir):
    """Compute paths to DW images to convert to BIDS.

    Args:
        source_dir: path to the ADNI directory
        csv_dir: path to the clinical data directory
        subjs_list: subjects list
        conversion_dir: path to the TSV files including the paths to original images

    Returns:
        images: pandas dataframe that contains the path to all the DW images to convert
    """
    from os import path

    import pandas as pd

    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        find_image_path,
        visits_to_timepoints,
    )

    dwi_col_df = [
        "Subject_ID",
        "VISCODE",
        "Visit",
        "Sequence",
        "Scan_Date",
        "Study_ID",
        "Series_ID",
        "Image_ID",
        "Field_Strength",
    ]
    dwi_df = pd.DataFrame(columns=dwi_col_df)
    dwi_dfs_list = []

    # Loading needed .csv files
    adni_merge = pd.read_csv(
        path.join(csv_dir, "ADNIMERGE.csv"), sep=",", low_memory=False
    )

    mayo_mri_qc = pd.read_csv(
        path.join(csv_dir, "MAYOADIRL_MRI_IMAGEQC_12_08_15.csv"),
        sep=",",
        low_memory=False,
    )
    mayo_mri_qc = mayo_mri_qc[mayo_mri_qc.series_type == "DTI"]

    mri_list = pd.read_csv(path.join(csv_dir, "MRILIST.csv"), sep=",", low_memory=False)

    # Selecting only DTI images that are not Multiband, processed or enhanced images
    mri_list = mri_list[mri_list.SEQUENCE.str.contains("dti", case=False, na=False)]
    unwanted_sequences = ["MB", "ADC", "FA", "TRACEW", "Enhanced", "Reg"]
    mri_list = mri_list[
        mri_list.SEQUENCE.map(
            lambda x: not any(subs in x for subs in unwanted_sequences)
        )
    ]

    for subj in subjs_list:

        # Filter ADNIMERGE, MRI_LIST and QC for only one subject and sort the rows/visits by examination date
        adnimerge_subj = adni_merge[adni_merge.PTID == subj]
        adnimerge_subj = adnimerge_subj.sort_values("EXAMDATE")

        mri_list_subj = mri_list[mri_list.SUBJECT == subj]
        mri_list_subj = mri_list_subj.sort_values("SCANDATE")

        mayo_mri_qc_subj = mayo_mri_qc[mayo_mri_qc.RID == int(subj[-4:])]

        # Obtain corresponding timepoints for the subject visits
        visits = visits_to_timepoints(subj, mri_list_subj, adnimerge_subj, "DWI")

        for visit_info in visits.keys():
            timepoint = visit_info[0]
            visit_str = visits[visit_info]

            visit_mri_list = mri_list_subj[mri_list_subj.VISIT == visit_str]
            axial = dwi_image(
                subj, timepoint, visits[visit_info], visit_mri_list, mayo_mri_qc_subj
            )

            if axial is not None:
                row_to_append = pd.DataFrame(
                    axial,
                    index=[
                        "i",
                    ],
                )
                dwi_dfs_list.append(row_to_append)

    if dwi_dfs_list:
        dwi_df = pd.concat(dwi_dfs_list, ignore_index=True)

    # Exceptions
    # ==========
    conversion_errors = [
        ("029_S_2395", "m60"),
        ("029_S_0824", "m108"),
        ("029_S_0914", "m108"),
        ("027_S_2219", "m36"),
        ("129_S_2332", "m12"),
        ("029_S_4384", "m48"),
        ("029_S_4385", "m48"),
        ("029_S_4585", "m48"),
        ("016_S_4591", "m24"),
        ("094_S_4630", "m06"),
        ("094_S_4649", "m06"),
        ("029_S_5219", "m24"),
        ("094_S_2238", "m48"),
        ("129_S_4287", "bl"),
        ("007_S_4611", "m03"),
        ("016_S_4638", "bl"),
        ("027_S_5118", "bl"),
        ("098_S_4018", "bl"),
        ("098_S_4003", "m12"),
        ("016_S_4584", "m24"),
        ("016_S_5007", "m12"),
        ("129_S_2347", "m06"),
        ("129_S_4220", "bl"),
        ("007_S_2058", "m12"),
        ("016_S_2007", "m06"),
        ("020_S_6358", "bl"),
        ("114_S_6039", "m12"),
        ("114_S_6057", "bl"),
        ("153_S_6274", "bl"),
        ("006_S_4485", "m84"),
        ("153_S_6237", "bl"),
        ("153_S_6336", "bl"),
        ("153_S_6450", "bl"),
        ("003_S_4441", "m12"),
        # Eq_1 Images
        ("006_S_6252", "m12"),
        ("010_S_0419", "m156"),
        ("127_S_0259", "m156"),
        ("129_S_6830", "bl"),
        ("135_S_6104", "m24"),
        ("021_S_5237", "m84"),
        ("027_S_2245", "m120"),
        ("099_S_6396", "m24"),
        ("126_S_4507", "m96"),
        ("126_S_4891", "m96"),
        ("126_S_4896", "m96"),
        ("126_S_6559", "m24"),
        ("126_S_6724", "m12"),
        ("127_S_6203", "m24"),
        ("127_S_6330", "m24"),
        ("127_S_6512", "m24"),
        ("127_S_6549", "m24"),
        ("129_S_6459", "m24"),
        ("129_S_6763", "m12"),
        ("129_S_6784", "m12"),
        ("129_S_6830", "m12"),
        ("135_S_6446", "m24"),
        ("301_S_6224", "m36"),
        # Missing .bval and .bvec
        ("020_S_5203", "m72"),
        ("020_S_6185", "m24"),
        ("020_S_6513", "m12"),
        ("021_S_0178", "m156"),
        ("153_S_6336", "m12"),
        ("153_S_6755", "bl"),
        ("020_S_6227", "m24"),
        ("020_S_6449", "m24"),
        ("020_S_6513", "m24"),
        # Volume mismatch between .nii and .bvec / .bval
        ("006_S_6610", "bl"),
        ("006_S_6682", "bl"),
        ("006_S_6696", "bl"),
        ("006_S_6770", "bl"),
        ("029_S_6289", "bl"),
        ("130_S_6043", "bl"),
        ("130_S_6329", "bl"),
        ("027_S_6183", "m24"),
        ("123_S_6891", "bl"),
        # Several output images
        ("029_S_2395", "m72"),
    ]

    # Removing known exceptions from images to convert
    if not dwi_df.empty:
        error_ind = dwi_df.index[
            dwi_df.apply(
                lambda x: ((x.Subject_ID, x.VISCODE) in conversion_errors), axis=1
            )
        ]
        dwi_df.drop(error_ind, inplace=True)

    # Checking for images paths in filesystem
    images = find_image_path(dwi_df, source_dir, "DWI", "S", "Series_ID")
    images.to_csv(path.join(conversion_dir, "dwi_paths.tsv"), sep="\t", index=False)

    return images


def dwi_image(subject_id, timepoint, visit_str, visit_mri_list, mri_qc_subj):
    """
    One image among those in the input list is chosen according to QC
    and then correspoding metadata is extracted to a dictionary.

    Args:
        subject_id: Subject identifier
        timepoint: Visit code
        visit_str: Visit name
        visit_mri_list: List of images metadata
        mri_qc_subj: Dataframe containing list of QC of scans for the subject

    Returns: dictionary - contains image metadata
    """
    from clinica.iotools.converters.adni_to_bids.adni_utils import (
        replace_sequence_chars,
        select_image_qc,
    )

    sel_image = select_image_qc(list(visit_mri_list.IMAGEUID), mri_qc_subj)
    if not sel_image:
        return None

    sel_scan = visit_mri_list[visit_mri_list.IMAGEUID == sel_image].iloc[0]

    image_dict = {
        "Subject_ID": subject_id,
        "VISCODE": timepoint,
        "Visit": visit_str,
        "Sequence": replace_sequence_chars(sel_scan.SEQUENCE),
        "Scan_Date": sel_scan["SCANDATE"],
        "Study_ID": str(int(sel_scan.STUDYID)),
        "Series_ID": str(int(sel_scan.SERIESID)),
        "Image_ID": str(int(sel_scan.IMAGEUID)),
        "Field_Strength": sel_scan.MAGSTRENGTH,
    }

    return image_dict
