"""
Convert the AIBL dataset (https://www.aibl.csiro.au/) into BIDS.
"""


def convert_images(path_to_dataset, path_to_csv, bids_dir, overwrite=False):

    # Conversion of the entire dataset in BIDS
    from os.path import exists

    from clinica.iotools.converters.aibl_to_bids.aibl_utils import paths_to_bids
    from clinica.utils.stream import cprint

    list_of_created_files = []

    for modality in ["t1", "av45", "flute", "pib"]:
        list_of_created_files.append(
            paths_to_bids(
                path_to_dataset, path_to_csv, bids_dir, modality, overwrite=overwrite
            )
        )

    error_string = ""
    for modality_list in list_of_created_files:
        for file in modality_list:
            if not exists(str(file)):
                error_string = error_string + str(file) + "\n"
    if error_string != "":
        cprint(
            msg=f"The following file were not converted: {error_string}", lvl="warning"
        )


def convert_clinical_data(bids_dir, path_to_csv):
    # clinical specifications in BIDS
    from os.path import join, realpath, split

    import clinica.iotools.bids_utils as bids
    from clinica.iotools.converters.aibl_to_bids.aibl_utils import (
        create_participants_df_AIBL,
        create_scans_dict_AIBL,
        create_sessions_dict_AIBL,
    )
    from clinica.utils.stream import cprint

    clinical_spec_path = join(
        split(realpath(__file__))[0], "../../data/clinical_specifications"
    )
    # if not exists(clinical_spec_path):
    #    raise FileNotFoundError(
    #        f"{clinical_spec_path} file not found ! This is an internal file of Clinica."
    #    )

    cprint("Creating modality agnostic files...")
    readme_data = {
        "link": "http://adni.loni.usc.edu/study-design/collaborative-studies/aibl/",
        "desc": (
            "The Australian Imaging, Biomarker & Lifestyle Flagship Study of Ageing (AIBL) seeks to discover which "
            "biomarkers, cognitive characteristics, and health and lifestyle factors determine the development of AD. "
            "Although AIBL and ADNI have many of the same goals, there are differences between the two projects."
        ),
    }
    bids.write_modality_agnostic_files(
        study_name="AIBL", readme_data=readme_data, bids_dir=bids_dir
    )

    cprint("Creating participants.tsv...")
    create_participants_df_AIBL(
        bids_dir, clinical_spec_path, path_to_csv, delete_non_bids_info=True
    )

    cprint("Creating sessions files...")
    create_sessions_dict_AIBL(bids_dir, path_to_csv, clinical_spec_path)

    cprint("Creating scans files...")
    create_scans_dict_AIBL(bids_dir, path_to_csv, clinical_spec_path)
